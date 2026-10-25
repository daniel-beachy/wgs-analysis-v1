"""Ancestry (ADR-020): which reference populations your DNA is most like, and your maternal and paternal lineages.

Genetic ancestry here means similarity to present-day reference people, not where ancestors lived or nationality.
Every number is calibrated on the reference panel itself, so the page can say how often the method is wrong.

1. **Continental shares.** Your position on the ancestry PCA (computed by the PGS stage with FRAPOSA on the 1000
   Genomes panel) is written as the closest mix of the five super-population centres (non-negative weights that sum
   to one), the projection approach used by rye (Conley et al. 2023). Each reference person is scored the same way
   with themselves left out; the share the method wrongly gives a group among reference people who are not from it
   (95th percentile) is that group's noise floor, and smaller shares are reported as "not distinguishable from zero".
2. **Nearest reference population** within your main group: the 10 nearest reference people. Leave-one-out shows
   which populations the method can tell apart (Finnish vs Tuscan: yes; British vs CEU: no), and populations it
   confuses with the nearest one, directly or through a chain of confusions, are reported together. When the
   high-depth panels are in the knowledge store (``ancestry_panels``), the comparison adds the Human Genome Diversity
   Project's populations: a PCA of everyone in your main group (1000 Genomes 30× + HGDP + you, so you are not
   projected) on about 40,000 positions, with the number of axes the leave-one-out test does best with. Otherwise it
   uses the 1000 Genomes PCA above.
3. **Maternal line** (mitochondrial DNA): haplogrep3 on your chrM calls (PhyloTree 17 Forensic Update).
4. **Paternal line** (Y chromosome, XY only): Yleaf on your Y reads (YFull tree), or on the VCF when no alignment is
   available (lower confidence: positions absent from a VCF are assumed to match the reference).
"""

from __future__ import annotations

import itertools
import json
import shutil
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from ..knowledge import Store
from ..knowledge.ancestry_panels import MISSING, read_bed
from ..pipeline import Context, Stage, run, write_json
from . import pgs, qc, reference, variants
from .pgx import _env_tool

SOURCES = ["AFR", "EUR", "EAS", "SAS", "AMR"]
ADMIXED = {"ASW", "ACB", "MXL", "PUR", "CLM", "PEL"}  # IGSR: recently admixed groups, not used to set noise floors
N_PCS = 10
K_NEIGHBOURS = 10
NOISE_QUANTILE = 0.95
CONFUSED = 0.20          # two populations are "not told apart" if ≥ 20% of either is classified as the other
MAIN_SHARE = 0.80        # below this, nearest-population within one group is not meaningful
MIN_PCA_VARIANTS = 50_000
REGION_PCS = (4, 6, 8, 10)   # axes tried for the finer comparison; the one with the best leave-one-out accuracy wins
MIN_REGION_SITES = 20_000    # positions in the finer PCA for a High call
# The same population sampled by both panels. If a pair is told apart, the panels differ for technical reasons
# (as phase 3's low-depth genotypes did) and the finer comparison can't be trusted (self-check).
TWINS = {"TuscanHGDP": "TSI", "JapaneseHGDP": "JPT", "DaiHGDP": "CDX", "YorubaHGDP": "YRI", "NorthernHanHGDP": "CHB"}
MT_HIGH, MT_MEDIUM = 0.90, 0.80   # haplogrep3 quality, 0–1 (thresholds are this project's, ADR-020)
Y_HIGH, Y_MEDIUM = 0.95, 0.90     # Yleaf QC-score (its default prediction threshold is 0.95)
MT_TREE = "phylotree-fu-rcrs@1.2"
YFULL = "https://www.yfull.com/tree/{hg}/"
YFULL_MT = "https://www.yfull.com/mtree/{hg}/"
SUPERPOP_NAME = {"AFR": "African", "EUR": "European", "EAS": "East Asian", "SAS": "South Asian",
                 "AMR": "Admixed American"}


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "ancestry"
    p.mkdir(parents=True, exist_ok=True)
    return p


def out_components(ctx: Context) -> Path:
    return out_dir(ctx) / "components.parquet"


def out_populations(ctx: Context) -> Path:
    return out_dir(ctx) / "populations.parquet"


def out_haplogroups(ctx: Context) -> Path:
    return out_dir(ctx) / "haplogroups.parquet"


def out_pca(ctx: Context) -> Path:
    return out_dir(ctx) / "pca.parquet"


def out_region_pca(ctx: Context) -> Path:
    return out_dir(ctx) / "region_pca.parquet"


def out_summary(ctx: Context) -> Path:
    return out_dir(ctx) / "ancestry.json"


# ---------------------------------------------------------------------------------------------------------------
# Continental shares
# ---------------------------------------------------------------------------------------------------------------


def mix(centres: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Weights w ≥ 0, Σw = 1 minimising |centresᵀw − x|². Exact: with ≤ 5 sources every support set is tried."""
    k = centres.shape[0]
    best, best_w = np.inf, np.full(k, np.nan)
    for r in range(1, k + 1):
        for s in itertools.combinations(range(k), r):
            a = centres[list(s)].T
            m = np.block([[2 * a.T @ a, np.ones((r, 1))], [np.ones((1, r)), np.zeros((1, 1))]])
            b = np.concatenate([2 * a.T @ x, [1.0]])
            try:
                sol = np.linalg.solve(m, b)[:r]
            except np.linalg.LinAlgError:
                continue
            if (sol < -1e-9).any():
                continue
            res = float(np.sum((a @ sol - x) ** 2))
            if res < best - 1e-12:
                best, best_w = res, np.zeros(k)
                best_w[list(s)] = np.clip(sol, 0, None)
    return best_w / best_w.sum()


def components(ref: dict, me: np.ndarray) -> tuple[list[dict], dict]:
    """Your shares, plus leave-one-out noise floors and accuracy on the reference panel."""
    x, sp, pop = ref["pcs"], ref["superpop"], ref["population"]
    centres = np.stack([x[sp == s].mean(0) for s in SOURCES])
    you = mix(centres, me)
    loo = np.zeros((len(x), len(SOURCES)))
    for i in range(len(x)):
        keep = np.arange(len(x)) != i
        c = np.stack([x[keep & (sp == s)].mean(0) for s in SOURCES])
        loo[i] = mix(c, x[i])
    pure = ~np.isin(pop, list(ADMIXED))
    floors, own = {}, {}
    for j, s in enumerate(SOURCES):
        others = pure & (sp != s)
        floors[s] = float(np.quantile(loo[others, j], NOISE_QUANTILE))
        mine = pure & (sp == s)
        own[s] = float(loo[mine, j].mean()) if mine.any() else None
    rows = [{"superpop": s, "share": float(you[j]), "noise_floor": floors[s],
             "detected": bool(you[j] > floors[s]), "reference_own_share": own[s]} for j, s in enumerate(SOURCES)]
    return rows, {"loo_own_share": own, "noise_floors": floors}


# ---------------------------------------------------------------------------------------------------------------
# Nearest reference population
# ---------------------------------------------------------------------------------------------------------------


def _knn(xs: np.ndarray, labels: np.ndarray, q: np.ndarray, k: int, skip: int | None = None) -> dict[str, int]:
    d = np.linalg.norm(xs - q, axis=1)
    if skip is not None:
        d[skip] = np.inf
    idx = np.argsort(d, kind="stable")[:k]
    out: dict[str, int] = {}
    for lab in labels[idx]:
        out[str(lab)] = out.get(str(lab), 0) + 1
    return out


def _vote(counts: dict[str, int]) -> str:
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def zoom_pcs(ref: dict, superpop: str) -> list[int]:
    """The two PCs that best separate the populations inside one group (largest between/within variance ratio)."""
    sel = ref["superpop"] == superpop
    x, y = ref["pcs"][sel], ref["population"][sel]
    pops = sorted(set(y.tolist()))
    within = sum(((x[y == p] - x[y == p].mean(0)) ** 2).sum(0) for p in pops)
    between = sum((y == p).sum() * (x[y == p].mean(0) - x.mean(0)) ** 2 for p in pops)
    ratio = between / np.where(within == 0, 1, within)
    return sorted(int(i) + 1 for i in np.argsort(-ratio)[:2])


def confusable(conf: dict[str, dict[str, float]], best: str) -> set[str]:
    """Populations that can't be told apart from `best`, directly or through a chain: if CHB is mixed up with Han
    and with CHS, a "Han" result can't rule out CHS either."""
    group, todo = {best}, [best]
    while todo:
        a = todo.pop()
        for b in conf:
            if b not in group and (conf[a][b] >= CONFUSED or conf[b][a] >= CONFUSED):
                group.add(b)
                todo.append(b)
    return group


def nearest(ref: dict, me: np.ndarray, superpop: str) -> tuple[list[dict], dict]:
    """kNN within one super-population, with a leave-one-out confusion matrix."""
    sel = ref["superpop"] == superpop
    x, y = ref["pcs"][sel], ref["population"][sel]
    sd = x.std(0)
    sd[sd == 0] = 1
    xs, q = x / sd, me / sd
    pops = sorted(set(y.tolist()))
    pred = np.array([_vote(_knn(xs, y, xs[i], K_NEIGHBOURS, skip=i)) for i in range(len(xs))])
    conf = {a: {b: float(np.mean(pred[y == a] == b)) for b in pops} for a in pops}
    votes = _knn(xs, y, q, K_NEIGHBOURS)
    best = _vote(votes)
    group = sorted(confusable(conf, best), key=lambda p: (-votes.get(p, 0), p))
    rows = [{"code": p, "neighbours": votes.get(p, 0), "loo_recall": conf[p][p],
             "confused_with": json.dumps({b: round(v, 3) for b, v in conf[p].items() if b != p and v > 0}),
             "nearest": p == best, "in_nearest_group": p in group} for p in pops]
    return rows, {"nearest": best, "group": group, "k": K_NEIGHBOURS, "accuracy": float(np.mean(pred == y)),
                  "confusion": conf}


def region_pca(x: np.ndarray, you: np.ndarray, n: int = max(REGION_PCS)) -> tuple[np.ndarray, np.ndarray, int]:
    """PCA of reference people (rows of x: copies of each position's alt allele, NaN = not called) plus you as one
    more row, so your point is not shrunk the way a projected sample is. Positions you weren't called at, or
    missing in more than 5% of people, or rarer than 1%, are left out. Returns (reference PCs, your PCs, positions)."""
    a = np.vstack([x, you])
    af = np.nanmean(a, 0) / 2
    keep = (np.isnan(a).mean(0) < 0.05) & (np.minimum(af, 1 - af) > 0.01) & ~np.isnan(you)
    a = a[:, keep]
    mu = np.nanmean(a, 0)
    a = np.where(np.isnan(a), mu, a)
    z = (a - mu) / np.sqrt(mu * (1 - mu / 2))  # sd of a 0/1/2 genotype under Hardy-Weinberg: sqrt(2p(1-p))
    u, sv, _ = np.linalg.svd(z, full_matrices=False)
    pcs = u[:, :n] * sv[:n]
    return pcs[:-1], pcs[-1], int(keep.sum())


def twins(conf: dict[str, dict[str, float]]) -> dict[str, dict]:
    """For each population sampled by both panels: was it confused with its counterpart, as it should be?"""
    return {h: {"twin": t, "confused": conf[h][t] >= CONFUSED or conf[t][h] >= CONFUSED}
            for h, t in TWINS.items() if h in conf and t in conf}


def regional(x: np.ndarray, labels: np.ndarray, you: np.ndarray, superpop: str) -> dict:
    """The finer nearest-population comparison: PCA within one continental group, then the kNN vote on the number of
    axes with the best leave-one-out accuracy."""
    ref_pcs, my_pcs, n_sites = region_pca(x, you)
    best = None
    for k in REGION_PCS:
        ref = {"pcs": ref_pcs[:, :k], "superpop": np.full(len(labels), superpop), "population": labels}
        rows, near = nearest(ref, my_pcs[:k], superpop)
        if best is None or near["accuracy"] > best[1]["accuracy"] + 1e-9:
            best = (rows, near, k)
    rows, near, k = best
    return {"rows": rows, "near": {**near, "pcs": k, "positions": n_sites, "twins": twins(near["confusion"])},
            "ref_pcs": ref_pcs, "my_pcs": my_pcs}


def _panels(ctx: Context, superpop: str) -> tuple[np.ndarray, list[dict], np.ndarray] | None:
    """The high-depth reference people of one continental group (genotypes, sample rows) and your genotypes at the
    same positions, all as copies of the GRCh37 alt allele."""
    st = Store(ctx.cfg)
    gt, sm = st.table("ancestry_panels", "genotypes"), st.table("ancestry_panels", "samples")
    stem = pgs.pgs_dir(ctx) / "pca" / "self"
    if not gt or not sm or not stem.with_suffix(".bed").exists():
        return None
    bim, _, mine = read_bed(stem)
    me: dict[str, float] = {}
    for row, g in zip(bim, mine[:, 0], strict=True):
        alt = row[1].split(":")[-1]
        if g != MISSING and alt in (row[4], row[5]):
            me[row[1]] = float(g if row[4] == alt else 2 - g)
    people = [r for r in pq.read_table(sm).to_pylist() if r["superpop"] == superpop and not r["related"]]
    if not people:
        return None
    t = pq.read_table(gt, columns=["id", "genotypes"])
    ids = t.column("id").to_pylist()
    at = [i for i, sid in enumerate(ids) if sid in me]
    blobs = t.column("genotypes")
    g = np.frombuffer(b"".join(blobs[i].as_py() for i in at), dtype=np.uint8).reshape(len(at), -1)
    x = g[:, [r["col"] for r in people]].T.astype(float)
    x[x == MISSING] = np.nan
    return x, people, np.array([me[ids[i]] for i in at])


def code(person: dict) -> str:
    """Population code as IGSR names it: 1000 Genomes codes (GBR), HGDP names plus "HGDP" (YorubaHGDP)."""
    return person["population"] + ("HGDP" if person["panel"] == "hgdp" else "")


# ---------------------------------------------------------------------------------------------------------------
# Haplogroups
# ---------------------------------------------------------------------------------------------------------------


def _level(q: float | None, high: float, medium: float) -> str:
    return "Not callable" if q is None else "High" if q >= high else "Medium" if q >= medium else "Low"


def _vcf_contig(ctx: Context, names: tuple[str, ...]) -> str | None:
    res = run([shutil.which("bcftools") or "bcftools", "view", "-h", variants.vcf_with_index(ctx)], capture=True)
    have = {ln.split("ID=")[1].split(",")[0] for ln in res.stdout.splitlines() if ln.startswith("##contig=<ID=")}
    return next((n for n in names if n in have), None)


def maternal(ctx: Context, d: Path) -> dict | None:
    hg3 = _env_tool("haplo", "haplogrep3")
    contig = _vcf_contig(ctx, ("chrM", "MT", "M", "chrMT"))
    if not hg3 or not contig:
        return None
    vcf, out = d / "mt.vcf.gz", d / "haplogrep.txt"
    run([shutil.which("bcftools") or "bcftools", "view", "-r", contig, "-f", "PASS", variants.vcf_with_index(ctx),
         "-Oz", "-o", str(vcf)], log=ctx.logs / "ancestry-mt.log")
    java = Path(hg3).parent.parent / "lib" / "jvm" / "bin"
    run([hg3, "classify", "--input", str(vcf), "--output", str(out), "--tree", MT_TREE, "--extend-report",
         "--hits", "3"], log=ctx.logs / "ancestry-mt.log", env={"PATH": f"{java}:{Path(hg3).parent}:/usr/bin:/bin"})
    hits = duckdb.sql(f"SELECT * FROM read_csv('{out}', delim='\t', header=true, all_varchar=true)").fetchdf()
    if hits.empty:
        return None
    hits = hits.fillna("")  # haplogrep writes empty fields, which pandas reads as NaN
    top = hits.iloc[0]
    quality = float(top["Quality"])
    split = lambda s: [t for t in str(s).split() if t and not t.startswith("(")]  # noqa: E731
    return {
        "line": "maternal", "haplogroup": top["Haplogroup"], "quality": quality,
        "call_level": _level(quality, MT_HIGH, MT_MEDIUM),
        "call_reasons": [f"haplogrep3 quality {quality:.2f} (High ≥ {MT_HIGH}, Medium ≥ {MT_MEDIUM})",
                         f"{len(split(top['Found_Polys']))} defining mutations found, "
                         f"{len(split(top['Not_Found_Polys']))} expected but missing"],
        "found": split(top["Found_Polys"]), "missing": split(top["Not_Found_Polys"]),
        "extra": [t for t in str(top["Remaining_Polys"]).split(") ") if t],
        "alternatives": [{"haplogroup": r["Haplogroup"], "quality": float(r["Quality"])}
                         for _, r in hits.iloc[1:].iterrows()],
        "path": None, "tool": "haplogrep3", "tree": "PhyloTree 17 – Forensic Update 1.2 (rCRS)", "input": "VCF",
        "url": YFULL_MT.format(hg=top["Haplogroup"]),
    }


def _y_bam(ctx: Context, d: Path) -> Path | None:
    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    fasta = reference.reference_fasta(ctx)
    if not aln or (aln.kind == "cram" and not fasta):
        return None
    contig = next((c for c in ("chrY", "Y") if c in reference.cram_contigs(ctx)), None)
    if not contig:
        return None
    bam = ctx.scratch / "ancestry-chrY.bam"
    samtools = shutil.which("samtools") or "samtools"
    ref = ["-T", str(fasta)] if aln.kind == "cram" else []
    run([samtools, "view", "-@", str(min(ctx.threads, 4)), "-b", *ref, "-o", str(bam), str(aln.path), contig],
        log=ctx.logs / "ancestry-y.log")
    run([samtools, "index", str(bam)], log=ctx.logs / "ancestry-y.log")
    return bam


def paternal(ctx: Context, d: Path) -> dict | None:
    yleaf = _env_tool("haplo", "Yleaf")
    fasta = reference.reference_fasta(ctx)
    if not yleaf or not fasta:
        return None
    out = d / "yleaf"
    shutil.rmtree(out, ignore_errors=True)
    bam = _y_bam(ctx, d)
    if bam:
        args, source = ["-bam", str(bam)], "aligned reads (CRAM/BAM)"
    else:
        contig = _vcf_contig(ctx, ("chrY", "Y"))
        if not contig:
            return None
        vy = ctx.scratch / "ancestry-chrY.vcf.gz"
        run([shutil.which("bcftools") or "bcftools", "view", "-r", contig, variants.vcf_with_index(ctx), "-Oz",
             "-o", str(vy)], log=ctx.logs / "ancestry-y.log")
        args, source = ["-vcf", str(vy)], "VCF (positions not in the VCF assumed to match the reference)"
    env = {"PATH": f"{Path(yleaf).parent}:/usr/bin:/bin"}
    run([yleaf, *args, "--ref-fasta", str(fasta), "-rg", "hg19", "-o", str(out), "-t", str(min(ctx.threads, 4)),
         "-force"], log=ctx.logs / "ancestry-y.log", env=env, cwd=ctx.scratch)
    for p in (ctx.scratch / "ancestry-chrY.bam", ctx.scratch / "ancestry-chrY.bam.bai",
              ctx.scratch / "ancestry-chrY.vcf.gz"):
        p.unlink(missing_ok=True)
    pred = next(out.glob("hg_prediction.hg"), None)
    if not pred:
        return None
    row = duckdb.sql(f"SELECT * FROM read_csv('{pred}', delim='\t', header=true, all_varchar=true)").fetchdf()
    if row.empty or row.iloc[0]["Hg"] in (None, "NA"):
        return None
    r = row.iloc[0]
    paths = list(out.glob("*/haplogroup_path_yfull.json"))
    path = json.loads(paths[0].read_text()) if paths else []
    derived = [p["haplogroup"] for p in path if p.get("state") == "D"]
    terminal = derived[-1] if derived else str(r["Hg"]).split("*")[0]
    quality = float(r["QC-score"])
    level = _level(quality, Y_HIGH, Y_MEDIUM)
    return {
        "line": "paternal", "haplogroup": terminal, "quality": quality,
        "call_level": level if bam or level != "High" else "Medium",
        "call_reasons": [f"Yleaf QC-score {quality:.2f} (High ≥ {Y_HIGH}, Medium ≥ {Y_MEDIUM})",
                         f"{int(r['Valid_markers']):,} lineage-tree positions checked in your {source.split(' (')[0]}"]
                        + ([] if bam else ["From a VCF only, so capped at Medium"]),
        "found": derived, "missing": [], "extra": [], "alternatives": [],
        "path": json.dumps([{"haplogroup": p["haplogroup"], "state": p["state"], "marker": p.get("marker")}
                            for p in path]),
        "yleaf_call": r["Hg"], "tool": "Yleaf", "tree": "YFull v14", "input": source,
        "url": YFULL.format(hg=terminal),
    }


# ---------------------------------------------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------------------------------------------


def _reference(ctx: Context) -> tuple[dict, np.ndarray] | None:
    pca = pgs.pgs_dir(ctx) / "pca"
    if not (pca / "reference.pcs").exists() or not (pca / "self.pcs").exists():
        return None
    psam = pgs._panel(ctx)["psam"]
    cols = ", ".join(f"r.PC{i}" for i in range(1, N_PCS + 1))
    con = duckdb.connect()
    ref = con.execute(f"""SELECT r.IID, {cols}, l.SuperPop, l.Population
        FROM read_csv('{pca / "reference.pcs"}', delim='\t', header=true) r
        JOIN read_csv('{psam}', delim='\t', header=true) l ON l."#IID" = r.IID
        ORDER BY r.IID""").fetchnumpy()
    me = con.execute(f"""SELECT {cols.replace("r.", "")} FROM read_csv('{pca / "self.pcs"}', delim='\t',
        header=true)""").fetchnumpy()
    pcs = np.column_stack([ref[f"PC{i}"] for i in range(1, N_PCS + 1)]).astype(float)
    return ({"pcs": pcs, "superpop": np.asarray(ref["SuperPop"], dtype=str),
             "population": np.asarray(ref["Population"], dtype=str), "iid": np.asarray(ref["IID"], dtype=str)},
            np.array([float(me[f"PC{i}"][0]) for i in range(1, N_PCS + 1)]))


def grade_components(rows: list[dict], pca_variants: int | None) -> list[dict]:
    call = "High" if (pca_variants or 0) >= MIN_PCA_VARIANTS else "Medium"
    call_reasons = [f"{pca_variants:,} ancestry-informative positions read from your genome" if pca_variants
                    else "Number of ancestry positions unknown"]
    for r in rows:
        floor = f"{r['noise_floor']:.1%}"
        if not r["detected"]:
            ev, why = "Limited", [f"Below the noise floor ({floor}): reference people with none of this ancestry "
                                  f"get up to this much by chance"]
        elif r["share"] >= 0.5:
            ev, why = "Strong", [f"Your main group; the method gives reference people of this group "
                                 f"{r['reference_own_share']:.0%} on average when they are left out"]
        else:
            ev, why = "Moderate", [f"Above the noise floor ({floor}); shares of mixed ancestry are estimates"]
        r.update({"name": SUPERPOP_NAME.get(r["superpop"], r["superpop"]), "evidence_level": ev,
                  "evidence_reasons": why, "call_level": call, "call_reasons": call_reasons,
                  "overall": _weaker(ev, call)})
    return rows


def _weaker(evidence: str, call: str) -> str:
    rank = {"Strong": 3, "High": 3, "Moderate": 2, "Medium": 2, "Limited": 1, "Low": 1, "Not callable": 0}
    return {3: "Strong", 2: "Moderate", 1: "Limited", 0: "Not callable"}[min(rank[evidence], rank[call])]


def _ancestry(ctx: Context) -> dict:
    d = out_dir(ctx)
    st = Store(ctx.cfg)
    stats = (json.loads(pgs.out_summary(ctx).read_text()).get("genotyping") or {}) if pgs.out_summary(ctx).exists() \
        else {}
    summary: dict = {"method": {"sources": SOURCES, "pcs": N_PCS, "k": K_NEIGHBOURS, "noise_quantile": NOISE_QUANTILE,
                                "confused": CONFUSED, "pca_variants": stats.get("pca_variants")}}
    comp_rows: list[dict] = []
    pop_rows: list[dict] = []
    loaded = _reference(ctx)
    if loaded:
        ref, me = loaded
        comp_rows, cal = components(ref, me)
        comp_rows = grade_components(comp_rows, stats.get("pca_variants"))
        main = max(comp_rows, key=lambda r: r["share"])
        summary.update(cal, main=main["superpop"], main_share=main["share"])
        pop_rows, near = nearest(ref, me, main["superpop"])
        summary["zoom_pcs"] = zoom_pcs(ref, main["superpop"])
        panel_people = None
        hd = _panels(ctx, main["superpop"])
        if hd:
            x, people, mine = hd
            labels = np.array([code(p) for p in people])
            reg = regional(x, labels, mine, main["superpop"])
            pop_rows, near = reg["rows"], reg["near"]
            panel_people = people
            counts = {c: int((labels == c).sum()) for c in set(labels.tolist())}
            panels = {code(p): p["panel_label"] for p in people}
            for r in pop_rows:
                r.update({"panel": panels[r["code"]], "panel_samples": counts[r["code"]],
                          "twin": (near["twins"].get(r["code"]) or {}).get("twin"),
                          "twin_confused": (near["twins"].get(r["code"]) or {}).get("confused")})
            rp = [{"iid": "you", "panel": None, "population": None, "you": True,
                   **{f"pc{i + 1}": float(reg["my_pcs"][i]) for i in range(4)}}]
            rp += [{"iid": p["iid"], "panel": p["panel_label"], "population": code(p), "you": False,
                    **{f"pc{i + 1}": float(v) for i, v in enumerate(row[:4])}}
                   for p, row in zip(people, reg["ref_pcs"], strict=True)]
            pq.write_table(pa.Table.from_pylist(rp), out_region_pca(ctx))
        else:
            out_region_pca(ctx).unlink(missing_ok=True)
        pts = [{"iid": "you", "superpop": None, "population": None, "you": True,
                **{f"pc{i + 1}": float(me[i]) for i in range(N_PCS)}}]
        pts += [{"iid": iid, "superpop": str(sp), "population": str(pp), "you": False,
                 **{f"pc{i + 1}": float(v) for i, v in enumerate(row)}}
                for iid, sp, pp, row in zip(ref["iid"], ref["superpop"], ref["population"], ref["pcs"], strict=True)]
        pq.write_table(pa.Table.from_pylist(pts), out_pca(ctx))
        names = {}
        tbl = st.table("igsr_populations", "populations")
        if tbl:
            names = {r["code"]: r for r in duckdb.sql(f"SELECT * FROM '{tbl}'").fetchdf().to_dict("records")}
        def label(c: str) -> str:
            n = names.get(c) or {}
            return f"{n['name']} ({'HGDP' if c.endswith('HGDP') else c})" if n.get("name") else c

        rec = near["confusion"][near["nearest"]][near["nearest"]]
        mixed = main["share"] < MAIN_SHARE
        broken = [h for h, t in (near.get("twins") or {}).items() if not t["confused"]]
        ev = "Limited" if mixed or broken or len(near["group"]) > 1 or rec < 0.9 else "Moderate"
        reasons = ([f"Your main group is only {main['share']:.0%} of your ancestry, so the nearest population inside "
                    "it is not very meaningful"] if mixed else []) + \
                  ([f"The method cannot tell {', '.join(label(c) for c in near['group'])} apart: in a leave-one-out "
                    "test on the "
                    "reference panel it mixes them up at least "
                    f"{CONFUSED:.0%} of the time"] if len(near["group"]) > 1 else
                   [f"In a leave-one-out test, {rec:.0%} of {label(near['nearest'])} reference people were placed in "
                    "their own population"]) + \
                  ([f"The same population sampled by both reference panels ({', '.join(map(label, broken))}) was "
                    "told apart, so "
                    "the panels may differ for technical reasons"] if broken else [])
        for r in pop_rows:
            n = names.get(r["code"], {})
            r.update({"superpop": main["superpop"], "name": n.get("name"), "description": n.get("description"),
                      "samples": r.get("panel_samples", n.get("samples")), "latitude": n.get("latitude"),
                      "longitude": n.get("longitude"), "url": n.get("url")})
            r.pop("panel_samples", None)
            for k, v in (("panel", "1000 Genomes (phase 3)"), ("twin", None), ("twin_confused", None)):
                r.setdefault(k, v)
        call = main["call_level"]
        if panel_people:
            call = "High" if near["positions"] >= MIN_REGION_SITES else "Medium"
            reasons.append(f"Compared with {len(panel_people):,} high-depth reference genomes (1000 Genomes 30× and "
                           f"the Human Genome Diversity Project) on {near['positions']:,} positions and "
                           f"{near['pcs']} ancestry axes")
        summary["nearest"] = {**{k: v for k, v in near.items() if k != "confusion"}, "evidence_level": ev,
                              "evidence_reasons": reasons, "call_level": call,
                              "overall": _weaker(ev, call), "confusion": near["confusion"],
                              "panels": "high_depth" if panel_people else "phase3"}
        summary["knowledge"] = {k: v for k in ("igsr_populations", "ancestry_panels")
                                if (v := (st.current(k) or {}).get("version"))}
    hap_rows = []
    sex = pgs._sex(ctx) or ""
    for fn, applies in ((maternal, True), (paternal, sex.upper().startswith("XY"))):
        if not applies:
            continue
        try:
            h = fn(ctx, d)
        except Exception as e:  # a missing optional tool or input should not sink the whole section
            summary.setdefault("errors", []).append(f"{fn.__name__}: {e}")
            h = None
        if h:
            h.update({"evidence_level": "Strong",
                      "evidence_reasons": [f"Placed on the {h['tree']}, a published, widely used lineage tree"],
                      "overall": _weaker("Strong", h["call_level"]),
                      **{k: json.dumps(h[k]) for k in ("found", "missing", "extra", "alternatives")}})
            hap_rows.append(h)
    summary["sex"] = sex
    if not sex.upper().startswith("XY"):
        summary["paternal_note"] = "No Y chromosome was detected, so there is no paternal (Y) lineage to read"
    pq.write_table(pa.Table.from_pylist(comp_rows) if comp_rows else pa.table({"superpop": pa.array([], pa.string())}),
                   out_components(ctx))
    pq.write_table(pa.Table.from_pylist(pop_rows) if pop_rows else pa.table({"code": pa.array([], pa.string())}),
                   out_populations(ctx))
    pq.write_table(pa.Table.from_pylist(hap_rows) if hap_rows else pa.table({"line": pa.array([], pa.string())}),
                   out_haplogroups(ctx))
    write_json(out_summary(ctx), summary)
    return {"main": summary.get("main"), "nearest": (summary.get("nearest") or {}).get("group"),
            "haplogroups": {h["line"]: h["haplogroup"] for h in hap_rows}}


LINEAGE_CONTIGS = {"chrM", "MT", "M", "chrMT", "chrY", "Y"}


def _available(ctx: Context) -> str | None:
    if (pgs.pgs_dir(ctx) / "pca" / "self.pcs").exists():
        return None
    if not variants.vcf_file(ctx):
        return "needs a VCF (haplogroups) or the PGS stage's ancestry PCA"
    if ctx.inv.first("cram") or ctx.inv.first("bam"):
        return None
    try:
        lineage = LINEAGE_CONTIGS & set(reference.vcf_contigs(ctx))
    except Exception:  # not indexed yet: let the stage find out
        return None
    # A truth set covering chromosomes 1-22 only (the demo) has nothing to place on a tree; writing empty tables
    # would leave the page and its self-checks querying columns that don't exist
    return None if lineage else ("needs the PGS stage's ancestry PCA, or mtDNA / Y-chromosome calls (this VCF has "
                                 "none)")


def _inputs(ctx: Context) -> list[Path]:
    pca = pgs.pgs_dir(ctx) / "pca"
    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    return [p for p in (pca / "reference.pcs", pca / "self.pcs", pca / "self.bed", pgs.out_summary(ctx),
                        variants.vcf_file(ctx),
                        aln.path if aln else None, qc.out_path(ctx)) if p and p.exists()]


def _params(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    return {"igsr_populations": (st.current("igsr_populations") or {}).get("version"),
            "ancestry_panels": (st.current("ancestry_panels") or {}).get("version"), "mt_tree": MT_TREE,
            "haplo_env": bool(_env_tool("haplo", "haplogrep3")), "sources": SOURCES, "k": K_NEIGHBOURS}


STAGE = Stage(
    name="ancestry", version="3", title="Ancestry and haplogroups", fn=_ancestry, inputs=_inputs,
    outputs=lambda ctx: [out_components(ctx), out_populations(ctx), out_haplogroups(ctx), out_pca(ctx),
                         out_summary(ctx)],  # region_pca.parquet is optional (needs ancestry_panels)
    available=_available, params=_params,
)
