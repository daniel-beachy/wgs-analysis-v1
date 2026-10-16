"""Ancestry (ADR-020): which reference populations your DNA is most like, and your maternal and paternal lineages.

Genetic ancestry here means similarity to present-day reference people, not where ancestors lived or nationality.
Every number is calibrated on the reference panel itself, so the page can say how often the method is wrong.

1. **Continental shares.** Your position on the ancestry PCA (computed by the PGS stage with FRAPOSA on the 1000
   Genomes panel) is written as the closest mix of the five super-population centres (non-negative weights that sum
   to one), the projection approach used by rye (Conley et al. 2023). Each reference person is scored the same way
   with themselves left out; the share the method wrongly gives a group among reference people who are not from it
   (95th percentile) is that group's noise floor, and smaller shares are reported as "not distinguishable from zero".
2. **Nearest reference population** within your main group: the 10 nearest reference people on the 10 PCs. Leave-
   one-out shows which populations the method can tell apart (Finnish vs Tuscan: yes; British vs CEU: no), and
   populations it confuses are reported together.
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
    group = sorted({best} | {p for p in pops if p != best and (conf[best][p] >= CONFUSED or conf[p][best] >= CONFUSED)},
                   key=lambda p: (-votes.get(p, 0), p))
    rows = [{"code": p, "neighbours": votes.get(p, 0), "loo_recall": conf[p][p],
             "confused_with": json.dumps({b: round(v, 3) for b, v in conf[p].items() if b != p and v > 0}),
             "nearest": p == best, "in_nearest_group": p in group} for p in pops]
    return rows, {"nearest": best, "group": group, "k": K_NEIGHBOURS, "accuracy": float(np.mean(pred == y)),
                  "confusion": conf}


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
        rec = near["confusion"][near["nearest"]][near["nearest"]]
        mixed = main["share"] < MAIN_SHARE
        ev = "Limited" if mixed or len(near["group"]) > 1 or rec < 0.9 else "Moderate"
        reasons = ([f"Your main group is only {main['share']:.0%} of your ancestry, so the nearest population inside "
                    "it is not very meaningful"] if mixed else []) + \
                  ([f"The method cannot tell {', '.join(near['group'])} apart: in a leave-one-out test on the "
                    "reference panel it mixes them up at least "
                    f"{CONFUSED:.0%} of the time"] if len(near["group"]) > 1 else
                   [f"In a leave-one-out test, {rec:.0%} of {near['nearest']} reference people were placed in "
                    "their own population"])
        for r in pop_rows:
            n = names.get(r["code"], {})
            r.update({"superpop": main["superpop"], "name": n.get("name"), "description": n.get("description"),
                      "samples": n.get("samples"), "latitude": n.get("latitude"), "longitude": n.get("longitude"),
                      "url": n.get("url")})
        summary["nearest"] = {**{k: v for k, v in near.items() if k != "confusion"}, "evidence_level": ev,
                              "evidence_reasons": reasons, "call_level": main["call_level"],
                              "overall": _weaker(ev, main["call_level"]), "confusion": near["confusion"]}
        summary["knowledge"] = {"igsr_populations": (st.current("igsr_populations") or {}).get("version")}
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


def _available(ctx: Context) -> str | None:
    if not variants.vcf_file(ctx) and not (pgs.pgs_dir(ctx) / "pca" / "self.pcs").exists():
        return "needs a VCF (haplogroups) or the PGS stage's ancestry PCA"
    return None


def _inputs(ctx: Context) -> list[Path]:
    pca = pgs.pgs_dir(ctx) / "pca"
    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    return [p for p in (pca / "reference.pcs", pca / "self.pcs", pgs.out_summary(ctx), variants.vcf_file(ctx),
                        aln.path if aln else None, qc.out_path(ctx)) if p and p.exists()]


def _params(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    return {"igsr_populations": (st.current("igsr_populations") or {}).get("version"), "mt_tree": MT_TREE,
            "haplo_env": bool(_env_tool("haplo", "haplogrep3")), "sources": SOURCES, "k": K_NEIGHBOURS}


STAGE = Stage(
    name="ancestry", version="2", title="Ancestry and haplogroups", fn=_ancestry, inputs=_inputs,
    outputs=lambda ctx: [out_components(ctx), out_populations(ctx), out_haplogroups(ctx), out_pca(ctx),
                         out_summary(ctx)],
    available=_available, params=_params,
)
