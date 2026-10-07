"""Pharmacogenomics (ADR-015): how your genes change the way medicines work for you.

Three stages:

* ``cyp2d6`` — CYP2D6 metabolises ~20 % of prescribed drugs but sits next to two look-alike pseudogenes and is
  often duplicated or deleted, so ordinary variant calls cannot type it. Cyrius (Illumina) reads the CRAM
  and calls star alleles *with copy number*. Slow on a slow drive, so it is its own cached stage.
* ``hla`` — HLA-A/HLA-B alleles (e.g. HLA-B*57:01 → abacavir hypersensitivity). T1K types them from the reads
  in the MHC region plus unmapped reads, against the IPD-IMGT/HLA catalogue (a knowledge source).
* ``pgx`` — PharmCAT. Its ~1,200 GRCh38 positions are lifted to your build (knowledge source ``pharmcat``),
  each is genotyped from your VCF *and* the callable map (absent + callable = reference; absent + not
  callable = missing, never assumed reference), REF bases are checked against your reference, and a
  GRCh38-coordinate VCF is handed to PharmCAT together with the outside calls above (CYP2D6, HLA, and
  MT-RNR1 read from your mitochondrial variants). PharmCAT returns diplotypes, phenotypes and the matching
  CPIC / DPWG / FDA guidance, which become graded claims.

Results of the same tools run elsewhere (e.g. on Linux) are ingested from ``external/`` when present.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path

import duckdb

from .. import evidence as ev
from ..knowledge import Store
from ..knowledge.hla import INDEX as HLA_INDEX
from ..knowledge.hla import t1k_tool
from ..knowledge.pharmcat import JAR, POSITIONS_VCF
from ..pipeline import Context, Stage, run, write_json
from . import coverage, qc, reference, variants

MIN_GQ = 20
# MHC (HLA) region, 1-based inclusive; T1K also needs unmapped reads (HLA alleles far from the reference).
MHC = {"GRCh37": (28477797, 33448354), "GRCh38": (28510120, 33480577)}
PAR_X37 = [(60001, 2699520), (154931044, 155260560)]
PAR_X38 = [(10001, 2781479), (155701383, 156030895)]
MT_RNR1 = (648, 1601)
PGX_HLA = ("HLA-A", "HLA-B")


def pgx_dir(ctx: Context) -> Path:
    p = ctx.work / "pgx"
    p.mkdir(parents=True, exist_ok=True)
    return p


def _build(ctx: Context) -> str | None:
    f = ctx.inv.first("cram") or ctx.inv.first("bam")
    return (f.meta or {}).get("build") if f else None


def _chr_prefixed(ctx: Context) -> bool:
    return "chr1" in reference.cram_contigs(ctx)


def _env_tool(env: str, name: str) -> str | None:
    """A tool from a sibling pixi environment (Cyrius lives in its own env to keep its pinned deps apart)."""
    p = Path(sys.prefix).parent / env / "bin" / name
    return str(p) if p.exists() else shutil.which(name)


def _java() -> str | None:
    p = Path(sys.prefix) / "lib" / "jvm" / "bin" / "java"
    return str(p) if p.exists() else shutil.which("java")


# ---------------------------------------------------------------------------------------------------------------
# CYP2D6 with Cyrius
# ---------------------------------------------------------------------------------------------------------------


def cyrius_out(ctx: Context) -> Path:
    return pgx_dir(ctx) / "cyrius.json"


def _cyrius(ctx: Context) -> dict:
    cram = ctx.inv.first("cram") or ctx.inv.first("bam")
    build = {"GRCh37": "19" if _chr_prefixed(ctx) else "37", "GRCh38": "38"}[_build(ctx)]
    tmp = ctx.scratch / "cyrius"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    (tmp / "manifest.txt").write_text(str(cram.path.resolve()) + "\n")
    run([_env_tool("cyrius", "cyrius"), "-m", str(tmp / "manifest.txt"), "-g", build, "-o", str(tmp), "-p",
         "cyrius", "-t", str(min(ctx.threads, 8)), "-r", str(reference.reference_fasta(ctx))],
        log=ctx.logs / "cyrius.log")
    raw = json.loads((tmp / "cyrius.json").read_text())
    (_, res), = raw.items()
    if not res.get("Genotype"):
        raise RuntimeError("Cyrius produced no genotype — see logs/cyrius.log")
    log = (ctx.logs / "cyrius.log").read_text()
    out = {"tool": "Cyrius", "genotype": res.get("Genotype"), "filter": res.get("Filter"),
           "total_cn": res.get("Total_CN"), "raw_star_alleles": res.get("Raw_star_allele"),
           "call_info": res.get("Call_info"), "median_depth": res.get("Median_depth"),
           "coverage_mad": res.get("Coverage_MAD"), "variants_called": res.get("Variants_called"),
           "warnings": ["Uneven coverage: Cyrius warns its copy-number calls may be unreliable"]
           if "uneven coverage" in log.lower() else []}
    write_json(cyrius_out(ctx), out)
    return {"genotype": out["genotype"], "filter": out["filter"]}


def _cyrius_available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "needs a CRAM/BAM"
    if not reference.reference_fasta(ctx):
        return "needs a verified reference"
    if _build(ctx) not in ("GRCh37", "GRCh38"):
        return "Cyrius supports GRCh37/GRCh38 alignments only"
    if not _env_tool("cyrius", "cyrius"):
        return "Cyrius is not installed (pixi install -e cyrius)"
    return None


CYRIUS = Stage(
    name="cyp2d6", version="1", title="CYP2D6 star alleles + copy number (Cyrius)", fn=_cyrius,
    inputs=lambda ctx: [f.path for f in (ctx.inv.of("cram") or ctx.inv.of("bam"))[:1]],
    outputs=lambda ctx: [cyrius_out(ctx)], available=_cyrius_available,
)


# ---------------------------------------------------------------------------------------------------------------
# HLA with T1K
# ---------------------------------------------------------------------------------------------------------------


def hla_out(ctx: Context) -> Path:
    return pgx_dir(ctx) / "hla.json"


def parse_t1k(path: Path) -> dict:
    """T1K *_genotype.tsv → {gene: {alleles: [...], quality: [...]}} at two-field (protein) resolution."""
    genes = {}
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    for row in rows:
        if len(row) < 8 or not row[0].startswith("HLA-"):
            continue
        alleles, quals = [], []
        for name, qual in ((row[2], row[4]), (row[5], row[7])):
            if name not in (".", ""):
                alleles.append("*" + ":".join(name.split("*", 1)[1].split(":")[:2]))
                quals.append(int(float(qual)))
        if int(row[1] or 0) == 1 and alleles:  # one allele reported = homozygous
            alleles, quals = alleles * 2, quals * 2
        genes[row[0]] = {"alleles": alleles, "quality": quals}
    return genes


def _hla(ctx: Context) -> dict:
    cram = ctx.inv.first("cram") or ctx.inv.first("bam")
    fasta = reference.reference_fasta(ctx)
    start, end = MHC[_build(ctx)]
    region = f"{'chr6' if _chr_prefixed(ctx) else '6'}:{start}-{end}"
    tmp = ctx.scratch / "hla"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    t = str(min(ctx.threads, 8))
    st = shutil.which("samtools", path=str(Path(sys.prefix) / "bin")) or "samtools"
    run(["bash", "-o", "pipefail", "-c",
         f"'{st}' view -@ {t} -T '{fasta}' -u -F 0x900 '{cram.path}' {region} '*' | "
         f"'{st}' collate -@ {t} -u -O - '{tmp}/collate' | "
         f"'{st}' fastq -@ {t} -1 '{tmp}/r1.fq.gz' -2 '{tmp}/r2.fq.gz' -s /dev/null -0 /dev/null -n -"],
        log=ctx.logs / "hla.log")
    index = Store(ctx.cfg).file("imgt_hla", HLA_INDEX)
    run([t1k_tool("run-t1k"), "-1", str(tmp / "r1.fq.gz"), "-2", str(tmp / "r2.fq.gz"), "--preset", "hla-wgs",
         "-f", str(index), "-t", t, "-o", "t1k", "--od", str(tmp / "out")], log=ctx.logs / "hla.log")
    shutil.copyfile(tmp / "out" / "t1k_genotype.tsv", pgx_dir(ctx) / "t1k_genotype.tsv")
    genes = parse_t1k(tmp / "out" / "t1k_genotype.tsv")
    out = {"tool": "T1K", "imgt_hla": Store(ctx.cfg).versions().get("imgt_hla"), "region": region, "genes": genes}
    write_json(hla_out(ctx), out)
    shutil.rmtree(tmp, ignore_errors=True)
    return {g: "/".join(genes.get(g, {}).get("alleles", [])) for g in PGX_HLA}


def _hla_available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "needs a CRAM/BAM"
    if not reference.reference_fasta(ctx):
        return "needs a verified reference"
    if _build(ctx) not in MHC:
        return "HLA region is defined for GRCh37/GRCh38 only"
    if not t1k_tool("run-t1k"):
        return "T1K is not installed on this platform (Apple Silicon only; or put a T1K result in external/)"
    if not Store(ctx.cfg).file("imgt_hla", HLA_INDEX):
        return "needs the IPD-IMGT/HLA knowledge source (wgs knowledge refresh --source imgt_hla)"
    return None


def _hla_inputs(ctx: Context) -> list[Path]:
    idx = Store(ctx.cfg).file("imgt_hla", HLA_INDEX)
    return [f.path for f in (ctx.inv.of("cram") or ctx.inv.of("bam"))[:1]] + ([idx] if idx else [])


HLA = Stage(
    name="hla", version="1", title="HLA-A/B/C typing (T1K)", fn=_hla, inputs=_hla_inputs,
    outputs=lambda ctx: [hla_out(ctx)], available=_hla_available,
)


# ---------------------------------------------------------------------------------------------------------------
# PharmCAT
# ---------------------------------------------------------------------------------------------------------------


def norm(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Minimal representation: trim shared trailing then leading bases (so CAT>C at 10 == AT>'' at 11)."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def genotype_position(p: dict, ref37: str | None, rows: list[dict], callable_: bool, haploid: bool) -> dict:
    """Genotype one PharmCAT position from your VCF rows overlapping it.

    `p`: {pos37, ref, alts}. `ref37`: your reference's bases over the PharmCAT REF span. `rows`: VCF records
    (pos, ref, alt, filter, gt, gq) overlapping the span. Returns {gt, status, note} with gt in PharmCAT's
    allele indices ('0/1', '1', './.')."""
    ref, alts, pos = p["ref"], list(p["alts"]), p["pos37"]
    swapped = None
    if ref37 is None:
        return {"gt": "./.", "status": "no_reference", "note": "position not on your reference"}
    if ref37.upper() != ref.upper():
        if len(ref) == 1 and ref37.upper() in [a.upper() for a in alts]:
            swapped = [a.upper() for a in alts].index(ref37.upper()) + 1
        else:
            return {"gt": "./.", "status": "ref_mismatch",
                    "note": f"PharmCAT expects {ref}, your reference has {ref37}"}
    home = swapped or 0  # the allele your reference carries, in PharmCAT's numbering
    wanted = {norm(pos, ref, a): i + 1 for i, a in enumerate(alts)}
    if swapped:
        wanted[norm(pos, ref37, ref)] = 0

    def missing(status: str, note: str) -> dict:
        return {"gt": "./.", "status": status, "note": note}

    matched: list[tuple[int, str]] = []
    refcall = False
    for r in rows:
        gt, flt = r["gt"] or "./.", r["filter"] or "."
        if flt == "RefCall" or gt in ("0/0", "0|0"):
            if (r["gq"] or 0) >= MIN_GQ:
                refcall = True
            continue
        if "." in gt or flt == "NoCall":
            return missing("no_call", "the variant caller could not decide here")
        key = norm(r["pos"], r["ref"], r["alt"])
        if key not in wanted:
            return missing("other_variant", f"you carry a different variant here ({r['ref']}>{r['alt']})")
        if flt not in ("PASS", "."):
            return missing("filtered", f"variant call flagged {flt}")
        if (r["gq"] or 0) < MIN_GQ:
            return missing("low_quality", f"genotype quality {r['gq']} < {MIN_GQ}")
        matched.append((wanted[key], gt))
    if not matched:
        if refcall or callable_:
            alleles = [home, home]
            status = "reference"
        else:
            return missing("not_callable", "too few good reads here to tell")
    elif len(matched) == 1:
        idx, gt = matched[0]
        alleles = [idx, idx] if gt.replace("|", "/") == "1/1" else sorted([home, idx])
        status = "variant"
    elif len(matched) == 2:
        alleles = sorted(i for i, _ in matched)
        status = "variant"
    else:
        return missing("complex", "more than two alleles reported")
    if haploid:
        if alleles[0] != alleles[1]:
            return missing("het_haploid", "two different alleles on a single-copy X chromosome — unreliable")
        return {"gt": str(alleles[0]), "status": status, "note": None}
    return {"gt": f"{alleles[0]}/{alleles[1]}", "status": status, "note": None}


def _sample_sex(ctx: Context) -> str | None:
    p = qc.out_path(ctx)
    return (json.loads(p.read_text()).get("inferred_sex") or {}).get("call") if p.exists() else None


def _haploid_x(ctx: Context, chrom: str, pos: int) -> bool:
    if chrom.removeprefix("chr") != "X" or not (_sample_sex(ctx) or "").startswith("XY"):
        return False
    par = PAR_X38 if _build(ctx) == "GRCh38" else PAR_X37
    return not any(a <= pos <= b for a, b in par)


def _genotype_all(ctx: Context, positions: Path) -> list[dict]:
    import pysam

    vp, cp = variants.variants_parquet(ctx), coverage.out_dir(ctx) / "callable.parquet"
    con = duckdb.connect()
    pos = con.execute(f"SELECT * FROM '{positions}' ORDER BY chrom38, pos38").fetchall()
    cols = [d[0] for d in con.description]
    pos = [dict(zip(cols, r, strict=True)) for r in pos]
    con.execute("CREATE TABLE q (i INTEGER, chrom VARCHAR, s INTEGER, e INTEGER)")
    con.executemany("INSERT INTO q VALUES (?, ?, ?, ?)",
                    [[i, (p["chrom37"] or "").removeprefix("chr"), p["pos37"] or 0,
                      (p["pos37"] or 0) + len(p["ref"]) - 1] for i, p in enumerate(pos)])
    rows: dict[int, list[dict]] = {}
    for r in con.execute(f"""SELECT q.i, v.pos, v.ref, v.alt, v.filter, v.gt, v.gq FROM q
            JOIN '{vp}' v ON v.chrom = q.chrom AND v.pos BETWEEN q.s - 60 AND q.e
             AND v.pos + length(v.ref) - 1 >= q.s ORDER BY q.i, v.pos""").fetchall():
        rows.setdefault(r[0], []).append(dict(zip(("pos", "ref", "alt", "filter", "gt", "gq"), r[1:], strict=True)))
    callable_ = {r[0] for r in con.execute(f"""SELECT q.i FROM q JOIN '{cp}' c ON c.chrom = q.chrom
            AND c.start < q.s AND c."end" >= q.e AND c.state = 'CALLABLE'""").fetchall()}
    fa = pysam.FastaFile(str(reference.reference_fasta(ctx)))
    prefixed = "chr1" in fa.references
    out = []
    for i, p in enumerate(pos):
        ref37 = None
        if p["chrom37"]:
            name = p["chrom37"] if prefixed else p["chrom37"].removeprefix("chr")
            name = "chrM" if name in ("chrMT",) else name
            if name in fa.references:
                ref37 = fa.fetch(name, p["pos37"] - 1, p["pos37"] - 1 + len(p["ref"]))
        call = genotype_position(p, ref37, rows.get(i, []), i in callable_, _haploid_x(ctx, p["chrom37"] or "",
                                                                                         p["pos37"] or 0))
        out.append({**p, **call})
    return out


def _write_vcf(path: Path, header_vcf: Path, calls: list[dict], sample: str) -> None:
    lines = []
    for line in header_vcf.read_text().splitlines():
        if line.startswith("##"):
            lines.append(line)
        elif line.startswith("#CHROM"):
            lines.append("\t".join(line.split("\t")[:9] + [sample]))
    lines.extend("\t".join([c["chrom38"], str(c["pos38"]), c["rsid"] or ".", c["ref"], ",".join(c["alts"]) or ".",
                            ".", "PASS", f"PX={c['gene']}" if c["gene"] else ".", "GT", c["gt"]]) for c in calls)
    path.write_text("\n".join(lines) + "\n")


def _mt_rnr1_names(jar: Path) -> dict[str, str]:
    with zipfile.ZipFile(jar) as z:
        d = json.loads(z.read("org/pharmgkb/pharmcat/phenotype/MT_RNR1.json"))
    return d["haplotypes"]


def _mt_rnr1(ctx: Context, jar: Path) -> dict | None:
    """MT-RNR1 from your mitochondrial variants (PharmCAT cannot read chrM from a VCF itself)."""
    names = _mt_rnr1_names(jar)
    con = duckdb.connect()
    vp, cp = variants.variants_parquet(ctx), coverage.out_dir(ctx) / "callable.parquet"
    a, b = MT_RNR1
    found = []
    for pos, ref, alt, gt, vaf in con.execute(f"""SELECT pos, ref, alt, gt, vaf FROM '{vp}'
            WHERE chrom = 'MT' AND pos BETWEEN {a - 5} AND {b} AND filter = 'PASS' AND gt NOT IN ('0/0', './.')
            ORDER BY pos""").fetchall():
        p, r, x = norm(pos, ref, alt)
        name = f"m.{p}{r}>{x}" if r and x else f"m.{p}{r}>del" if r else None
        found.append({"name": name, "known": name in names, "vaf": vaf, "gt": gt})
    known = [f for f in found if f["known"]]
    covered = con.execute(f"""SELECT coalesce(sum(least("end", {b}) - greatest(start, {a - 1})), 0) FROM '{cp}'
            WHERE chrom = 'MT' AND state IN ('CALLABLE', 'HIGH') AND start < {b} AND "end" > {a - 1}""").fetchone()[0]
    frac = covered / (b - a + 1)
    if known:
        known.sort(key=lambda f: ("Increased" not in names[f["name"]], f["name"]))
        call = known[0]["name"]
    elif frac >= 0.95:
        call = "Reference"
    else:
        return None
    return {"call": call, "function": names.get(call), "variants_in_gene": found, "callable_fraction": round(frac, 3),
            "heteroplasmic": [f["name"] for f in known if (f["vaf"] or 1) < 0.9]}


def _external(ctx: Context, tool: str, suffix: str) -> Path | None:
    hits = [f.path for f in ctx.inv.of(f"external:{tool}") if f.path.name.endswith(suffix)]
    return sorted(hits, key=lambda p: p.stat().st_mtime)[-1] if hits else None


def _external_cyrius(path: Path) -> dict | None:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        if row.get("Genotype") and row["Genotype"] != "None":
            return {"tool": "Cyrius (external)", "genotype": row["Genotype"], "filter": row.get("Filter"),
                    "warnings": [], "file": path.name}
    return None


def outside_calls(ctx: Context, jar: Path) -> dict[str, dict]:
    """{gene: {call, tool, detail}} — calls PharmCAT cannot make from a VCF."""
    out: dict[str, dict] = {}
    cy = json.loads(cyrius_out(ctx).read_text()) if cyrius_out(ctx).exists() else None
    if not cy and (ext := _external(ctx, "cyrius", ".tsv")):
        cy = _external_cyrius(ext)
    if cy and cy.get("genotype") and cy["genotype"] != "None":
        # Cyrius writes CN variants like *1/*2x2 or *4+*68; PharmCAT accepts the same notation.
        out["CYP2D6"] = {"call": cy["genotype"], "tool": cy["tool"], "detail": cy}
    hla = json.loads(hla_out(ctx).read_text()) if hla_out(ctx).exists() else None
    if not hla and (ext := _external(ctx, "t1k", "_genotype.tsv")):
        hla = {"tool": "T1K (external)", "genes": parse_t1k(ext), "file": ext.name}
    if hla:
        for g in PGX_HLA:
            al = hla["genes"].get(g, {}).get("alleles", [])
            if len(al) == 2:
                out[g] = {"call": "/".join(al), "tool": hla["tool"],
                          "detail": {"quality": hla["genes"][g].get("quality"), "imgt_hla": hla.get("imgt_hla")}}
    mt = _mt_rnr1(ctx, jar)
    if mt:
        out["MT-RNR1"] = {"call": mt["call"], "tool": "mitochondrial variants (this pipeline)", "detail": mt}
    return out


def out_paths(ctx: Context) -> dict[str, Path]:
    d = pgx_dir(ctx)
    return {"genes": d / "pgx_genes.parquet", "drugs": d / "pgx_drugs.parquet", "claims": d / "claims_pgx.parquet",
            "positions": d / "pgx_positions.parquet", "summary": d / "pgx.json", "report": d / "pharmcat.report.html"}


def _run_pharmcat(ctx: Context, jar: Path, vcf: Path, oc: Path | None, outdir: Path) -> Path:
    # The JVM reads the jar with many small random reads; on slow external drives that alone takes minutes.
    local = ctx.scratch / "pharmcat.jar"
    if not local.exists() or local.stat().st_size != jar.stat().st_size:
        shutil.copyfile(jar, local)
    jar = local
    cmd = [_java(), "-jar", str(jar), "-vcf", str(vcf), "-o", str(outdir), "-bf", "pharmcat", "-reporterJson",
           "-reporterHtml"]
    if oc:
        cmd[5:5] = ["-po", str(oc)]
    outdir.mkdir(parents=True, exist_ok=True)
    run(cmd, log=ctx.logs / "pharmcat.log", cwd=outdir)  # PharmCAT also writes a log into its cwd
    return outdir / "pharmcat.report.json"


def _pgx(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    jar, positions, header = st.file("pharmcat", JAR), st.table("pharmcat", "positions"), \
        st.file("pharmcat", POSITIONS_VCF)
    paths = out_paths(ctx)
    tmp = ctx.scratch / "pharmcat"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    calls = _genotype_all(ctx, positions)
    vcf = tmp / "pharmcat_input.vcf"
    _write_vcf(vcf, header, calls, ctx.sample)
    oc = outside_calls(ctx, jar)
    oc_file = None
    if oc:
        oc_file = tmp / "outside_calls.tsv"
        oc_file.write_text("".join(f"{g}\t{v['call']}\n" for g, v in oc.items()))
    report = _run_pharmcat(ctx, jar, vcf, oc_file, tmp / "out")
    shutil.copyfile(tmp / "out" / "pharmcat.report.html", paths["report"])
    with gzip.open(pgx_dir(ctx) / "pharmcat.report.json.gz", "wt") as fh:
        fh.write(report.read_text())
    versions = {"pharmcat": st.versions().get("pharmcat")}
    if "HLA-A" in oc or "HLA-B" in oc:
        versions["imgt_hla"] = st.versions().get("imgt_hla")
    return write_tables(ctx, json.loads(report.read_text()), calls, oc, versions)


# --- report.json → tables and graded claims -------------------------------------------------------------------

SOURCE_LABEL = {"CPIC Guideline Annotation": "CPIC", "DPWG Guideline Annotation": "DPWG",
                "FDA Label Annotation": "FDA label", "FDA PGx Association": "FDA association"}
NORMAL = re.compile(r"^(normal|.*negative|n/a)", re.I)
STANDARD = re.compile(r"^(use .* per standard|initiate therapy with (standard )?(recommended )?(starting )?dose|"
                      r"initiate standard dosing|use ([\w-]+ )?label recommended|no recommendation|"
                      r"it is not possible to (offer|provide)|"
                      r"prescribe desired starting dose|no adjustments? needed|no reason to avoid|"
                      r"there is no need to avoid|label recommended|no action|"
                      r"the guideline does not provide a recommendation|"
                      r"clinical findings, family history, further genetic testing)", re.I)
NO_RESULT = {"No Result", "n/a", "Indeterminate", ""}
GENOTYPE_ONLY = "Genotype reported (guidelines use it directly)"

PHENO_CLASS = [("poor", "Poor"), ("intermediate", "Intermediate"), ("ultrarapid", "Ultrarapid"), ("rapid", "Rapid"),
               ("normal", "Normal"), ("decreased", "Decreased"), ("increased", "Increased"), ("positive", "Positive"),
               ("negative", "Negative"), ("possible", "Possible"), ("deficien", "Deficient"),
               ("indeterminate", "Indeterminate"), ("risk", "Risk")]


def is_standard(rec: str) -> bool:
    """True when a recommendation just says "use as normal" (a leading "For PAIN:"-style scope is ignored)."""
    return bool(STANDARD.search(re.sub(r"^for [a-z ,]+:\s*", "", rec.strip(), flags=re.I)))


def _pheno_class(phenotype: str) -> str:
    low = phenotype.lower()
    if not phenotype or phenotype in NO_RESULT:
        return "no_result"
    if phenotype == GENOTYPE_ONLY:
        return "genotype"
    if "increased risk" in low or "positive" in low and "negative" not in low:
        return "risk"
    for key, _ in PHENO_CLASS:
        if key in low:
            return key
    return "other"


def _gene_call_level(gene: str, g: dict, oc: dict, pos: dict) -> tuple[str, list[str]]:
    """Call confidence for a diplotype — how sure we are it is really yours."""
    if gene in oc:
        o = oc[gene]
        if gene == "CYP2D6":
            reasons = ["Called by Cyrius, which handles CYP2D6's copy-number changes and look-alike pseudogene",
                       "Cyrius was validated on Illumina reads; tellmeGen used DNBSEQ (similar short reads, but not "
                       "formally validated)"]
            d = o["detail"]
            if d.get("filter") not in ("PASS", None):
                return "Low", reasons + [f"Cyrius flagged the call ({d.get('filter')})"]
            if d.get("total_cn") is not None:
                reasons.append(f"{d['total_cn']} copies of CYP2D6+CYP2D7 counted (2 CYP2D6 = no deletion or "
                               "duplication)" if d["total_cn"] == 4 else f"Copy number: CYP2D6+CYP2D7 total "
                                                                         f"{d['total_cn']}")
            reasons.extend(d.get("warnings", []))
            if len(d.get("raw_star_alleles") or []) > 1:
                reasons.append("More than one star-allele combination fits the reads; Cyrius picked the most common "
                               f"({', '.join(d['raw_star_alleles'])})")
            return "Medium", reasons
        if gene.startswith("HLA-"):
            q = o["detail"].get("quality") or []
            reasons = ["Typed by T1K from reads in the HLA region against the IPD-IMGT/HLA catalogue",
                       "HLA genes are the most variable in the genome; short-read typing is ~95–99 % accurate at "
                       "this (two-field) resolution, below a clinical HLA lab"]
            if q and min(q) < 20:
                return "Low", reasons + [f"Low T1K allele quality ({min(q)})"]
            return "Medium", reasons
        d = o["detail"]
        reasons = [f"Read directly from your mitochondrial DNA (callable over {d['callable_fraction']:.0%} of "
                   "MT-RNR1)"]
        if d.get("heteroplasmic"):
            return "Medium", reasons + ["Variant present in only some of your mitochondria (heteroplasmy)"]
        return "High", reasons
    total = pos.get("total", 0)
    miss = pos.get("missing", 0)
    reasons = [f"{total - miss} of {total} defining positions genotyped from your data"]
    if not total:
        return "Not callable", ["No positions for this gene"]
    if g.get("uncalled"):
        reasons.append("Some alleles could not be ruled out because positions are missing: "
                       + ", ".join(g["uncalled"][:8]))
    if g.get("ambiguous"):
        reasons.append("More than one diplotype fits your genotypes (your VCF is not phased)")
    if miss == total:
        return "Not callable", reasons
    if miss / total > 0.2 or g.get("ambiguous"):
        return "Low", reasons
    if miss or g.get("uncalled"):
        return "Medium", reasons
    reasons.append("Every defining position called with good quality")
    return "High", reasons


def _cpic_evidence(classification: str | None, source: str) -> tuple[str, list[str]]:
    c = (classification or "").strip()
    if source == "CPIC":
        level = {"Strong": "Strong", "Moderate": "Moderate"}.get(c, "Limited")
        return level, [f"CPIC guideline, recommendation strength “{c or 'unspecified'}”"]
    if source == "DPWG":
        return "Moderate", ["DPWG (Dutch Pharmacogenetics Working Group) guideline — evidence-graded expert advice"]
    if source == "FDA label":
        return "Moderate", ["FDA-approved drug label pharmacogenomic information"]
    return "Limited", ["FDA table of pharmacogenetic associations"]


def write_tables(ctx: Context, report: dict, calls: list[dict], oc: dict, versions: dict) -> dict:
    paths = out_paths(ctx)
    pos_by_gene: dict[str, dict] = {}
    for c in calls:
        g = pos_by_gene.setdefault(c["gene"], {"total": 0, "missing": 0, "variant": 0, "missing_rsids": []})
        g["total"] += 1
        if c["gt"].startswith("."):
            g["missing"] += 1
            g["missing_rsids"].append(c["rsid"] or f"{c['chrom38']}:{c['pos38']}")
        elif c["status"] == "variant":
            g["variant"] += 1

    genes, claims = [], []
    drug_rows = []
    for src_key, drugs in report.get("drugs", {}).items():
        src = SOURCE_LABEL.get(src_key, src_key)
        for name, d in drugs.items():
            for gl in d.get("guidelines", []):
                for a in gl.get("annotations", []):
                    dips = [dp for gt in a.get("genotypes", []) for dp in gt.get("diplotypes", [])]
                    gene_list = sorted({dp["gene"] for dp in dips} | {k for lk in a.get("lookupKey") or []
                                                                      for k in (lk or {})})
                    phenos = [f"{k}: {v}" for lk in a.get("lookupKey") or [] for k, v in (lk or {}).items()]
                    rec = (a.get("drugRecommendation") or "").strip()
                    normal = all(NORMAL.match(str(v)) for lk in a.get("lookupKey") or [] for v in (lk or {}).values())
                    # CFTR guidance only applies to people who have cystic fibrosis; reference CFTR is not a finding.
                    cf_only = gene_list == ["CFTR"] and all("Reference" in (dp.get("label") or "") for dp in dips)
                    action = bool(rec) and not cf_only and not is_standard(rec) and (
                        a.get("dosingInformation") or a.get("alternateDrugAvailable") or not normal)
                    drug_rows.append({
                        "drug": name, "drug_id": d.get("id"), "source": src, "guideline": gl.get("name"),
                        "url": gl.get("url") or (d.get("urls") or [None])[0], "genes": gene_list,
                        "phenotypes": phenos, "diplotypes": sorted({dp.get("label") for dp in dips}),
                        "population": a.get("population"), "classification": a.get("classification"),
                        "recommendation": rec or None, "implications": a.get("implications") or [],
                        "action": action, "matched": bool(rec),
                        # Unmatched rows (e.g. CPIC warfarin: "use the flowchart") carry their advice as drug notes.
                        "messages": [m.get("message") for m in (a.get("messages") or ([] if rec else d.get("messages"))
                                                                 or []) if m.get("message")][:5],
                        "citations": json.dumps([{"pmid": c.get("pmid"), "title": c.get("title"), "year": c.get("year")}
                                                 for c in d.get("citations") or []][:6]),
                    })
                if not gl.get("annotations"):
                    drug_rows.append({"drug": name, "drug_id": d.get("id"), "source": src, "guideline": gl.get("name"),
                                      "url": gl.get("url"), "genes": [], "phenotypes": [], "diplotypes": [],
                                      "population": None, "classification": None, "recommendation": None,
                                      "implications": [], "action": False, "matched": False,
                                      "messages": [m.get("message") for m in d.get("messages") or []][:5],
                                      "citations": "[]"})
    cpic_genes = {g for r in drug_rows if r["source"] == "CPIC" for g in r["genes"]}

    for gene, g in sorted(report.get("genes", {}).items()):
        dips = g.get("recommendationDiplotypes") or g.get("sourceDiplotypes") or []
        d0 = dips[0] if dips else {}
        label = " or ".join(dp.get("label") for dp in dips) if dips else "Unknown"
        phen = "; ".join(p for dp in dips[:1] for p in dp.get("phenotypes") or [])
        no_result = phen in NO_RESULT - {""} or not dips or label.startswith("Unknown")
        if no_result:
            phen = "No Result"
        elif not phen:
            # VKORC1, CYP4F2, IFNL3…: guidelines act on the genotype itself; there is no phenotype word.
            phen = GENOTYPE_ONLY
        info = {"uncalled": g.get("uncalledHaplotypes") or [], "ambiguous": len(dips) > 1}
        if no_result:
            call_level, call_reasons = "Not callable", (
                ["Not typed: " + ("needs a dedicated caller (see Medicines notes)" if gene in ("CYP2D6", "HLA-A",
                 "HLA-B", "MT-RNR1") else "too many defining positions missing")])
        else:
            call_level, call_reasons = _gene_call_level(gene, info, oc, pos_by_gene.get(gene, {}))
        evidence = "Strong" if gene in cpic_genes else "Moderate"
        ev_reasons = (["CPIC standardised how these alleles translate to a phenotype (allele function from "
                       "PharmVar/CPIC tables)"] if evidence == "Strong" else
                      ["Allele function from DPWG/FDA resources (no CPIC guideline for this gene)"])
        overall = ev.overall(evidence, call_level)
        a1, a2 = d0.get("allele1") or {}, d0.get("allele2") or {}
        related = [r["name"] for r in g.get("relatedDrugs") or []]
        variants_ = [{"rsid": v.get("dbSnpId"), "pos": v.get("position"), "call": v.get("call"),
                      "alleles": v.get("alleles"), "ref": v.get("referenceAllele")}
                     for v in g.get("variants") or []]
        variants_of_interest = [v for v in variants_ if v["call"] and v["ref"] and
                                any(b != v["ref"] for b in re.split(r"[/|]", v["call"]) if b)]
        source = oc[gene]["tool"] if gene in oc else ("PharmCAT (from your VCF)" if not no_result else None)
        row = {"gene": gene, "diplotype": None if no_result else label, "allele1": a1.get("name"),
               "allele2": a2.get("name"), "allele1_function": a1.get("function"),
               "allele2_function": a2.get("function"), "phenotype": phen,
               "phenotype_class": _pheno_class(phen), "activity_score": d0.get("activityScore"),
               "call_source": source, "chrom": g.get("chr"),
               "positions_total": pos_by_gene.get(gene, {}).get("total", 0),
               "positions_missing": pos_by_gene.get(gene, {}).get("missing", 0),
               "positions_variant": pos_by_gene.get(gene, {}).get("variant", 0),
               "missing_rsids": pos_by_gene.get(gene, {}).get("missing_rsids", [])[:50],
               "uncalled_haplotypes": info["uncalled"][:50],
               "variants_of_interest": json.dumps(variants_of_interest),
               "messages": [m.get("message") for m in g.get("messages") or [] if m.get("message")][:8],
               "related_drugs": related, "evidence_level": evidence, "evidence_reasons": ev_reasons,
               "call_level": call_level, "call_reasons": call_reasons, "overall": overall,
               "outside_detail": json.dumps(oc[gene]["detail"], default=str) if gene in oc else None}
        genes.append(row)
        claims.append({
            "claim_id": f"pgx:gene:{gene}", "section": "pgx", "group": "pgx_gene", "kind": "pgx_gene",
            "category": row["phenotype_class"], "category_label": phen, "subject": gene, "gene": gene,
            "genotype": row["diplotype"],
            "statement": (f"{gene} {label}: {phen}." if not no_result else f"{gene} could not be typed from "
                                                                          "your data."),
            "evidence_level": evidence, "evidence_reasons": ev_reasons, "call_level": call_level,
            "call_reasons": call_reasons, "overall": overall, "sensitive": False,
            "sources": json.dumps([{"source": s, "version": v, "record": gene,
                                    "url": f"https://www.clinpgx.org/gene/{gene}" if s == "pharmcat" else ""}
                                   for s, v in versions.items() if v and (s != "imgt_hla" or gene.startswith("HLA"))]),
        })

    gene_level = {r["gene"]: r["call_level"] for r in genes}
    for r in drug_rows:
        if not r["matched"] or r["source"] not in ("CPIC", "DPWG"):
            continue
        levels = [gene_level.get(g, "Not callable") for g in r["genes"]] or ["Not callable"]
        call = min(levels, key=ev.CALL.index)
        evidence, ev_reasons = _cpic_evidence(r["classification"], r["source"])
        r_id = f"pgx:drug:{r['source']}:{r['drug']}:{r['population'] or ''}"
        claims.append({
            "claim_id": r_id, "section": "pgx", "group": "pgx_drug", "kind": "pgx_drug",
            "category": "action" if r["action"] else "standard",
            "category_label": "Guidance differs for you" if r["action"] else "Standard use",
            "subject": r["drug"], "gene": ", ".join(r["genes"]) or None, "genotype": "; ".join(r["phenotypes"]),
            "conditions": r["population"],
            "statement": f"{r['source']}: {r['recommendation']}", "evidence_level": evidence,
            "evidence_reasons": ev_reasons, "call_level": call,
            "call_reasons": [f"Depends on {g}: {gene_level.get(g, 'not typed')} call confidence" for g in r["genes"]],
            "overall": ev.overall(evidence, call), "sensitive": False,
            "sources": json.dumps([{"source": "pharmcat", "version": versions.get("pharmcat"), "record": r["drug"],
                                    "url": r["url"] or ""}]),
        })
        r["overall"] = claims[-1]["overall"]

    con = duckdb.connect()

    def save(rows: list[dict], dest: Path) -> None:
        tmp = dest.with_suffix(".json")
        tmp.write_text("\n".join(json.dumps(r, default=str) for r in rows))
        con.execute(f"COPY (SELECT * FROM read_json('{tmp}', format='newline_delimited', "
                    f"sample_size=-1)) TO '{dest}' (FORMAT parquet)")
        tmp.unlink()

    save(genes, paths["genes"])
    save(drug_rows, paths["drugs"])
    save([{k: c.get(k) for k in ("chrom38", "pos38", "rsid", "ref", "alts", "gene", "chrom37", "pos37", "gt",
                                  "status", "note")} for c in calls], paths["positions"])
    save(claims, paths["claims"])
    summary = {
        "pharmcat_version": report.get("pharmcatVersion"), "data_version": report.get("dataVersion"),
        "knowledge": versions, "positions": len(calls),
        "positions_missing": sum(1 for c in calls if c["gt"].startswith(".")),
        "positions_status": _count(c["status"] for c in calls),
        "outside_calls": {g: {"call": v["call"], "tool": v["tool"]} for g, v in oc.items()},
        "genes_called": sum(1 for g in genes if g["diplotype"]), "genes": len(genes),
        "drugs_action": sorted({r["drug"] for r in drug_rows if r["action"]}),
        "sample_sex": _sample_sex(ctx), "build": _build(ctx) or (ctx.inv.first("vcf").meta or {}).get("build"),
        "messages": [m.get("message") for m in report.get("messages") or []][:10],
    }
    write_json(paths["summary"], summary)
    return {"genes_called": summary["genes_called"], "drugs_action": len(summary["drugs_action"]),
            "positions_missing": summary["positions_missing"]}


def _count(it) -> dict:
    out: dict = {}
    for x in it:
        out[x] = out.get(x, 0) + 1
    return out


def _pgx_available(ctx: Context) -> str | None:
    st = Store(ctx.cfg)
    if not st.file("pharmcat", JAR) or not st.table("pharmcat", "positions"):
        return "needs the PharmCAT knowledge source (wgs knowledge refresh --source pharmcat)"
    if not _java():
        return "needs Java 17+ (provided by pixi)"
    if not variants.variants_parquet(ctx).exists():
        return "needs a VCF"
    if not (coverage.out_dir(ctx) / "callable.parquet").exists():
        return "needs the callable map (coverage stage) to tell reference from missing"
    if not reference.reference_fasta(ctx):
        return "needs a verified reference to check REF bases"
    vb = ((ctx.inv.first("vcf") or ctx.inv.first("gvcf")).meta or {}).get("build")
    if vb not in ("GRCh37", "GRCh38"):
        return f"PharmCAT positions are lifted to GRCh37/GRCh38 only (VCF is {vb})"
    return None


def _pgx_inputs(ctx: Context) -> list[Path]:
    st = Store(ctx.cfg)
    ins = [variants.variants_parquet(ctx), coverage.out_dir(ctx) / "callable.parquet", qc.out_path(ctx),
           st.file("pharmcat", JAR), st.table("pharmcat", "positions"), cyrius_out(ctx), hla_out(ctx)]
    ins += [f.path for t in ("cyrius", "t1k") for f in ctx.inv.of(f"external:{t}")]
    return [p for p in ins if p and p.exists()]


PGX = Stage(
    name="pgx", version="2", title="Medicines: PharmCAT diplotypes, phenotypes and guidelines", fn=_pgx,
    inputs=_pgx_inputs, outputs=lambda ctx: [out_paths(ctx)["summary"], out_paths(ctx)["claims"]],
    available=_pgx_available,
)
