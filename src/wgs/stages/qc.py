"""QC summary: VCF metrics, sex check, provider concordance, plus read/coverage stats when available.

Each metric is compared against an expected range for 30x short-read WGS called with DeepVariant, so the dashboard
can show pass / warn with an explanation rather than bare numbers. Ranges come from published population-scale WGS
(e.g. gnomAD / 1000 Genomes high-coverage callsets) and are deliberately generous; they flag gross problems, not
subtle ones.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import duckdb

from ..pipeline import Context, Stage, write_json
from . import coverage, genotypes, reads, reference, variants

# GRCh37 chrX pseudo-autosomal regions (1-based)
PAR_X = [(60001, 2699520), (154931044, 155260560)]

EXPECTED = {
    "snv_count": (3.3e6, 4.6e6, "≈3.5–4.4 M SNVs per genome (higher in African ancestry)"),
    "indel_count": (5.5e5, 1.2e6, "≈0.6–1.0 M small indels with DeepVariant at 30x"),
    "titv": (1.95, 2.15, "Genome-wide Ti/Tv ≈ 2.0–2.1; lower values suggest false-positive SNVs"),
    "het_hom_ratio": (1.2, 2.2, "Het/hom-alt ≈ 1.3–1.7 (non-African) or ≈ 2 (African)"),
    "het_vaf_median": (0.45, 0.55, "Heterozygous calls should centre on 50% allele fraction"),
    "autosomal_mean_depth": (25, 60, "Ordered as 30x"),
    "fraction_20x": (0.85, 1.0, "≥ 85% of the genome at ≥ 20x is typical for 30x WGS"),
    "duplicate_pct": (0, 15, "PCR-free DNBSEQ libraries usually show < 10% duplicates"),
    "mapped_pct": (97, 100.1, "≥ 97% of reads should map"),
    "error_rate": (0, 0.01, "Mismatch rate against the reference (includes true variants)"),
    "concordance_pct": (99.0, 100.1, "Agreement between the VCF and the provider's genotype file"),
}


def out_path(ctx: Context) -> Path:
    return ctx.work / "qc" / "qc.json"


def _assess(metrics: dict) -> list[dict]:
    rows = []
    for key, (lo, hi, why) in EXPECTED.items():
        v = metrics.get(key)
        if v is None:
            continue
        rows.append({"metric": key, "value": v, "low": lo, "high": hi, "status": "pass" if lo <= v <= hi else "warn",
                     "explanation": why})
    return rows


def _vcf_metrics(con, vp: Path) -> dict:
    m = {}
    q = con.execute(f"""
      SELECT count(*) FILTER (WHERE vtype='SNV') snv, count(*) FILTER (WHERE vtype IN ('INS','DEL')) indel,
             count(*) FILTER (WHERE vtype='INS') ins, count(*) FILTER (WHERE vtype='DEL') del,
             count(*) FILTER (WHERE vtype='MNV') mnv,
             count(*) FILTER (WHERE vtype='SNV' AND zygosity='het') snv_het,
             count(*) FILTER (WHERE vtype='SNV' AND zygosity='hom') snv_hom,
             count(*) FILTER (WHERE vtype='SNV' AND (ref||alt) IN ('AG','GA','CT','TC')) ts,
             count(*) FILTER (WHERE multiallelic) multi,
             median(vaf) FILTER (WHERE zygosity='het' AND vtype='SNV' AND chrom NOT IN ('X','Y','MT')) vaf_med,
             quantile_cont(vaf, 0.25) FILTER (WHERE zygosity='het' AND vtype='SNV' AND chrom NOT IN ('X','Y','MT')) q1,
             quantile_cont(vaf, 0.75) FILTER (WHERE zygosity='het' AND vtype='SNV' AND chrom NOT IN ('X','Y','MT')) q3,
             median(dp) dp_med, median(gq) gq_med
      FROM '{vp}' WHERE filter='PASS' AND zygosity IN ('het','hom','hemi')""").fetchone()
    (snv, indel, ins, dele, mnv, het, hom, ts, multi, vaf_med, q1, q3, dp_med, gq_med) = q
    m.update(snv_count=snv, indel_count=indel, insertion_count=ins, deletion_count=dele, mnv_count=mnv,
             multiallelic_rows=multi, titv=round(ts / max(snv - ts, 1), 3), het_hom_ratio=round(het / max(hom, 1), 3),
             het_vaf_median=round(vaf_med, 3), het_vaf_iqr=[round(q1, 3), round(q3, 3)], depth_median=dp_med,
             gq_median=gq_med)
    filters = dict(con.execute(f"SELECT filter, count(*) FROM '{vp}' GROUP BY 1").fetchall())
    m["records_by_filter"] = filters

    par = " OR ".join(f"pos BETWEEN {a} AND {b}" for a, b in PAR_X)
    xh, xo, y = con.execute(f"""
      SELECT count(*) FILTER (WHERE chrom='X' AND NOT ({par}) AND zygosity='het'),
             count(*) FILTER (WHERE chrom='X' AND NOT ({par}) AND zygosity IN ('hom','hemi')),
             count(*) FILTER (WHERE chrom='Y')
      FROM '{vp}' WHERE filter='PASS' AND vtype='SNV'""").fetchone()
    m["chrX_nonpar_het_fraction"] = round(xh / max(xh + xo, 1), 4)
    m["chrY_pass_snvs"] = y

    m["per_chrom"] = [dict(zip(("chrom", "snv", "indel", "het", "hom"), r, strict=True)) for r in con.execute(f"""
      SELECT chrom, count(*) FILTER (WHERE vtype='SNV'), count(*) FILTER (WHERE vtype IN ('INS','DEL')),
             count(*) FILTER (WHERE zygosity='het'), count(*) FILTER (WHERE zygosity IN ('hom','hemi'))
      FROM '{vp}' WHERE filter='PASS' AND zygosity IN ('het','hom','hemi')
      GROUP BY chrom, chrom_order ORDER BY chrom_order""").fetchall()]
    m["hist"] = {
        "het_vaf": con.execute(f"""SELECT round(floor(vaf*50)/50, 2) b, count(*) FROM '{vp}'
            WHERE filter='PASS' AND zygosity='het' AND vtype='SNV' AND chrom NOT IN ('X','Y','MT')
            GROUP BY b ORDER BY b""").fetchall(),
        "depth": con.execute(f"""SELECT least(dp, 100) b, count(*) FROM '{vp}' WHERE filter='PASS'
            GROUP BY b ORDER BY b""").fetchall(),
        "gq": con.execute(f"""SELECT gq b, count(*) FROM '{vp}' WHERE filter='PASS'
            GROUP BY b ORDER BY b""").fetchall(),
        "indel_length": con.execute(f"""SELECT greatest(least(length(alt)-length(ref), 20), -20) b, count(*)
            FROM '{vp}' WHERE filter='PASS' AND vtype IN ('INS','DEL') AND zygosity <> 'ref'
            GROUP BY b ORDER BY b""").fetchall(),
        "substitutions": con.execute(f"""SELECT ref||'>'||alt b, count(*) FROM '{vp}'
            WHERE filter='PASS' AND vtype='SNV' AND zygosity <> 'ref' GROUP BY b ORDER BY b""").fetchall(),
    }
    return m


def _concordance(ctx: Context, con, vp: Path, gp: Path) -> dict:
    """Compare provider genotypes to the VCF. Hom-ref expectations need the reference base."""
    df = con.execute(f"""
      WITH v AS (
        SELECT chrom, pos, any_value(ref) AS "ref",
               list(alt ORDER BY alt) FILTER (WHERE filter='PASS' AND zygosity IN ('het','hom','hemi')) alts,
               list(zygosity ORDER BY alt) FILTER (WHERE filter='PASS' AND zygosity IN ('het','hom','hemi')) zyg,
               bool_or(filter='NoCall') nocall, bool_or(filter='RefCall') refcall
        FROM '{vp}' WHERE vtype='SNV' GROUP BY chrom, pos)
      SELECT g.chrom, g.pos, g.genotype, g.in_panel, v."ref", v.alts, v.zyg, v.nocall, v.refcall
      FROM '{gp}' g LEFT JOIN v USING (chrom, pos)
      WHERE regexp_full_match(g.genotype, '[ACGT]{{2}}')
      ORDER BY g.chrom, g.pos""").fetchall()
    cols = ("chrom", "pos", "genotype", "in_panel", "ref", "alts", "zyg", "nocall", "refcall")
    fasta = reference.reference_fasta(ctx)
    fa = None
    if fasta:
        import pysam
        fa = pysam.FastaFile(str(fasta))
    cache: dict[str, str] = {}

    def ref_base(chrom: str, pos: int) -> str | None:
        if fa is None:
            return None
        name = "chrM" if chrom == "MT" else (f"chr{chrom}" if f"chr{chrom}" in fa.references else chrom)
        if name not in cache:
            cache.clear()
            cache[name] = fa.fetch(name).upper() if name in fa.references else ""
        seq = cache[name]
        return seq[pos - 1] if 0 < pos <= len(seq) else None

    counts = {"concordant_variant": 0, "concordant_reference": 0, "discordant": 0, "vcf_nocall": 0,
              "unverifiable": 0}
    by_panel = {True: [0, 0], False: [0, 0]}
    examples = []
    for row in df:
        r = SimpleNamespace(**dict(zip(cols, row, strict=True)))
        raw = "".join(sorted(r.genotype))
        alts = list(r.alts or [])
        if alts:
            zyg = list(r.zyg)
            exp = (r.ref + alts[0] if zyg[0] == "het" else alts[0] * 2) if len(alts) == 1 else alts[0] + alts[1]
            kind = "concordant_variant"
        elif r.nocall and not r.refcall:
            counts["vcf_nocall"] += 1
            continue
        else:
            base = r.ref or ref_base(r.chrom, int(r.pos))
            if not base:
                counts["unverifiable"] += 1
                continue
            exp = base * 2
            kind = "concordant_reference"
        exp = "".join(sorted(exp))
        ok = exp == raw
        counts[kind if ok else "discordant"] += 1
        by_panel[bool(r.in_panel)][0 if ok else 1] += 1
        if not ok and len(examples) < 25:
            examples.append({"chrom": r.chrom, "pos": int(r.pos), "provider": r.genotype, "vcf_expected": exp,
                             "vcf_alts": alts})
    if fa:
        fa.close()
    compared = counts["concordant_variant"] + counts["concordant_reference"] + counts["discordant"]
    return {
        **counts,
        "compared": compared,
        "concordance_pct": round(100 * (compared - counts["discordant"]) / max(compared, 1), 4),
        "panel_concordance_pct": round(100 * by_panel[True][0] / max(sum(by_panel[True]), 1), 4),
        "discordant_examples": examples,
    }


def build(ctx: Context) -> dict:
    con = duckdb.connect()
    vp = variants.variants_parquet(ctx)
    out: dict = {"sample": ctx.sample, "sections": {}}
    metrics: dict = {}
    if vp.exists():
        vm = _vcf_metrics(con, vp)
        out["sections"]["vcf"] = vm
        metrics.update({k: v for k, v in vm.items() if not isinstance(v, (dict, list))})
    gp = genotypes.out_path(ctx)
    if vp.exists() and gp.exists():
        conc = _concordance(ctx, con, vp, gp)
        out["sections"]["concordance"] = conc
        metrics["concordance_pct"] = conc["concordance_pct"]
    for name, path in (("coverage", coverage.outputs(ctx)[2]), ("alignment", reads.aln_out(ctx)),
                       ("fastq", reads.fq_out(ctx))):
        if path.exists():
            data = json.loads(path.read_text())
            out["sections"][name] = data
            if name == "coverage":
                metrics["autosomal_mean_depth"] = data["autosomal_mean_depth"]
                metrics["fraction_20x"] = data["fraction_at_least"].get("20")
            if name == "alignment":
                for k in ("duplicate_pct", "mapped_pct", "error_rate"):
                    metrics[k] = data.get(k)
    rf = reference.ready_file(ctx)
    if rf.exists():
        out["sections"]["reference"] = {k: v for k, v in json.loads(rf.read_text()).items() if k != "md5"}

    # Genetic sex: chrX heterozygosity is decisive; coverage ratios confirm when available.
    xh = metrics.get("chrX_nonpar_het_fraction")
    cov = out["sections"].get("coverage", {})
    sex = None
    if xh is not None:
        sex = "XY" if xh < 0.1 else "XX"
        if cov.get("y_ratio") is not None and sex == "XY" and cov["y_ratio"] < 0.1:
            sex = "XY? (low chrY coverage)"
    out["inferred_sex"] = {"call": sex, "chrX_nonpar_het_fraction": xh, "x_ratio": cov.get("x_ratio"),
                           "y_ratio": cov.get("y_ratio"), "chrY_pass_snvs": metrics.get("chrY_pass_snvs")}
    out["metrics"] = metrics
    out["assessment"] = _assess(metrics)
    write_json(out_path(ctx), out)
    warns = [a["metric"] for a in out["assessment"] if a["status"] == "warn"]
    return {"inferred_sex": sex, "metrics_checked": len(out["assessment"]), "warnings": warns}


def _inputs(ctx: Context) -> list[Path]:
    cands = [variants.variants_parquet(ctx), genotypes.out_path(ctx), coverage.outputs(ctx)[2],
             reads.aln_out(ctx), reads.fq_out(ctx), reference.ready_file(ctx)]
    return [p for p in cands if p.exists()]


STAGE = Stage(
    name="qc", version="1", title="QC summary", fn=build, inputs=_inputs, outputs=lambda ctx: [out_path(ctx)],
    available=lambda ctx: None if _inputs(ctx) else "nothing to summarise yet",
)
