"""A second variant caller on the owner's reads (ADR-023).

The VCF comes from one variant caller (tellmeGen's, DeepVariant). bcftools, a different and simpler caller, is run on
the same reads in the ACMG genes' exons and at the medicine-defining positions (chromosomes 1-22: the regions the
HG002 benchmark measures), and RTG vcfeval compares the two sets of calls, matching variants written in different
but equivalent ways. Every variant only bcftools finds is re-counted from the reads (ADR-022), checked against
GIAB's map of hard-to-sequence regions, looked up in ClinVar, and looked up in gnomAD and 1000 Genomes (a change gnomAD
flags as failing its own quality checks usually comes from reads of a look-alike copy). The HG002 benchmark says how
often such a variant is real. Nothing here changes a finding: a variant the VCF lacks has no annotation or grade yet.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import duckdb
import pandas as pd

from ..knowledge import Store
from ..pipeline import Context, Stage, run, write_json
from . import benchmark as bm
from . import recount as rc
from . import reference, variants

VERDICTS = ["reads_support", "reads_flagged", "unclear", "reads_disagree"]
NEARBY = 10   # letters: a confident VCF call this close is often the same change written another way


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "secondcall"
    p.mkdir(parents=True, exist_ok=True)
    return p


def out_table(ctx: Context) -> Path:
    return out_dir(ctx) / "second_caller.parquet"


def out_summary(ctx: Context) -> Path:
    return out_dir(ctx) / "secondcall.json"


def verdict(r: dict) -> str:
    """What the reads say about a variant only bcftools called."""
    clear = (r["reads"] or 0) >= rc.MIN_READS and (r["lod"] or 0) >= rc.MIN_LOD
    if clear and (r["observed"] or 0) >= 1:
        return "reads_flagged" if len(r["flags"]) else "reads_support"
    if clear and r["observed"] == 0:
        return "reads_disagree"
    return "unclear"


def _vcf_prefix(vcf: str) -> str:
    import pysam

    with pysam.VariantFile(vcf) as vf:
        return "chr" if "chr1" in vf.header.contigs else ""


def _deepvariant_records(ctx: Context, cand: pd.DataFrame) -> pd.DataFrame:
    """The VCF's own record at each candidate position, if any (RefCall = the caller looked and saw no variant), and
    the distance to its nearest confident call within NEARBY letters (often the same change written differently)."""
    con = duckdb.connect()
    con.register("c", cand[["chrom", "pos"]].drop_duplicates())
    return con.execute(f"""SELECT c.chrom, c.pos,
          first(v.filter) FILTER (WHERE v.pos = c.pos) AS vcf_filter,
          first(v.gt) FILTER (WHERE v.pos = c.pos) AS vcf_gt,
          first(v.gq) FILTER (WHERE v.pos = c.pos) AS vcf_gq, first(v.dp) FILTER (WHERE v.pos = c.pos) AS vcf_dp,
          min(abs(v.pos - c.pos)) FILTER (WHERE v.pos <> c.pos AND v.filter = 'PASS') AS vcf_pass_nearby
        FROM c LEFT JOIN '{variants.variants_parquet(ctx)}' v
          ON v.chrom = c.chrom AND v.pos BETWEEN c.pos - {NEARBY} AND c.pos + {NEARBY}
        GROUP BY 1, 2""").df()


def _annotate(st: Store, cand: pd.DataFrame) -> pd.DataFrame:
    """Gene (Ensembl) and ClinVar record of each candidate."""
    con = duckdb.connect()
    con.register("c", cand[["chrom", "pos", "ref", "alt"]].drop_duplicates())
    genes = con.execute(f"""SELECT c.chrom, c.pos, string_agg(DISTINCT g.symbol, ', ' ORDER BY g.symbol) AS gene
        FROM c JOIN '{st.table("ensembl", "genes")}' g ON g.chrom = c.chrom AND g.start <= c.pos AND g."end" >= c.pos
        WHERE g.biotype = 'protein_coding' GROUP BY 1, 2""").df()
    cv = con.execute(f"""SELECT c.chrom, c.pos, c.ref, c.alt, first(v.significance) AS clinvar_significance,
        first(v.stars) AS clinvar_stars, first(v.variation_id) AS clinvar_id, first(v.sig_class) AS clinvar_class
        FROM c JOIN '{st.table("clinvar", "variants")}' v USING (chrom, pos, ref, alt) GROUP BY ALL""").df()
    return (cand.merge(genes, on=["chrom", "pos"], how="left")
                .merge(cv, on=["chrom", "pos", "ref", "alt"], how="left"))


def _secondcall(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    g = bm.giab_files(st)
    fasta = reference.reference_fasta(ctx)
    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    log = ctx.logs / "secondcall.log"
    work = ctx.scratch / "secondcall"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    iv = duckdb.connect().execute(f"SELECT chrom, start, \"end\" FROM '{g['regions']}' ORDER BY 1, 2").fetchall()
    aln_prefix = "chr" if "chr1" in reference.cram_contigs(ctx) else ""
    vcf = variants.vcf_with_index(ctx)
    vcf_prefix = _vcf_prefix(vcf)
    if aln_prefix != vcf_prefix:
        raise RuntimeError("the CRAM and the VCF name chromosomes differently; compare by hand")
    bed = bm.bed_file(iv, work / "regions.bed", aln_prefix)
    calls = bm.bcftools_calls(aln.path, fasta, bed, work / "bcftools.vcf.gz", log, ctx.threads)
    dv = work / "vcf.vcf.gz"
    run(f"bcftools view -f PASS,. -R '{bed}' '{vcf}' -Oz -o '{dv}'", shell=True, log=log)
    run(["bcftools", "index", "-f", "-t", str(dv)], log=log)
    ev_out = bm.vcfeval(dv, calls, bm.reference_sdf(ctx, fasta), bed, work / "vcfeval", log, ctx.threads,
                        squash=True)
    recs = bm.vcfeval_records(ev_out)
    cand = recs[recs["label"] == "fp"].drop(columns=["label"]).reset_index(drop=True)
    if len(cand):
        counted = pd.DataFrame(bm.recount_calls(aln.path, fasta, cand, lambda ch: f"{aln_prefix}{ch}"))
        cand = pd.concat([cand, counted], axis=1)
        cand["verdict"] = [verdict(r) for r in cand.to_dict("records")]
        cand["difficult"] = rc.difficult_reasons(cand, g["difficult"])
        cand = _annotate(st, cand).merge(_deepvariant_records(ctx, cand), on=["chrom", "pos"], how="left")
        cand = bm.population(st, cand, work, log)
    cand.to_parquet(out_table(ctx), index=False)
    n = recs["label"].value_counts().to_dict()
    summary = {
        "regions": len(iv), "bases": sum(e - s for _, s, e in iv),
        "tool": bm.tool_version(["bcftools", "--version"]), "scorer": bm.tool_version(["rtg", "version"]),
        "min_qual": bm.MIN_QUAL, "matched": n.get("tp", 0), "only_bcftools": n.get("fp", 0),
        "only_vcf": n.get("fn", 0),
        "only_bcftools_by_kind": {k: int(v) for k, v in cand["kind"].value_counts().items()} if len(cand) else {},
        "verdicts": {v: int((cand["verdict"] == v).sum()) for v in VERDICTS} if len(cand) else {},
        "population": {v: {p: int(((cand["verdict"] == v) & (cand["population"] == p)).sum())
                               for p in bm.POPULATION} for v in VERDICTS} if len(cand) else {},
        "only_vcf_by_kind": {k: int(v) for k, v in recs[recs["label"] == "fn"]["kind"].value_counts().items()},
        "clinvar_pathogenic": int(cand["clinvar_class"].isin(["pathogenic", "likely_pathogenic"]).sum())
        if len(cand) else 0,
    }
    write_json(out_summary(ctx), summary)
    return {k: summary[k] for k in ("matched", "only_bcftools", "only_vcf")}


def _inputs(ctx: Context) -> list[Path]:
    st = Store(ctx.cfg)
    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    g = bm.giab_files(st)
    return [p for p in (aln.path if aln else None, variants.vcf_file(ctx), g["regions"], g["difficult"],
                        st.table("clinvar", "variants"), st.table("gnomad", "fields"), st.table("1000genomes", "af"),
                        reference.reference_fasta(ctx)) if p]


def _available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "no CRAM/BAM found (a second caller needs the reads)"
    if not variants.vcf_file(ctx):
        return "no VCF to compare with"
    if not reference.reference_fasta(ctx):
        return "CRAM reference not verified yet (reference stage)"
    g = bm.giab_files(Store(ctx.cfg))
    if not g["regions"]:
        return "Genome in a Bottle regions missing: run `wgs knowledge refresh --source giab`"
    if not shutil.which("rtg"):
        return "RTG Tools not installed (pixi environment 'default' provides it)"
    return None


STAGE = Stage(
    name="secondcall", version="2", title="Second variant caller (bcftools) on your reads", fn=_secondcall,
    inputs=_inputs, outputs=lambda ctx: [out_table(ctx), out_summary(ctx)], available=_available,
    params=lambda ctx: {"min_qual": bm.MIN_QUAL, "nearby": NEARBY, "recount": rc.SETTINGS},
)
