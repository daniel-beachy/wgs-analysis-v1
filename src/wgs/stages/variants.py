"""VCF → columnar variant store (Parquet)."""

from __future__ import annotations

from pathlib import Path

import duckdb

from ..pipeline import CANONICAL, Context, Stage, run

QUERY_FMT = r"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/OLD\t[%GT\t%GQ\t%DP\t%AD\t%VAF]\n"
COLUMNS = {
    "chrom_raw": "VARCHAR", "pos": "BIGINT", "vid": "VARCHAR", "ref": "VARCHAR", "alt": "VARCHAR",
    "qual": "VARCHAR", "filter": "VARCHAR", "old": "VARCHAR", "gt": "VARCHAR", "gq": "VARCHAR",
    "dp": "VARCHAR", "ad": "VARCHAR", "vaf": "VARCHAR",
}


def vcf_file(ctx: Context) -> Path | None:
    f = ctx.inv.first("vcf") or ctx.inv.first("gvcf")
    return f.path if f else None


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "variants"
    p.mkdir(parents=True, exist_ok=True)
    return p


def variants_parquet(ctx: Context) -> Path:
    return out_dir(ctx) / "variants.parquet"


def vcf_index(ctx: Context) -> Path:
    """Index stored in the workspace so input folders stay read-only."""
    return out_dir(ctx) / "input.vcf.gz.tbi"


def vcf_with_index(ctx: Context) -> str:
    vcf = vcf_file(ctx)
    own = Path(str(vcf) + ".tbi")
    return str(vcf) if own.exists() else f"{vcf}##idx##{vcf_index(ctx)}"


def _sql_num(col: str, typ: str) -> str:
    return f"TRY_CAST(NULLIF(NULLIF({col}, '.'), '') AS {typ})"


def build(ctx: Context) -> dict:
    vcf = vcf_file(ctx)
    log = ctx.logs / "variants.log"
    if not Path(str(vcf) + ".tbi").exists():
        run(["bcftools", "index", "-f", "-t", "--threads", "2", "-o", str(vcf_index(ctx)), str(vcf)], log=log)
    tsv = ctx.scratch / "variants.tsv"
    run(f"bcftools norm -m -any --old-rec-tag OLD --threads 2 -Ou '{vcf}' | "
        f"bcftools query -f '{QUERY_FMT}' > '{tsv}'", shell=True, log=log)

    chrom_order = "CASE c " + " ".join(f"WHEN '{c}' THEN {i}" for i, c in enumerate(CANONICAL)) + " ELSE 99 END"
    out = variants_parquet(ctx)
    tmp = out.with_suffix(".parquet.tmp")
    con = duckdb.connect()
    con.execute(f"SET threads={ctx.threads}; SET preserve_insertion_order=false;")
    con.execute(f"""
        CREATE TEMP VIEW raw AS SELECT * FROM read_csv('{tsv}', delim='\t', header=false, quote='', escape='',
            columns={COLUMNS!r}, auto_detect=false);
    """)
    con.execute(f"""
        COPY (
          WITH t AS (
            SELECT *, regexp_replace(regexp_replace(chrom_raw, '^chr', ''), '^M$', 'MT') AS c,
                   string_split(ad, ',') AS ads FROM raw)
          SELECT
            c AS chrom, ({chrom_order})::UTINYINT AS chrom_order, pos::INTEGER AS pos,
            NULLIF(vid, '.') AS vcf_id, ref, alt,
            CASE WHEN length(ref)=1 AND length(alt)=1 THEN 'SNV'
                 WHEN length(ref)=length(alt) THEN 'MNV'
                 WHEN length(ref)<length(alt) THEN 'INS' ELSE 'DEL' END AS vtype,
            {_sql_num('qual', 'FLOAT')} AS qual, filter,
            gt,
            CASE WHEN gt IN ('./.', '.', './0', '0/.') THEN 'nocall'
                 WHEN gt IN ('0/0', '0', '0|0') THEN 'ref'
                 WHEN gt IN ('1', '2') THEN 'hemi'
                 WHEN split_part(replace(gt, '|', '/'), '/', 1) = split_part(replace(gt, '|', '/'), '/', 2)
                      THEN 'hom'
                 ELSE 'het' END AS zygosity,
            {_sql_num('gq', 'SMALLINT')} AS gq, {_sql_num('dp', 'SMALLINT')} AS dp,
            TRY_CAST(ads[1] AS SMALLINT) AS ad_ref, TRY_CAST(ads[2] AS SMALLINT) AS ad_alt,
            {_sql_num('vaf', 'FLOAT')} AS vaf,
            old IS NOT NULL AND old <> '.' AS multiallelic
          FROM t
          -- Splitting e.g. a 0/2 site yields a meaningless 0/0 row for allele 1; drop those.
          WHERE NOT (old IS NOT NULL AND old <> '.' AND gt IN ('0/0', '0|0') AND filter = 'PASS')
          ORDER BY chrom_order, chrom, pos, ref, alt
        ) TO '{tmp}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);
    """)
    tmp.replace(out)
    stats = con.execute(f"""
        SELECT count(*) AS records,
               count(*) FILTER (WHERE filter='PASS') AS pass,
               count(*) FILTER (WHERE filter='RefCall') AS refcall,
               count(*) FILTER (WHERE filter='NoCall') AS nocall
        FROM '{out}'""").fetchone()
    tsv.unlink(missing_ok=True)
    return {"records": stats[0], "pass": stats[1], "refcall": stats[2], "nocall": stats[3],
            "parquet_bytes": out.stat().st_size}


STAGE = Stage(
    name="variants", version="2", title="Variant store (VCF → Parquet)", fn=build,
    inputs=lambda ctx: [vcf_file(ctx)],
    outputs=lambda ctx: [variants_parquet(ctx)],
    available=lambda ctx: None if vcf_file(ctx) else "no VCF found",
)
