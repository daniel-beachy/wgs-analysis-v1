"""Provider raw-genotype files (23andMe-style `rsid chrom pos genotype`) → Parquet."""

from __future__ import annotations

import zipfile
from pathlib import Path

import duckdb

from ..pipeline import Context, Stage


def out_path(ctx: Context) -> Path:
    p = ctx.work / "genotypes"
    p.mkdir(parents=True, exist_ok=True)
    return p / "raw_genotypes.parquet"


def _extract(ctx: Context, f) -> Path:
    if not f.meta.get("zipped"):
        return f.path
    dest = ctx.scratch / Path(f.meta["member"]).name
    with zipfile.ZipFile(f.path) as zf, zf.open(f.meta["member"]) as src, open(dest, "wb") as out:
        while chunk := src.read(8 << 20):
            out.write(chunk)
    return dest


def build(ctx: Context) -> dict:
    files = ctx.inv.of("raw_genotypes")
    con = duckdb.connect()
    parts, extracted = [], []
    for i, f in enumerate(sorted(files, key=lambda x: x.size)):
        path = _extract(ctx, f)
        if path != f.path:
            extracted.append(path)
        con.execute(f"""
            CREATE TEMP TABLE g{i} AS SELECT column0 AS rsid, column1 AS chrom_raw, column2::INTEGER AS pos,
                   column3 AS genotype, '{f.path.name}' AS source_file, {i} AS file_rank
            FROM read_csv('{path}', delim='\t', header=false, comment='#', quote='', escape='',
                          columns={{'column0':'VARCHAR','column1':'VARCHAR','column2':'VARCHAR','column3':'VARCHAR'}},
                          auto_detect=false)""")
        parts.append(f"SELECT * FROM g{i}")
    union = " UNION ALL ".join(parts)
    out = out_path(ctx)
    tmp = out.with_suffix(".parquet.tmp")
    # Smallest file is the provider's array-compatible panel; positions in it are flagged in_panel.
    con.execute(f"""
        COPY (
          WITH u AS ({union}),
          d AS (
            SELECT regexp_replace(regexp_replace(chrom_raw, '^chr', ''), '^(M|26)$', 'MT') AS chrom, pos,
                   arg_min(rsid, file_rank) AS rsid, arg_min(genotype, file_rank) AS genotype,
                   bool_or(file_rank = 0) AS in_panel, count(DISTINCT genotype) > 1 AS conflicting
            FROM u GROUP BY ALL)
          SELECT * FROM d ORDER BY chrom, pos
        ) TO '{tmp}' (FORMAT parquet, COMPRESSION zstd)""")
    tmp.replace(out)
    n, panel, conflicts = con.execute(
        f"SELECT count(*), count(*) FILTER (WHERE in_panel), count(*) FILTER (WHERE conflicting) FROM '{out}'"
    ).fetchone()
    for p in extracted:
        p.unlink(missing_ok=True)
    return {"positions": n, "panel_positions": panel, "conflicting_between_files": conflicts,
            "files": [f.path.name for f in files]}


STAGE = Stage(
    name="genotypes", version="1", title="Provider raw genotypes → Parquet", fn=build,
    inputs=lambda ctx: [f.path for f in ctx.inv.of("raw_genotypes")],
    outputs=lambda ctx: [out_path(ctx)],
    available=lambda ctx: None if ctx.inv.of("raw_genotypes") else "no raw genotype file found",
)
