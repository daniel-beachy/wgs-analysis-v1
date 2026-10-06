"""Read-level coverage with mosdepth (one pass over the CRAM).

Outputs:
  * coverage_windows.parquet  — mean depth in 1 kb windows (genome-wide plots, CNV-ish eyeballing)
  * callable.parquet          — intervals binned by depth (0, 1-4, 5-9, 10-149, 150+). The VCF is not a gVCF, so
                                this is what lets us say "reference" (covered, no variant) vs "unknown" (not covered).
  * coverage.json             — summary + cumulative depth distribution
"""

from __future__ import annotations

from pathlib import Path

import duckdb

from ..pipeline import Context, Stage, run
from .reference import reference_fasta

QUANTIZE = "0:1:5:10:150:"
BINS = {"0": "NO_COVERAGE", "1": "LOW_1_4", "2": "LOW_5_9", "3": "CALLABLE", "4": "HIGH"}


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "coverage"
    p.mkdir(parents=True, exist_ok=True)
    return p


def outputs(ctx: Context) -> list[Path]:
    d = out_dir(ctx)
    return [d / "coverage_windows.parquet", d / "callable.parquet", d / "coverage.json"]


def build(ctx: Context) -> dict:
    cram = ctx.inv.first("cram") or ctx.inv.first("bam")
    fasta = reference_fasta(ctx)
    prefix = ctx.scratch / "mosdepth" / ctx.sample
    prefix.parent.mkdir(parents=True, exist_ok=True)
    env = {f"MOSDEPTH_Q{k}": v for k, v in BINS.items()}
    run(["mosdepth", "-t", "4", "-n", "--fast-mode", "--by", "1000", "--quantize", QUANTIZE,
         "--fasta", str(fasta), str(prefix), str(cram.path)], log=ctx.logs / "coverage.log", env=env)

    con = duckdb.connect()
    d = out_dir(ctx)
    norm = "regexp_replace(regexp_replace(column0, '^chr', ''), '^M$', 'MT')"
    con.execute(f"""COPY (SELECT {norm} AS chrom, column1::INTEGER AS start, column2::INTEGER AS "end",
                         column3::FLOAT AS depth
                  FROM read_csv('{prefix}.regions.bed.gz', delim='\t', header=false))
                  TO '{d / 'coverage_windows.parquet'}' (FORMAT parquet, COMPRESSION zstd)""")
    con.execute(f"""COPY (SELECT {norm} AS chrom, column1::INTEGER AS start, column2::INTEGER AS "end",
                         column3 AS state
                  FROM read_csv('{prefix}.quantized.bed.gz', delim='\t', header=false))
                  TO '{d / 'callable.parquet'}' (FORMAT parquet, COMPRESSION zstd)""")
    summary = con.execute(f"""SELECT * FROM read_csv('{prefix}.mosdepth.summary.txt', delim='\t', header=true)
                          """).fetchdf()
    dist = con.execute(f"""SELECT column0 AS chrom, column1::INTEGER AS depth, column2::DOUBLE AS frac
                       FROM read_csv('{prefix}.mosdepth.global.dist.txt', delim='\t', header=false)""").fetchdf()
    states = con.execute(f"""SELECT state, sum("end"-start) AS bp FROM '{d / 'callable.parquet'}'
                          WHERE chrom IN ({",".join(f"'{c}'" for c in [*map(str, range(1, 23)), 'X', 'Y'])})
                          GROUP BY state""").fetchall()
    per_chrom = {r["chrom"].removeprefix("chr").replace("M", "MT") if r["chrom"] != "total" else "total":
                 {"length": int(r["length"]), "mean": float(r["mean"]), "min": int(r["min"]), "max": int(r["max"])}
                 for _, r in summary.iterrows() if not str(r["chrom"]).endswith("_region")}
    total = dist[dist.chrom == "total"]
    cum = {int(x): float(total[total.depth == x].frac.iloc[0]) for x in (1, 5, 10, 15, 20, 30, 40, 50)
           if (total.depth == x).any()}
    auto = [per_chrom[str(i)]["mean"] for i in range(1, 23) if str(i) in per_chrom]
    auto_mean = sum(auto) / len(auto) if auto else 0
    result = {
        "autosomal_mean_depth": round(auto_mean, 2),
        "x_ratio": round(per_chrom.get("X", {}).get("mean", 0) / auto_mean, 3) if auto_mean else None,
        "y_ratio": round(per_chrom.get("Y", {}).get("mean", 0) / auto_mean, 3) if auto_mean else None,
        "mt_mean_depth": per_chrom.get("MT", {}).get("mean"),
        "fraction_at_least": cum,
        "callable_bp": {s: int(bp) for s, bp in states},
        "per_chrom": per_chrom,
    }
    from ..pipeline import write_json
    write_json(d / "coverage.json", result)
    return {k: v for k, v in result.items() if k != "per_chrom"}


def _available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "no CRAM/BAM found"
    if ctx.inv.first("cram") and not reference_fasta(ctx):
        return "CRAM reference not verified yet (reference stage)"
    return None


STAGE = Stage(
    name="coverage", version="1", title="Read coverage (mosdepth)", fn=build,
    inputs=lambda ctx: [f.path for f in (ctx.inv.of("cram") or ctx.inv.of("bam"))[:1]],
    outputs=outputs, available=_available,
)
