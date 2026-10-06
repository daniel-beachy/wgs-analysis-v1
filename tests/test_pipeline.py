"""End-to-end run of the cheap stages on a tiny synthetic genome (no real data)."""

import json
import subprocess
import zipfile
from pathlib import Path

import duckdb
import pytest

from wgs.config import load_config
from wgs.discovery import discover
from wgs.pipeline import Context, canon_chrom, run_stage
from wgs.stages import genotypes, qc, variants

HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=RefCall,Description="ref">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="VAF">
##contig=<ID=chr1,length=249250621>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrM,length=16569>
##DeepVariant_version=1.10.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH
"""
ROWS = [
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:AD:VAF\t0/1:40:30:15,15:0.5",
    "chr1\t200\t.\tC\tT\t60\tPASS\t.\tGT:GQ:DP:AD:VAF\t1/1:50:32:0,32:1",
    "chr1\t300\t.\tG\tA,T\t55\tPASS\t.\tGT:GQ:DP:AD:VAF\t0/2:45:30:15,0,15:0,0.5",
    "chr1\t400\t.\tAT\tA\t40\tPASS\t.\tGT:GQ:DP:AD:VAF\t0/1:30:28:14,14:0.5",
    "chr1\t500\t.\tC\tG\t5\tRefCall\t.\tGT:GQ:DP:AD:VAF\t0/0:10:30:28,2:0.07",
    "chrX\t5000000\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:AD:VAF\t1/1:50:15:0,15:1",
    "chrM\t73\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP:AD:VAF\t1/1:50:300:0,300:1",
]


@pytest.fixture
def ctx(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    plain = data / "calls.vcf"
    plain.write_text(HEADER + "\n".join(ROWS) + "\n")
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    rows = ["rs1\t1\t100\tAG", "rs2\t1\t200\tTT", "rs3\t1\t300\tGT", "rs4\t1\t500\tCC", "rs5\t1\t600\tDI"]
    with zipfile.ZipFile(data / "SYNTH.zip", "w") as zf:
        zf.writestr("SYNTH.txt", "# rsid\tchromosome\tposition\tgenotype\n" + "\n".join(rows) + "\n")
    cfg = load_config(data=[str(data)], workspace=str(tmp_path / "ws"))
    return Context(cfg=cfg, inv=discover(cfg), threads=2)


def test_canon_chrom():
    assert [canon_chrom(c) for c in ("chr1", "chrM", "MT", "X", "chrY")] == ["1", "MT", "MT", "X", "Y"]


def test_variants_stage_splits_and_drops_artifacts(ctx):
    assert run_stage(ctx, variants.STAGE)["status"] == "ran"
    con = duckdb.connect()
    rows = con.execute(f"SELECT chrom, pos, alt, vtype, zygosity, filter FROM '{variants.variants_parquet(ctx)}'"
                       " ORDER BY chrom_order, pos").fetchall()
    # the 0/2 site keeps only its real allele (T); the split 0/0 artefact for A is dropped
    assert ("1", 300, "T", "SNV", "het", "PASS") in rows
    assert not any(r[1] == 300 and r[2] == "A" for r in rows)
    assert ("1", 400, "A", "DEL", "het", "PASS") in rows
    assert rows[-1][0] == "MT"
    # second run is a no-op thanks to the stamp
    assert run_stage(ctx, variants.STAGE)["status"] == "cached"


def test_qc_concordance_and_sex(ctx):
    for st in (variants.STAGE, genotypes.STAGE, qc.STAGE):
        assert run_stage(ctx, st)["status"] == "ran"
    out = json.loads(Path(qc.STAGE.outputs(ctx)[0]).read_text())
    c = out["sections"]["concordance"]
    # rs1 het, rs2 hom-alt, rs3 multiallelic het, rs4 RefCall → reference via VCF REF; DI is not an SNV
    assert c["discordant"] == 0 and c["concordant_variant"] == 3 and c["concordant_reference"] == 1
    assert out["inferred_sex"]["call"] == "XY"
