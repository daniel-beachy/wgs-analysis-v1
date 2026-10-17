"""Demo mode (Step 9c): a truth-set VCF with no reads runs through the pipeline (synthetic data only)."""

import json
import subprocess
from pathlib import Path

import duckdb
import pytest
from typer.testing import CliRunner

from wgs import demo
from wgs.cli import app
from wgs.config import load_config
from wgs.discovery import discover
from wgs.knowledge import Store, refresh
from wgs.pipeline import Context, Stage, run_stage
from wgs.stages import ancestry, coverage, qc, reference, variants

SEQ = "ACGTACGTAC" * 100  # 1,000 letters

TRUTH = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=platforms,Number=1,Type=Integer,Description="platforms">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="summed over datasets">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##contig=<ID=1,length=1000>
##contig=<ID=hs37d5,length=500>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
1\t11\t.\tA\tT\t50\tPASS\tplatforms=3\tGT:PS:DP:AD:GQ\t0/1:.:300:150,150:99
1\t101\t.\tA\tG\t50\tPASS\tplatforms=3\tGT:PS:DP:AD:GQ\t1/1:.:280:0,280:120
1\t901\t.\tA\tC\t50\tPASS\tplatforms=2\tGT:PS:DP:AD:GQ\t0/1:.:250:120,130:80
"""


def _bgzip(text: str, path: Path) -> Path:
    path.write_text(text)
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    gz = Path(str(path) + ".gz")
    subprocess.run(["tabix", "-f", "-p", "vcf", str(gz)], check=True)
    return gz


@pytest.fixture
def demo_cfg(tmp_path):
    giab = tmp_path / "giab"
    giab.mkdir()
    truth = _bgzip(TRUTH, giab / "truth.vcf")
    (giab / "confident.bed").write_text("1\t0\t500\n")  # the record at 901 lies outside
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "\n".join(SEQ[i:i + 60] for i in range(0, len(SEQ), 60)) + "\n")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    ws = tmp_path / "ws"
    (ws / "input").mkdir(parents=True)
    conf = tmp_path / "demo.toml"
    conf.write_text(f"""[paths]
workspace = "{ws}"
search = ["{ws / 'input'}"]
knowledge = "{tmp_path / 'main' / 'knowledge'}"
[sample]
id = "HG002"
[inputs]
reference = "{fasta}"
[truth]
label = "test truth set"
vcf = "{truth}"
regions = "{giab / 'confident.bed'}"
""")
    return load_config(str(conf))


def _ctx(cfg) -> Context:
    return Context(cfg=cfg, inv=discover(cfg), threads=2)


def test_prepare_keeps_confident_records_renamed_and_slimmed(demo_cfg):
    out = demo.prepare(demo_cfg)
    lines = subprocess.run(["bcftools", "view", "-H", str(out)], capture_output=True, text=True,
                           check=True).stdout.splitlines()
    assert [ln.split("\t")[:2] for ln in lines] == [["chr1", "11"], ["chr1", "101"]]
    assert lines[0].split("\t")[7:] == [".", "GT:GQ", "0/1:99"]  # summed depths would read as one run's reads
    assert demo.prepare(demo_cfg) == out  # up to date: no rewrite


def test_variants_stage_reads_a_vcf_without_depth_fields(demo_cfg):
    demo.prepare(demo_cfg)
    ctx = _ctx(demo_cfg)
    assert run_stage(ctx, variants.STAGE)["status"] == "ran"
    rows = duckdb.connect().execute(f"SELECT chrom, pos, zygosity, gq, dp, ad_alt, vaf FROM "
                                    f"'{variants.variants_parquet(ctx)}' ORDER BY pos").fetchall()
    assert rows == [("1", 11, "het", 99, None, None, None), ("1", 101, "hom", 120, None, None, None)]


def test_reference_without_reads_is_checked_against_the_vcf(demo_cfg):
    demo.prepare(demo_cfg)
    ctx = _ctx(demo_cfg)
    run_stage(ctx, variants.STAGE)
    res = run_stage(ctx, reference.STAGE)
    assert res["status"] == "ran" and res["vcf_ref_mismatches"] == 0
    assert reference.reference_fasta(ctx) == Path(demo_cfg.explicit["reference"][0])


def test_reference_without_reads_rejects_a_mismatching_fasta(demo_cfg, tmp_path):
    demo.prepare(demo_cfg)
    fasta = Path(demo_cfg.explicit["reference"][0])
    fasta.write_text(">chr1\n" + SEQ.replace("ACGTACGTAC", "TTTTTTTTTT", 2) + "\n")  # position 11 is no longer A
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    ctx = _ctx(demo_cfg)
    run_stage(ctx, variants.STAGE)
    with pytest.raises(RuntimeError, match="disagree"):
        run_stage(ctx, reference.STAGE)
    assert reference.reference_fasta(ctx) is None


def test_callable_map_comes_from_the_truth_regions(demo_cfg):
    ctx = _ctx(demo_cfg)
    res = run_stage(ctx, coverage.STAGE)
    assert res["status"] == "ran" and res["callable_bp"] == {"CALLABLE": 500}
    rows = duckdb.connect().execute(f"SELECT * FROM '{coverage.out_dir(ctx) / 'callable.parquet'}'").fetchall()
    assert rows == [("1", 0, 500, "CALLABLE")]


def test_borrowed_knowledge_store_is_read_only(demo_cfg):
    assert demo_cfg.knowledge_borrowed and demo_cfg.knowledge_dir.name == "knowledge"
    Store(demo_cfg)
    assert not demo_cfg.knowledge_dir.exists()  # never created in someone else's workspace
    with pytest.raises(RuntimeError, match="borrows"):
        refresh(demo_cfg, demo_cfg.workspace / "scratch")


def test_qc_runs_on_a_truth_set_without_reads(demo_cfg):
    demo.prepare(demo_cfg)
    ctx = _ctx(demo_cfg)
    for st in (variants.STAGE, coverage.STAGE):
        run_stage(ctx, st)
    assert run_stage(ctx, qc.STAGE)["status"] == "ran"  # no read fractions, and only the callable map
    out = json.loads(qc.out_path(ctx).read_text())
    vcf = out["sections"]["vcf"]
    assert vcf["het_vaf_median"] is None and vcf["hist"]["depth"] == []
    assert out["inferred_sex"]["call"] is None  # no chrX records: unknown, not XY


def test_demo_runs_from_the_sample_scratch_folder_and_leaves_out_pgs(demo_cfg, tmp_path, monkeypatch):
    import wgs.stages.registry as registry

    seen = []

    def stage(name):
        return Stage(name=name, version="1", title=name, fn=lambda ctx: seen.append((name, Path.cwd())) or {},
                     inputs=lambda ctx: [], outputs=lambda ctx: [])

    monkeypatch.setattr(registry, "STAGES", [stage("pgs"), stage("probe")])
    monkeypatch.setenv("WGS_SCRATCH", str(tmp_path / "scratch"))
    home = Path.cwd()
    res = CliRunner().invoke(app, ["--config", str(tmp_path / "demo.toml"), "demo"])
    assert res.exit_code == 0, res.output
    # DuckDB spills into the working folder under fixed names: each sample's run gets its own
    scratch = (tmp_path / "scratch" / "wgs-scratch" / "HG002").resolve()
    assert seen == [("probe", scratch)] and Path.cwd() == home
    assert "--pgs" in res.output  # the hours-long reference scoring only on request
    seen.clear()
    assert CliRunner().invoke(app, ["--config", str(tmp_path / "demo.toml"), "demo", "--pgs",
                                    "--force", "probe"]).exit_code == 0
    assert [name for name, _ in seen] == ["pgs", "probe"]


def test_ancestry_skips_a_truth_set_with_nothing_to_place(demo_cfg):
    demo.prepare(demo_cfg)
    ctx = _ctx(demo_cfg)
    run_stage(ctx, variants.STAGE)
    # no PCA, no reads and no chrM/chrY records: skip with a reason instead of writing tables without columns
    assert run_stage(ctx, ancestry.STAGE) == {
        "status": "skipped", "reason": "needs the PGS stage's ancestry PCA, or mtDNA / Y-chromosome calls (this VCF "
                                       "has none)"}
    assert not ancestry.out_components(ctx).exists()
