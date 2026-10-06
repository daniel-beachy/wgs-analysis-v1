import gzip
import zipfile
from pathlib import Path

from wgs.config import load_config
from wgs.discovery import discover
from wgs.modules import evaluate

VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
##contig=<ID=chrM,length=16569>
##DeepVariant_version=1.10.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLEX
chr1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1
"""


def make_vcf(path: Path, sample="SAMPLEX", build_len=249250621):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = VCF_HEADER.replace("SAMPLEX", sample).replace("249250621", str(build_len))
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def make_raw(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = "\n".join(f"rs{i}\t1\t{1000 + i}\tAG" for i in range(10))
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("raw.txt", "# rsid\tchromosome\tposition\tgenotype\n" + rows + "\n")


def test_discovers_by_content_in_any_layout(tmp_path):
    make_vcf(tmp_path / "weird" / "nested" / "calls.vcf.gz")
    make_raw(tmp_path / "elsewhere" / "geno.zip")
    (tmp_path / "elsewhere" / "._geno.zip").write_bytes(b"\x00\x05\x16\x07")  # AppleDouble noise
    with gzip.open(tmp_path / "r_1.fq.gz", "wt") as fh:
        fh.write("@read1\nACGT\n+\nIIII\n")
    cfg = load_config(data=[str(tmp_path)], workspace=str(tmp_path / "ws"))
    inv = discover(cfg)
    kinds = sorted(f.kind for f in inv.files)
    assert kinds == ["fastq", "raw_genotypes", "vcf"]
    vcf = inv.first("vcf")
    assert vcf.meta["build"] == "GRCh37" and vcf.meta["samples"] == ["SAMPLEX"]
    assert inv.sample == "SAMPLEX"
    assert inv.first("fastq").meta["mate"] == 1


def test_graceful_degradation(tmp_path):
    make_raw(tmp_path / "geno.zip")
    cfg = load_config(data=[str(tmp_path)], workspace=str(tmp_path / "ws"))
    rows = {r["id"]: r["status"] for r in evaluate(discover(cfg), cfg)}
    assert rows["traits"] == "partial"  # single-SNP traits only, no polygenic scores
    assert rows["ancestry"] == "ready"
    assert rows["health"] == "unavailable"
    assert rows["pgx"] == "unavailable"


def test_empty_folder_is_not_an_error(tmp_path):
    cfg = load_config(data=[str(tmp_path)], workspace=str(tmp_path / "ws"))
    inv = discover(cfg)
    assert inv.files == [] and inv.sample is None
    assert all(r["status"] == "unavailable" for r in evaluate(inv, cfg))


def test_multiple_samples_select(tmp_path):
    make_vcf(tmp_path / "a" / "HG002.vcf.gz", sample="HG002", build_len=248956422)
    make_vcf(tmp_path / "b" / "ME.vcf.gz", sample="ME")
    cfg = load_config(data=[str(tmp_path)], workspace=str(tmp_path / "ws"), sample="HG002")
    inv = discover(cfg)
    assert [f.meta["samples"] for f in inv.of("vcf")] == [["HG002"]]
    assert inv.first("vcf").meta["build"] == "GRCh38"


def test_config_file_relative_paths_and_host_override(tmp_path, monkeypatch):
    make_vcf(tmp_path / "drive" / "data" / "x.vcf.gz")
    cfg_file = tmp_path / "drive" / "repo" / "wgs.local.toml"
    cfg_file.parent.mkdir(parents=True)
    cfg_file.write_text('[paths]\nworkspace = "../out"\nsearch = ["../nowhere"]\n'
                        '[host."testhost"]\nsearch = ["../data"]\n')
    monkeypatch.setattr("socket.gethostname", lambda: "testhost.local")
    cfg = load_config(config_path=cfg_file)
    assert cfg.workspace == (tmp_path / "drive" / "out").resolve()
    assert cfg.search == [(tmp_path / "drive" / "data").resolve()]
    assert discover(cfg).first("vcf") is not None


def test_zip_of_pdfs_is_a_report_bundle(tmp_path):
    with zipfile.ZipFile(tmp_path / "reports.zip", "w") as zf:
        zf.writestr("a_report.pdf", b"%PDF-1.4 fake")
        zf.writestr("b_report.pdf", b"%PDF-1.4 fake")
    make_raw(tmp_path / "geno.zip")
    cfg = load_config(data=[str(tmp_path)], workspace=str(tmp_path / "ws"))
    kinds = {f.path.name: (f.kind, f.meta.get("pdfs")) for f in discover(cfg).files}
    assert kinds["reports.zip"] == ("report", 2)
    assert kinds["geno.zip"][0] == "raw_genotypes"
