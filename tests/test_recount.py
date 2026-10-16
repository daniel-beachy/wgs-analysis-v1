"""Read re-count (ADR-022): counting, the genotype model, verdicts and how they lower a call."""

from __future__ import annotations

import random

import pandas as pd
import pysam
import pytest

from wgs import evidence as ev
from wgs.stages import recount as rc

UNIT = "GATGATTCTGTCCCAGGAGC"           # a 20-letter repeat unit, two copies in the reference


def _reference() -> str:
    rng = random.Random(7)
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    return flank(400) + UNIT + UNIT + flank(800)                      # repeat at [400, 440)


def _write(tmp_path, ref: str, reads: list[tuple[int, str, str]]):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fa))
    bam = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(ref)}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i, (start, seq, cigar) in enumerate(sorted(reads)):
            a = pysam.AlignedSegment()
            a.query_name, a.query_sequence, a.flag = f"r{i}", seq, 16 if i % 2 else 0
            a.reference_id, a.reference_start, a.mapping_quality, a.cigarstring = 0, start, 60, cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            out.write(a)
    pysam.index(str(bam))
    return pysam.AlignmentFile(str(bam), "rb"), pysam.FastaFile(str(fa))


def test_softclipped_insertion_is_counted(tmp_path):
    """A 20-letter insertion in a tandem repeat: the aligner soft-clips the reads that carry it. Counting only reads
    aligned across the site missed every one of them and called a real heterozygous insertion absent."""
    ref = _reference()
    alt_hap = ref[:420] + UNIT + ref[420:]                            # three copies
    reads = [(s, ref[s:s + 150], "150M") for s in range(300, 380, 8)]
    reads += [(s, alt_hap[s:s + 150], f"{440 - s}M{150 - (440 - s)}S") for s in range(300, 380, 8)]
    af, fa = _write(tmp_path, ref, reads)
    t = rc.Target("k", "1", 420, ref[419], ref[419] + UNIT, expected=1)
    c = rc.count_site(af, fa, "chr1", t)
    assert c.alt >= 8 and c.ref >= 8
    assert rc.classify(c, t)["status"] == "confirmed"


def test_snv_counts_and_strand(tmp_path):
    ref = _reference()
    p0 = 700
    alt = "A" if ref[p0] != "A" else "C"
    reads = [(s, ref[s:s + 150], "150M") for s in range(600, 690, 6)]
    reads += [(s, ref[s:p0] + alt + ref[p0 + 1:s + 150], "150M") for s in range(603, 693, 6)]
    af, fa = _write(tmp_path, ref, reads)
    t = rc.Target("k", "1", p0 + 1, ref[p0], alt, expected=2)
    c = rc.count_site(af, fa, "chr1", t)
    assert (c.ref, c.alt) == (15, 15)
    v = rc.classify(c, t)
    assert v["observed"] == 1 and v["status"] == "disagrees"


def test_indel_span_slides_through_repeat():
    seq = "CGAAAAAT"
    assert rc.indel_span(seq, 1, "GA", "G") == (2, 7)      # one A deleted: any of the five A's
    assert rc.indel_span(seq, 1, "G", "GA") == (2, 7)      # one A inserted: anywhere in the run


@pytest.mark.parametrize("k,n,ploidy,copies", [(0, 20, 2, 0), (10, 20, 2, 1), (20, 20, 2, 2), (10, 10, 1, 1),
                                               (1, 30, 2, 0)])
def test_genotype_model(k, n, ploidy, copies):
    assert rc.call(k, n, ploidy)[0] == copies


def test_fisher_two_sided():
    assert rc.fisher_p(10, 0, 0, 10) == pytest.approx(2 / 184756)
    assert rc.fisher_p(5, 5, 5, 5) == pytest.approx(1.0)


@pytest.mark.parametrize("ref,alt,expected,status", [
    (12, 11, 1, "confirmed"), (3, 2, 1, "weak"), (25, 0, 1, "disagrees"), (19, 3, 1, "unclear"),
])
def test_classify(ref, alt, expected, status):
    c = rc.Counts(ref=ref, alt=alt, ref_fwd=ref // 2, alt_fwd=alt // 2)
    assert rc.classify(c, rc.Target("k", "1", 1, "A", "G", expected))["status"] == status


def test_look_alike_flag():
    c = rc.Counts(ref=10, alt=10, ref_fwd=5, alt_fwd=5, low_mapq=12)
    assert "look_alike" in rc.classify(c, rc.Target("k", "1", 1, "A", "G", 1))["flags"]


def test_recount_lowers_never_raises():
    base = {"ref": "A", "alt": "G", "ref_reads": 25, "alt_reads": 0, "expected": 1, "observed": 0, "flags": []}
    assert ev.recount_call("High", [], {**base, "status": "disagrees"})[0] == "Low"
    assert ev.recount_call("High", [], {**base, "status": "unclear"})[0] == "Medium"
    assert ev.recount_call("Low", [], {**base, "status": "confirmed"})[0] == "Low"
    level, why = ev.recount_call("High", [], {**base, "status": "confirmed", "flags": ["look_alike"],
                                              "low_mapq_reads": 9})
    assert level == "Medium" and any("look-alike" in w for w in why)
    gene = {"kind": "gene", "status": "disagrees", "positions": 40, "problems": ["rs1 (disagrees)"]}
    assert ev.recount_call("High", [], gene)[0] == "Low"


def test_pgx_target_uses_grch37_reference(tmp_path):
    """PharmCAT lists GRCh38 alleles; where GRCh37 carries the other letter (Factor V Leiden), the expected copies
    must be counted against GRCh37's reference, or a typical genotype reads as a disagreement."""
    ref = _reference()
    p0 = 500
    _, fa = _write(tmp_path, ref, [])
    other = "C" if ref[p0] != "C" else "G"
    pos = pd.DataFrame([{"chrom37": "chr1", "pos37": p0 + 1, "ref": other, "alts": [ref[p0]], "gt": "0/0",
                         "gene": "F5", "rsid": "rs6025", "status": "variant"}])
    (t,) = rc.pgx_targets(pos, fa, lambda ch: f"chr{ch}", None)
    assert (t.ref, t.alt, t.expected) == (ref[p0], other, 2)


def test_gene_verdict_worst_position():
    sites = pd.DataFrame([
        {"source": "pgx", "gene": "CYP2C19", "rsid": "rs1", "chrom": "10", "pos": 1, "status": "confirmed", "flags": [],
         "note": None, "alt": "A", "alt_reads": 0, "reads": 30, "observed": 0, "lod": 9.0},
        {"source": "pgx", "gene": "CYP2C19", "rsid": "rs2", "chrom": "10", "pos": 2, "status": "disagrees",
         "flags": [], "note": None, "alt": "A", "alt_reads": 14, "reads": 30, "observed": 1, "lod": 9.0},
        {"source": "pgx", "gene": "CYP2D6", "rsid": "rs3", "chrom": "22", "pos": 3, "status": "disagrees", "flags": [],
         "note": None, "alt": "A", "alt_reads": 14, "reads": 30, "observed": 1, "lod": 9.0},
    ])
    genes = pd.DataFrame([{"gene": "CYP2C19", "call_source": "PharmCAT (from your VCF)"},
                          {"gene": "CYP2D6", "call_source": "Cyrius"}])
    v = rc.gene_verdicts(sites, genes)
    assert v["CYP2C19"]["status"] == "disagrees" and "rs2 (disagrees)" in v["CYP2C19"]["problems"]
    assert "CYP2D6" not in v                     # Cyrius reads the raw reads itself


def test_long_deletion_is_not_called_a_repeat(tmp_path):
    """A 150-letter deletion can't be shown whole by 150-letter reads; the note said it sat in a long repeat."""
    ref = _reference()
    af, fa = _write(tmp_path, ref, [(s, ref[s:s + 150], "150M") for s in range(400, 500, 10)])
    t = rc.Target("k", "1", 501, ref[500:651], ref[500], expected=1)
    c = rc.count_site(af, fa, "chr1", t)
    assert "150 letters long" in c.note and "repeat" not in c.note
