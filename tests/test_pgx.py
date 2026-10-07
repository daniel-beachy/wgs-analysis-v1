"""Pharmacogenomics: per-position genotyping for PharmCAT, chain lifting, T1K parsing and guideline classification."""

import gzip

import pytest

from wgs.knowledge.pharmcat import Chain
from wgs.stages import pgx


def row(pos, ref, alt, gt, gq=50, flt="PASS"):
    return {"pos": pos, "ref": ref, "alt": alt, "gt": gt, "gq": gq, "filter": flt}


P = {"pos37": 100, "ref": "C", "alts": ["T"]}


def test_norm_trims_shared_bases():
    assert pgx.norm(10, "CAT", "C") == (11, "AT", "")
    assert pgx.norm(10, "CAT", "CGT") == (11, "A", "G")
    assert pgx.norm(10, "A", "G") == (10, "A", "G")


@pytest.mark.parametrize(("rows", "callable_", "gt", "status"), [
    ([], True, "0/0", "reference"),                                   # callable and no variant → confident reference
    ([], False, "./.", "not_callable"),                               # absence of evidence is not "normal"
    ([row(100, "C", "C", "0/0", flt="RefCall")], False, "0/0", "reference"),
    ([row(100, "C", "T", "0/1")], True, "0/1", "variant"),
    ([row(100, "C", "T", "1/1")], True, "1/1", "variant"),
    ([row(100, "C", "T", "0/1", gq=5)], True, "./.", "low_quality"),
    ([row(100, "C", "T", "0/1", flt="LowQual")], True, "./.", "filtered"),
    ([row(100, "C", "G", "0/1")], True, "./.", "other_variant"),
    ([row(100, "C", "T", "./.")], True, "./.", "no_call"),
])
def test_genotype_position(rows, callable_, gt, status):
    out = pgx.genotype_position(P, "C", rows, callable_, haploid=False)
    assert (out["gt"], out["status"]) == (gt, status)


def test_genotype_position_reference_swap():
    # Your reference carries PharmCAT's ALT: no variant called means you are homozygous for allele 1.
    assert pgx.genotype_position(P, "T", [], True, haploid=False)["gt"] == "1/1"
    # A variant back to PharmCAT's REF is a het 0/1.
    assert pgx.genotype_position(P, "T", [row(100, "T", "C", "0/1")], True, haploid=False)["gt"] == "0/1"
    assert pgx.genotype_position(P, "G", [], True, haploid=False)["status"] == "ref_mismatch"


def test_genotype_position_multiallelic_and_haploid():
    p = {"pos37": 100, "ref": "C", "alts": ["T", "A"]}
    assert pgx.genotype_position(p, "C", [row(100, "C", "A", "0/1")], True, haploid=False)["gt"] == "0/2"
    assert pgx.genotype_position(p, "C", [row(100, "C", "T", "0/1"), row(100, "C", "A", "0/1")], True,
                                 haploid=False)["gt"] == "1/2"
    assert pgx.genotype_position(P, "C", [row(100, "C", "T", "1/1")], True, haploid=True)["gt"] == "1"
    assert pgx.genotype_position(P, "C", [row(100, "C", "T", "0/1")], True, haploid=True)["status"] == "het_haploid"


def test_chain_lift(tmp_path):
    # chr1: 0-100 → chrA 1000-1100 (+); chr2: 0-50 → chrB (-) in a 500-long query.
    text = ("chain 100 chr1 1000 + 0 100 chrA 5000 + 1000 1100 1\n60 10 10\n30\n\n"
            "chain 90 chr2 1000 + 0 50 chrB 500 - 0 50 2\n50\n\n")
    path = tmp_path / "t.chain.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    c = Chain(path)
    assert c.lift("chr1", 1) == ("chrA", 1001, "+")
    assert c.lift("chr1", 60) == ("chrA", 1060, "+")
    assert c.lift("chr1", 65) is None                 # in the gap
    assert c.lift("chr1", 71) == ("chrA", 1071, "+")
    assert c.lift("chr1", 59, length=3) is None       # span crosses into the gap
    assert c.lift("chr2", 1) == ("chrB", 500, "-")
    assert c.lift("chr9", 1) is None


def test_parse_t1k(tmp_path):
    path = tmp_path / "s_genotype.tsv"
    path.write_text("HLA-A\t2\tA*03:01:01:01\t100\t60\tA*02:01:01:01\t90\t55\n"
                    "HLA-B\t1\tB*44:02:01\t100\t60\t.\t0\t0\n"
                    "KIR2DL1\t2\tx\t1\t1\ty\t1\t1\n")
    got = pgx.parse_t1k(path)
    assert got["HLA-A"]["alleles"] == ["*03:01", "*02:01"]
    assert got["HLA-B"]["alleles"] == ["*44:02", "*44:02"]   # one allele reported = homozygous
    assert "KIR2DL1" not in got


@pytest.mark.parametrize(("phenotype", "cls"), [
    ("Poor Metabolizer", "poor"),
    ("Intermediate Metabolizer", "intermediate"),
    ("Ultrarapid Metabolizer", "ultrarapid"),
    ("Rapid Metabolizer", "rapid"),
    ("Normal Function", "normal"),
    ("Decreased Function", "decreased"),
    ("*57:01 negative", "negative"),
    ("*58:01 positive", "risk"),
    ("No Result", "no_result"),
    (pgx.GENOTYPE_ONLY, "genotype"),
])
def test_pheno_class(phenotype, cls):
    assert pgx._pheno_class(phenotype) == cls


def _a(cls="Strong", dose=False, alt=False, other=False):
    return {"classification": cls, "dosingInformation": dose, "alternateDrugAvailable": alt,
            "otherPrescribingGuidance": other}


@pytest.mark.parametrize(("ann", "rec", "cf", "expected"), [
    (_a(dose=True), "Consider a 50% reduction.", False, "change"),
    (_a(alt=True), "Avoid clopidogrel if possible.", False, "change"),
    # DPWG venlafaxine: the wording opens "It is not possible to offer..." but the curated flags say avoid/adjust.
    (_a("Unspecified", dose=True, alt=True), "It is not possible to offer adequately substantiated advice. Avoid.",
     False, "change"),
    (_a(other=True), "Use standard doses with therapeutic dose monitoring.", False, "note"),
    (_a(), "Initiate therapy with recommended starting dose.", False, "standard"),
    (_a("No recommendation", dose=True), "No recommendation based on insufficient evidence.", False, "standard"),
    (_a(alt=True), "Ivacaftor is not recommended", True, "standard"),   # CFTR guidance only applies with CF
    (_a(dose=True), "", False, "standard"),
])
def test_tier_uses_curated_flags_not_wording(ann, rec, cf, expected):
    assert pgx.tier(ann, rec, cf) == expected
