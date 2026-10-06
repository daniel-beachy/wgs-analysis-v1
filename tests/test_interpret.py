"""Evidence model v2: inheritance, carrier/affected roles, ACMG reporting rules, computational predictions."""

import duckdb
import pytest

from wgs import evidence as ev
from wgs import inheritance as inh
from wgs.interpret import Kb


def test_from_text_and_cues():
    assert inh.from_text("Autosomal recessive") == "AR"
    assert inh.from_text("Semidominant") == "SD"
    assert inh.from_text("X-linked") == "XL"
    assert inh.from_text("Unknown") is None
    assert inh.cues("Autosomal recessive nonsyndromic hearing loss 1A") == {"AR"}
    assert inh.cues("Duchenne, X-linked") == {"XL"}
    assert inh.cues("Breast cancer") == set()


@pytest.mark.parametrize(("modes", "zyg", "chrom", "expected"), [
    ({"AR"}, "het", "13", "carrier"),
    ({"AR"}, "hom", "13", "affected"),
    ({"AD"}, "het", "1", "possible"),
    ({"AD", "AR"}, "het", "5", "possible"),         # SDHA-like: dominant and recessive conditions
    ({"XLR"}, "hemi", "X", "affected"),
    ({"XLR"}, "het", "X", "carrier"),
    ({"DG"}, "het", "1", "unknown"),
    (set(), "het", "1", "unknown"),
])
def test_role(modes, zyg, chrom, expected):
    assert ev.role(modes, zyg, chrom)[0] == expected


def test_acmg_reportable_rules():
    kw = {"consequence": "missense_variant", "aa_change": "p.Cys282Tyr"}
    assert ev.acmg_reportable("Any P/LP", copies=1, zygosity="het", **kw)
    assert not ev.acmg_reportable(None, copies=1, zygosity="het", **kw)
    assert not ev.acmg_reportable("Biallelic P/LP only", copies=1, zygosity="het", **kw)
    assert ev.acmg_reportable("Biallelic P/LP only", copies=1, zygosity="het", compound=True, **kw)
    assert ev.acmg_reportable("Biallelic P/LP only", copies=2, zygosity="hom", **kw)
    hfe = "p.C282Y homozygous only"
    assert not ev.acmg_reportable(hfe, copies=1, zygosity="het", consequence="missense_variant", aa_change="C282Y")
    assert ev.acmg_reportable(hfe, copies=2, zygosity="hom", consequence="missense_variant", aa_change="C282Y")
    assert not ev.acmg_reportable(hfe, copies=2, zygosity="hom", consequence="missense_variant", aa_change="H63D")
    ttn = "P/LP truncating variants only"
    assert ev.acmg_reportable(ttn, copies=1, zygosity="het", consequence="stop_gained", aa_change=None)
    assert not ev.acmg_reportable(ttn, copies=1, zygosity="het", consequence="missense_variant", aa_change=None)


def _pred(**kw):
    base = {"consequence": "missense_variant", "am_score": None, "revel": None, "loeuf": None,
            "popmax_af": None, "kg_af": None, "recessive_gene": False}
    return ev.predicted_evidence(**(base | kw))


def test_predicted_missense_thresholds():
    assert _pred(revel=0.95)[0] == "Limited"
    assert _pred(revel=0.70, am_score=0.8)[0] == "Limited"
    assert _pred(revel=0.70, am_score=0.3)[0] is None           # predictors disagree
    assert _pred(revel=0.50, am_score=0.99)[0] is None           # REVEL present but weak
    assert _pred(am_score=0.95)[0] == "Limited"                  # no REVEL: AlphaMissense alone, high bar
    assert _pred(am_score=0.80)[0] is None
    assert _pred(revel=0.99, popmax_af=0.01)[0] is None          # common → not flagged


def test_predicted_lof_needs_constraint_or_recessive_gene():
    assert _pred(consequence="stop_gained", loeuf=0.3)[0] == "Limited"
    assert _pred(consequence="stop_gained", loeuf=1.5)[0] is None
    assert _pred(consequence="frameshift&splice_region", recessive_gene=True)[0] == "Limited"
    level, reasons = _pred(consequence="stop_gained", loeuf=0.3)
    assert any("capped at Limited" in r for r in reasons)


@pytest.fixture
def kb(tmp_path):
    con = duckdb.connect()

    def table(name, sql):
        p = tmp_path / f"{name}.parquet"
        con.execute(f"COPY ({sql}) TO '{p}' (FORMAT parquet)")
        return p

    return Kb({
        "hpo_inh": table("hpo_inh", "SELECT * FROM (VALUES ('OMIM:220290', 'DFNB1A', ['AR', 'DG']), "
                                    "('OMIM:601544', 'DFNA3A', ['AD'])) t(disease_id, disease_name, modes)"),
        "hpo_gd": table("hpo_gd", "SELECT * FROM (VALUES ('GJB2', 'OMIM:220290', 'MENDELIAN'), "
                                  "('GJB2', 'OMIM:601544', 'MENDELIAN')) t(gene, disease_id, association_type)"),
        "mondo": table("mondo", "SELECT 'MONDO:0009076' mondo_id, 'autosomal recessive nonsyndromic hearing loss 1A' "
                                "AS name, 'Hearing loss present at birth.' AS definition, ['OMIM:220290'] AS xrefs"),
        "orpha_dis": None, "orpha_gd": None, "validity": None, "actionability": None, "constraint": None,
        "acmg": table("acmg", "SELECT 'HFE' gene, 'Miscellaneous' category, 'Hemochromatosis' condition, "
                              "'AR' inheritance, 'p.C282Y homozygous only' report"),
    })


def test_kb_resolves_condition_inheritance(kb):
    conds = kb.conditions("Nonsyndromic hearing loss 1A|not provided", "MONDO:MONDO:0009076,OMIM:220290|")
    assert len(conds) == 1
    c = conds[0]
    assert c.modes_from == "HPO" and "AR" in c.modes
    assert c.definition.startswith("Hearing loss")
    modes, basis = kb.modes_for("GJB2", conds)
    assert basis == "condition" and ev.role(modes, "het", "13")[0] == "carrier"


def test_kb_gene_fallback(kb):
    modes, basis = kb.modes_for("GJB2", [])
    assert basis == "gene" and modes >= {"AR", "AD"}
    assert ev.role(modes, "het", "13")[0] == "possible"
    assert kb.disease_gene("GJB2") and not kb.disease_gene("FOO")
    names = [c.name for c in kb.gene_conditions("GJB2")]
    assert any("hearing loss 1A" in n for n in names)
    assert kb.acmg["HFE"]["report"].startswith("p.C282Y")


def test_asserted_conditions_discount_gene_wide_submissions():
    from wgs.interpret import Condition
    from wgs.stages.annotate import _asserted

    conds = [Condition("DFNB1A", modes=["AR"]), Condition("KID syndrome", modes=["AD"])]
    subs = [{"class": "pathogenic", "conditions": ["DFNB1A"]} for _ in range(6)]
    subs.append({"class": "pathogenic", "conditions": ["DFNB1A", "KID syndrome", "A", "B"]})   # gene-wide bundle
    subs.append({"class": "benign", "conditions": ["KID syndrome"]})
    keep, other = _asserted(conds, subs, "pathogenic")
    assert [c.name for c in keep] == ["DFNB1A"] and [c.name for c in other] == ["KID syndrome"]
    assert _asserted(conds, [], "pathogenic") == (conds, [])
