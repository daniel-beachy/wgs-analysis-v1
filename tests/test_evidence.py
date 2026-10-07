import json
from pathlib import Path

import duckdb

from wgs import evidence as ev
from wgs.diff import diff, summarise


def test_stars_map_to_evidence_and_reasons():
    level, reasons = ev.clinvar_evidence(3, "pathogenic", popmax_af=0.0001)
    assert level == "Strong"
    assert any("expert panel" in r for r in reasons)
    assert ev.clinvar_evidence(1, "pathogenic")[0] == "Limited"
    assert ev.clinvar_evidence(2, "likely_pathogenic", popmax_af=0.0001)[0] == "Moderate"


def test_common_pathogenic_variant_is_downgraded():
    level, reasons = ev.clinvar_evidence(4, "pathogenic", popmax_af=0.08)
    assert level == "Limited" and any("BA1" in r for r in reasons)
    assert ev.clinvar_evidence(3, "pathogenic", popmax_af=0.02)[0] == "Moderate"


def test_disputed_gene_and_conflicts_cap_evidence():
    assert ev.clinvar_evidence(3, "pathogenic", gene_validity="Refuted")[0] == "Limited"
    assert ev.clinvar_evidence(3, "pathogenic", gene_validity="Limited")[0] == "Moderate"
    assert ev.clinvar_evidence(2, "conflicting")[0] == "Limited"


def test_drug_response_never_strong_from_clinvar_alone():
    level, reasons = ev.clinvar_evidence(3, "drug_response")
    assert level == "Moderate" and any("CPIC" in r for r in reasons)


def test_call_confidence_and_overall_is_weaker_grade():
    assert ev.call_confidence(filter="PASS", gq=40, dp=35, vaf=0.5, zygosity="het")[0] == "High"
    assert ev.call_confidence(filter="PASS", gq=25, dp=12, vaf=0.5, zygosity="het")[0] == "Medium"
    assert ev.call_confidence(filter="PASS", gq=40, dp=35, vaf=0.1, zygosity="het")[0] == "Medium"
    assert ev.call_confidence(filter="LowQual", gq=40, dp=35, vaf=0.5, zygosity="het")[0] == "Low"
    assert ev.call_confidence(filter="PASS", gq=40, dp=35, vaf=0.5, zygosity="het", provider="disagree")[0] == "Low"
    assert ev.overall("Strong", "Medium") == "Moderate"
    assert ev.overall("Limited", "High") == "Limited"
    assert ev.overall("Strong", "Not callable") == "Not callable"


def _claims(path: Path, rows: list[dict]) -> Path:
    con = duckdb.connect()
    con.execute("CREATE TABLE c(claim_id VARCHAR, section VARCHAR, subject VARCHAR, gene VARCHAR, rsid VARCHAR, "
                "category_label VARCHAR, category VARCHAR, clinvar_stars INT, evidence_level VARCHAR, "
                "call_level VARCHAR, overall VARCHAR, gene_validity VARCHAR, statement VARCHAR, kind VARCHAR)")
    for r in rows:
        con.execute("INSERT INTO c VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", [
            r["id"], "health", r.get("gene", "G"), r.get("gene", "G"), None, r["cat"], r["cat"], r["stars"],
            r["ev"], "High", r["ev"], None, r.get("s", "s"), "clinvar"])
    con.execute(f"COPY c TO '{path}' (FORMAT parquet)")
    return path


def test_diff_reports_added_removed_changed_and_versions(tmp_path):
    _claims(tmp_path / "old.parquet", [
        {"id": "a", "cat": "Pathogenic", "stars": 1, "ev": "Limited"},
        {"id": "b", "cat": "Risk factor", "stars": 0, "ev": "Limited"},
        {"id": "d", "cat": "Drug response", "stars": 1, "ev": "Moderate", "s": "old text"},
    ])
    new = _claims(tmp_path / "new.parquet", [
        {"id": "a", "cat": "Pathogenic", "stars": 3, "ev": "Strong"},
        {"id": "d", "cat": "Drug response", "stars": 1, "ev": "Moderate", "s": "new text"},
        {"id": "c", "cat": "Drug response", "stars": 2, "ev": "Moderate"},
    ])
    prev = {"id": "r1", "tables": {"claims": {"path": "old.parquet"}},
            "knowledge": {"versions": {"clinvar": "2025-10-06"}, "evidence_model": "1"}}
    ch = diff(prev, tmp_path, {"claims": new}, {"clinvar": "2026-10-04"}, "1")
    assert ch["knowledge"] == {"clinvar": ["2025-10-06", "2026-10-04"]}
    assert [r["claim_id"] for r in ch["claims"]["added"]] == ["c"]
    assert [r["claim_id"] for r in ch["claims"]["removed"]] == ["b"]
    kinds = {r["claim_id"]: r["kind"] for r in ch["claims"]["changed"]}
    assert kinds == {"a": "regraded", "d": "reworded"}
    assert ch["claims"]["counts"]["changed"] == 2 and ch["claims"]["counts"]["regraded"] == 1
    changed = ch["claims"]["changed"][0]
    assert changed["claim_id"] == "a"
    assert changed["fields"]["clinvar_stars"] == [1, 3] and changed["fields"]["evidence_level"] == ["Limited", "Strong"]
    assert "clinvar 2025-10-06 → 2026-10-04" in summarise(ch)
    assert "1 regraded, 1 reworded" in summarise(ch)
    json.dumps(ch)


def test_first_release_has_no_diff(tmp_path):
    ch = diff(None, tmp_path, {}, {"clinvar": "x"}, "1")
    assert ch["first"] and summarise(ch) == "First release."


def test_majority_allele_note():
    from wgs.stages.annotate import _alt_is_major, _major_note
    assert _alt_is_major({"af_1kg": 0.994, "gnomad_popmax_af": 1.0}) == (0.994, "1000 Genomes, all populations")
    assert _alt_is_major({"af_1kg": 0.26, "gnomad_popmax_af": 0.6}) is None
    assert _alt_is_major({"af_1kg": None, "gnomad_popmax_af": 0.7})[1].startswith("gnomAD")
    assert _alt_is_major({"af_1kg": None, "gnomad_popmax_af": None}) is None
    note = _major_note((0.99, "x"), 2, ["no_sequence_alteration"])
    # Only sourced facts: never infer which allele "the effect" belongs to (it contradicted risk-factor claims)
    assert "do not carry" not in note and "most common one" in note and "no sequence alteration" in note
    assert "most common" not in _major_note((0.6, "x"), 1)


def test_ba1_withholds_rare_disease_roles():
    from wgs import evidence as ev
    assert ev.too_common(popmax_af=0.07, kg_af=0.02, stars=1).startswith("Carried on up to 7%")
    assert ev.too_common(popmax_af=0.04, kg_af=0.02, stars=1) is None
    assert ev.too_common(popmax_af=0.08, kg_af=None, stars=3) is None   # expert panel already weighed frequency


def test_predicted_wording_respects_existing_clinvar_entry():
    from wgs import evidence as ev
    kw = dict(consequence="missense", am_score=0.9, revel=0.95, loeuf=None, popmax_af=None, kg_af=None,
              recessive_gene=False)
    assert "no lab" in ev.predicted_evidence(**kw)[1][-1]
    assert "uncertain significance (1★)" in ev.predicted_evidence(**kw, clinvar="uncertain significance (1★)")[1][-1]


def test_conflict_and_condition_text_is_stable():
    from wgs.stages.annotate import _conditions, _conflicts
    assert _conflicts("Pathogenic(2)|Likely_pathogenic(1)") == "Pathogenic 2, Likely pathogenic 1"
    assert _conflicts("Pathogenic (2)|Likely pathogenic (1)") == "Pathogenic 2, Likely pathogenic 1"
    assert _conditions("B disease|a disease|not provided") == _conditions("a disease|B disease")
