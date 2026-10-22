"""What changed (Step 10): a simulated knowledge update between two releases."""

from __future__ import annotations

import json

import pandas as pd

from wgs.diff import diff, summarise


def _src(*ids):
    return json.dumps([{"source": i, "version": "x", "record": "r", "url": ""} for i in ids])


def _release(root, name, tables: dict[str, list[dict]]):
    d = root / name
    d.mkdir()
    out = {}
    for t, rows in tables.items():
        pd.DataFrame(rows).to_parquet(d / f"{t}.parquet")
        out[t] = {"path": f"{name}/{t}.parquet"}
    return out


def _claim(cid, overall="Moderate", stars=1, sources=("clinvar",), **kw):
    return {"claim_id": cid, "section": "health", "subject": cid, "gene": "G" + cid, "rsid": None,
            "category_label": "Pathogenic", "category": "pathogenic", "clinvar_stars": stars,
            "evidence_level": overall, "call_level": "High", "overall": overall, "statement": f"About {cid}.",
            "sensitive": False, "sources": _src(*sources), **kw}


def _simulate(tmp_path, evidence_model=("4", "4")):
    old = _release(tmp_path, "old", {
        "claims": [_claim("a"), _claim("b", sources=("pharmcat",)), _claim("gone")],
        "pgx_genes": [{"gene": "TPMT", "diplotype": "*1/*3A", "phenotype": "Intermediate", "activity_score": None},
                      {"gene": "CYP2C9", "diplotype": "*1/*2", "phenotype": "Intermediate", "activity_score": None}],
        "pgs_traits": [{"trait_id": "t1", "label": "Height", "percentile": 50.0, "agreement": "consistent",
                        "overall": "Limited"},
                       {"trait_id": "t2", "label": "BMI", "percentile": 40.0, "agreement": "consistent",
                        "overall": "Limited"}],
    })
    new = _release(tmp_path, "new", {
        "claims": [_claim("a", overall="Strong", stars=2), _claim("b", sources=("pharmcat",)),
                   _claim("new", sources=("clinvar", "clingen"))],
        "pgx_genes": [{"gene": "TPMT", "diplotype": "*1/*3A", "phenotype": "Likely intermediate",
                       "activity_score": None},
                      {"gene": "CYP2C9", "diplotype": "*1/*2", "phenotype": "Intermediate", "activity_score": None}],
        "pgs_traits": [{"trait_id": "t1", "label": "Height", "percentile": 50.4, "agreement": "consistent",
                        "overall": "Limited"},          # within 1 point: float noise, not news
                       {"trait_id": "t2", "label": "BMI", "percentile": 43.0, "agreement": "consistent",
                        "overall": "Limited"}],
    })
    prev = {"id": "old", "tables": old,
            "knowledge": {"versions": {"clinvar": "2026-10-04", "pharmcat": "3.4.0", "clingen": "2026-10-06"},
                          "evidence_model": evidence_model[0]}}
    return diff(prev, tmp_path, {k: tmp_path / v["path"] for k, v in new.items()},
                {"clinvar": "2026-10-18", "pharmcat": "3.5.0", "clingen": "2026-10-06"}, evidence_model[1])


def test_each_change_names_the_source_update_behind_it(tmp_path):
    ch = _simulate(tmp_path)
    assert ch["knowledge"] == {"clinvar": ["2026-10-04", "2026-10-18"], "pharmcat": ["3.4.0", "3.5.0"]}
    changed = {r["claim_id"]: r for r in ch["claims"]["changed"]}
    assert set(changed) == {"a"} and changed["a"]["kind"] == "regraded" and changed["a"]["why"] == ["clinvar"]
    assert changed["a"]["fields"]["overall"] == ["Moderate", "Strong"] and "_sources" not in changed["a"]
    assert [(r["claim_id"], r["why"]) for r in ch["claims"]["added"]] == [("new", ["clinvar"])]  # clingen unchanged
    assert [(r["claim_id"], r["why"]) for r in ch["claims"]["removed"]] == [("gone", ["clinvar"])]


def test_results_that_are_not_claims_are_compared(tmp_path):
    ch = _simulate(tmp_path)
    items = {(r["table"], r["key"]): r for r in ch["results"]["items"]}
    assert set(items) == {("pgx_genes", "TPMT"), ("pgs_traits", "t2")}
    assert items[("pgx_genes", "TPMT")]["fields"] == {"phenotype": ["Intermediate", "Likely intermediate"]}
    assert items[("pgx_genes", "TPMT")]["why"] == ["pharmcat"]
    assert items[("pgs_traits", "t2")]["why"] == [] and items[("pgs_traits", "t2")]["name"] == "BMI"
    assert "2 other results changed" in summarise(ch)


def test_grading_rule_change_is_named(tmp_path):
    ch = _simulate(tmp_path, evidence_model=("4", "5"))
    assert ch["evidence_model"] == ["4", "5"]
    assert next(r for r in ch["claims"]["changed"] if r["claim_id"] == "a")["why"] == ["clinvar", "evidence_model"]


def test_first_release_has_nothing_to_compare(tmp_path):
    ch = diff(None, tmp_path, {}, {}, "5")
    assert ch["first"] and summarise(ch) == "First release."


def test_a_new_table_is_named_and_unchanged_findings_are_not_counted(tmp_path):
    old = _release(tmp_path, "old", {"claims": [_claim("a")]})
    new = _release(tmp_path, "new", {"claims": [_claim("a")], "highlights": [{"rank": 1}]})
    prev = {"id": "old", "tables": old, "knowledge": {"versions": {"clinvar": "1"}, "evidence_model": "5"}}
    ch = diff(prev, tmp_path, {k: tmp_path / v["path"] for k, v in new.items()}, {"clinvar": "1"}, "5")
    assert ch["tables"] == {"added": ["highlights"], "removed": []}
    assert summarise(ch) == "new: highlights"
