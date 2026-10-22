"""Highlights (Step 11): fixed rules pick rows other stages already graded, and never change their grade."""

from __future__ import annotations

import json

import duckdb
import pandas as pd

from wgs import highlights
from wgs.checks import run as run_checks


def _src(i):
    return json.dumps([{"source": i, "version": "x", "record": "r", "url": f"https://example.org/{i}/r"}])


def _claim(cid, section, group, overall, **kw):
    return {"claim_id": cid, "section": section, "group": group, "subject": cid, "gene": "G" + cid,
            "category_label": "Pathogenic", "statement": f"About {cid}.", "evidence_level": overall,
            "evidence_reasons": [f"reason for {cid}"], "call_level": "High", "overall": overall,
            "sensitive": section != "pgx", "sources": _src("clinvar"), **kw}


def _drug(drug, source, tier="change", matched=True, gene="TPMT"):
    return {"drug": drug, "source": source, "tier": tier, "matched": matched, "genes": [gene]}


def _tables(tmp_path):
    rows = {
        "claims": [_claim("mono", "health", "monogenic", "Moderate"),
                   _claim("carr", "carrier", "carrier", "Limited"),
                   _claim("risk_weak", "health", "risk", "Limited"),          # Limited risk factor: not a highlight
                   _claim("risk_ok", "health", "risk", "Moderate"),
                   _claim("vus", "health", "uncertain", "Limited"),           # a VUS is not an answer
                   _claim("TPMT", "pgx", "pgx_gene", "Strong", gene="TPMT", statement="Slower.")],
        "pgx_genes": [{"gene": "TPMT", "diplotype": "*1/*3A", "phenotype": "Intermediate Metabolizer",
                       "phenotype_class": "intermediate", "overall": "Strong", "evidence_level": "Strong",
                       "call_level": "High"},
                      {"gene": "DPYD", "diplotype": "Reference/Reference", "phenotype": "Normal Metabolizer",
                       "phenotype_class": "normal", "overall": "Strong", "evidence_level": "Strong",
                       "call_level": "High"}],
        "pgx_drugs": [_drug("azathioprine", "CPIC"), _drug("mercaptopurine", "CPIC"),
                      _drug("mercaptopurine", "DPWG"),                   # CPIC and DPWG: still one medicine
                      _drug("thioguanine", "FDA label"),                 # FDA rows never set the bucket
                      _drug("cisplatin", "CPIC", matched=False),         # a guideline that didn't match you
                      _drug("fluorouracil", "CPIC", tier="standard", gene="DPYD")],
        "pgs_traits": [{"trait_id": "t_hi", "label": "High", "section": "health", "percentile": 97.2,
                        "agreement": "consistent", "agreement_text": "Studies agree.", "overall": "Moderate",
                        "evidence_level": "Moderate", "call_confidence": "High", "risk_pgs_id": None,
                        "best_pgs_id": "PGS000001", "abs_risk": None},
                       {"trait_id": "t_mid", "label": "Middle", "section": "health", "percentile": 60.0,
                        "agreement": "consistent", "agreement_text": "Studies agree.", "overall": "Moderate",
                        "evidence_level": "Moderate", "call_confidence": "High", "risk_pgs_id": None,
                        "best_pgs_id": "PGS000002", "abs_risk": None},
                       {"trait_id": "t_dis", "label": "Disputed", "section": "health", "percentile": 2.0,
                        "agreement": "disagree", "agreement_text": "Studies disagree.", "overall": "Limited",
                        "evidence_level": "Limited", "call_confidence": "High", "risk_pgs_id": None,
                        "best_pgs_id": "PGS000003", "abs_risk": None}],
        "trait_snps": [{"rsid": "rs1", "label": "Rare look", "effect_text": "Per copy…", "genotype_share": 0.1,
                        "overall": "Strong", "evidence_level": "Strong", "call_confidence": "High", "pmid": "123"},
                       {"rsid": "rs2", "label": "Common look", "effect_text": "Per copy…", "genotype_share": 0.5,
                        "overall": "Strong", "evidence_level": "Strong", "call_confidence": "High", "pmid": "456"}],
    }
    out = {}
    for name, r in rows.items():
        out[name] = tmp_path / f"{name}.parquet"
        pd.DataFrame(r).to_parquet(out[name])
    return out


def test_rules_pick_notable_rows_and_keep_their_grades(tmp_path):
    t = _tables(tmp_path)
    df = duckdb.sql(highlights.sql(t)).df()
    assert set(df.ref_key) == {"mono", "carr", "risk_ok", "TPMT", "t_hi", "rs1"}
    assert df["rank"].tolist() == list(range(1, 7))
    assert df.overall.tolist() == ["Strong", "Strong", "Moderate", "Moderate", "Moderate", "Limited"]
    pgx = df[df.kind == "pgx"].iloc[0]
    assert pgx.statement == "Slower." and pgx.notable_because.endswith("differs for 2 medicines")
    assert df[df.kind == "pgs"].iloc[0].notable_because == "In the top 5%: the 97th percentile (median across studies)"
    assert pgx.title == "TPMT"
    assert df.set_index("ref_key").sensitive.to_dict()["mono"]


def test_self_checks_pass_on_highlights_and_catch_an_upgraded_grade(tmp_path):
    t = _tables(tmp_path)
    t["highlights"] = tmp_path / "highlights.parquet"
    duckdb.sql(f"COPY ({highlights.sql({k: v for k, v in t.items() if k != 'highlights'})}) TO '{t['highlights']}'")
    status = {c["id"]: c["status"] for c in run_checks(t)["checks"]}
    assert status["highlights.match_source"] == "pass" and status["highlights.complete"] == "pass"
    bad = tmp_path / "bad.parquet"
    duckdb.sql(f"COPY (SELECT * REPLACE (CASE WHEN ref_key = 'carr' THEN 'Strong' ELSE overall END AS overall) "
               f"FROM '{t['highlights']}') TO '{bad}'")
    status = {c["id"]: c for c in run_checks({**t, "highlights": bad})["checks"]}
    assert status["highlights.match_source"]["status"] == "fail"
    assert status["highlights.match_source"]["examples"] == ["claims:carr"]


def test_nothing_to_highlight_without_inputs():
    assert highlights.sql({}) is None


def test_percentiles_are_worded_like_the_dashboard():
    words = duckdb.sql(f"SELECT {highlights._percentile('x')} FROM (VALUES (0.2), (2.0), (3.4), (97.4), (99.9)) t(x)")
    assert [w for (w,) in words.fetchall()] == ["below the 1st percentile", "the 2nd percentile", "the 3rd percentile",
                                                "the 97th percentile", "above the 99th percentile"]
