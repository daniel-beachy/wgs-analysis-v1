"""Release self-checks and the sourced-content rules (no medical text written into the viewer)."""

import json
import re
from pathlib import Path

import duckdb

from wgs import checks
from wgs.knowledge.genes import _short_date, parse_medlineplus

ROOT = Path(__file__).resolve().parents[1]


def _parquet(tmp_path: Path, name: str, rows: list[dict]) -> Path:
    src, dest = tmp_path / f"{name}.json", tmp_path / f"{name}.parquet"
    src.write_text("\n".join(json.dumps(r) for r in rows))
    duckdb.execute(f"COPY (SELECT * FROM read_json('{src}', format='newline_delimited')) TO '{dest}' (FORMAT parquet)")
    return dest


def _drug(**kw):
    base = {"drug": "x", "source": "CPIC", "genes": ["CYP2D6"], "matched": True, "tier": "standard",
            "recommendation": "Use as usual.", "flags": [], "classification": "Strong", "messages": []}
    return {**base, **kw}


def _status(result: dict) -> dict:
    return {c["id"]: c["status"] for c in result["checks"]}


def test_checks_pass_on_consistent_data(tmp_path):
    t = {"pgx_drugs": _parquet(tmp_path, "pgx_drugs", [
        _drug(), _drug(drug="y", tier="change", flags=["dosingInformation"], recommendation="Halve the dose."),
        _drug(drug="w", matched=False, tier="note", recommendation=None, messages=["Use the flowchart."])])}
    st = _status(checks.run(t))
    assert st["pgx.tier_from_flags"] == st["pgx.matched_has_text"] == st["pgx.notes_are_guidance"] == "pass"
    assert st["claims.sourced"] == "skipped"     # missing tables skip, never fail


def test_checks_catch_the_bugs_they_were_written_for(tmp_path):
    t = {"pgx_drugs": _parquet(tmp_path, "pgx_drugs", [
        _drug(drug="venlafaxine", source="DPWG", flags=["alternateDrugAvailable"]),       # flagged but "standard"
        _drug(drug="z", recommendation=""),                                                # matched, no text
        _drug(drug="bupivacaine", source="FDA label", matched=False, tier="note", messages=["data caveat"])]),
        "claims": _parquet(tmp_path, "claims", [
            {"claim_id": "a", "group": "drug_response", "sensitive": True, "sources": "[]", "evidence_level": "Limited",
             "call_level": "High", "overall": "Limited"},
            {"claim_id": "a", "group": "variant", "sensitive": False, "sources": '[{"source": "clinvar"}]',
             "evidence_level": "Strong", "call_level": "High", "overall": "Strong"}])}
    r = checks.run(t)
    st = _status(r)
    for cid in ("pgx.tier_from_flags", "pgx.matched_has_text", "pgx.notes_are_guidance", "claims.unique",
                "claims.sourced", "claims.drug_response_not_sensitive"):
        assert st[cid] == "fail", cid
    assert r["counts"]["fail"] == 6


def test_medlineplus_parse(tmp_path):
    p = tmp_path / "ghr.xml"
    p.write_text("""<?xml version="1.0"?><summaries>
      <gene-summary><name>thiopurine S-methyltransferase</name><gene-symbol>TPMT</gene-symbol>
        <ghr-page>https://medlineplus.gov/genetics/gene/tpmt/</ghr-page><reviewed>2015-04</reviewed>
        <text-list><text><text-role>function</text-role><html xmlns="http://www.w3.org/1999/xhtml">
          <p>The TPMT gene   makes an enzyme.</p><p>Second paragraph.</p></html></text></text-list>
        <db-key-list><db-key><db>NCBI Gene</db><key>7172</key></db-key></db-key-list>
        <related-health-condition-list><related-health-condition><name>Thiopurine toxicity</name>
        </related-health-condition></related-health-condition-list>
      </gene-summary></summaries>""")
    genes, conds = parse_medlineplus(p)
    assert genes == [{"gene": "TPMT", "name": "thiopurine S-methyltransferase",
                      "url": "https://medlineplus.gov/genetics/gene/tpmt/", "reviewed": "2015-04", "ncbi_id": "7172",
                      "function": "The TPMT gene makes an enzyme.\nSecond paragraph.",
                      "conditions": ["Thiopurine toxicity"]}]
    assert conds == []


def test_short_date_is_safe_as_a_folder_name():
    assert _short_date("Tue, 06 Oct 2026 05:15:58 GMT") == "2026-10-06"
    assert _short_date("2026-10-06T05:15:58Z") == "2026-10-06"
    assert ":" not in _short_date("garbage: 1:2")


def test_viewer_has_no_hand_written_gene_or_drug_guidance():
    """Medical statements come from sourced tables. The viewer may *mention* genes/drugs as examples in teaching
    text, but must not hold lookup tables keyed by gene or drug (the old ROLE / DRUG_NOTE dictionaries)."""
    keyed = re.compile(r"^\s*'?(CYP\d\w*|HLA-[A-Z]+|[A-Z][A-Z0-9]{2,}\d|warfarin|codeine|clopidogrel|tamoxifen)'?"
                       r"\s*:\s*['`(]", re.M)
    files = [p for ext in ("*.svelte", "*.ts") for p in (ROOT / "viewer/src").rglob(ext) if not p.name.startswith("._")]
    offenders = [f"{p.relative_to(ROOT)}: {m.group(0).strip()}" for p in files for m in keyed.finditer(p.read_text())]
    assert not offenders, offenders
