"""Claude's take (Step 11): every paragraph cites, every citation resolves, and a changed citation marks it stale."""

from __future__ import annotations

import pandas as pd

from wgs import take

TAKE = """---
author: Claude (test)
written: 2026-10-22
written_for: old
---
## Medicines

You process some medicines **more slowly** [[pgx_genes:TPMT]], which matters *only* when prescribed one.

- A carrier result [[claims:c1]]
- Another [[claims:c1]]

This paragraph forgot to cite anything.

And this one cites a result that isn't there [[claims:nope]] and a table that can't be cited [[variants:1]].
"""


def _tables(tmp_path, name, overall):
    d = tmp_path / name
    d.mkdir()
    pd.DataFrame([{"claim_id": "c1", "gene": "G1", "subject": "c1", "section": "carrier", "overall": overall,
                   "statement": "One copy."}]).to_parquet(d / "claims.parquet")
    pd.DataFrame([{"gene": "TPMT", "diplotype": "*1/*3A", "phenotype": "Intermediate", "overall": "Strong"}]
                 ).to_parquet(d / "pgx_genes.parquet")
    return {"claims": d / "claims.parquet", "pgx_genes": d / "pgx_genes.parquet"}


def test_take_is_parsed_and_its_citations_checked(tmp_path):
    f = tmp_path / "take.md"
    f.write_text(TAKE)
    old = _tables(tmp_path, "old", "Limited")
    doc = take.build(f, old, lambda rid: old if rid == "old" else None)
    assert [b["type"] for b in doc["blocks"]] == ["h", "p", "ul", "p", "p"]
    assert doc["blocks"][1]["segments"][1] == {"t": "b", "v": "more slowly"}
    assert doc["blocks"][1]["segments"][3] == {"t": "cite", "ref": "pgx_genes:TPMT"}
    cites = {c["ref"]: c for c in doc["citations"]}
    assert cites["claims:c1"]["found"] and cites["claims:c1"]["page"] == "carrier"
    assert cites["pgx_genes:TPMT"]["title"] == "TPMT *1/*3A" and cites["pgx_genes:TPMT"]["page"] == "pgx"
    assert len(doc["problems"]) == 3 and doc["problems"][0].startswith("Cites nothing")
    assert "Not in this release: claims:nope" in doc["problems"]
    assert "Can't cite table variants: variants:1" in doc["problems"]
    assert not doc["stale"]


def test_take_goes_stale_when_a_cited_grade_changes(tmp_path):
    f = tmp_path / "take.md"
    f.write_text(TAKE)
    old, new = _tables(tmp_path, "old", "Limited"), _tables(tmp_path, "new", "Moderate")
    doc = take.build(f, new, lambda rid: old if rid == "old" else None)
    assert doc["stale"] and {c["ref"] for c in doc["citations"] if c["changed"]} >= {"claims:c1"}
    assert not next(c for c in doc["citations"] if c["ref"] == "pgx_genes:TPMT")["changed"]


def test_no_take_written(tmp_path):
    assert take.build(tmp_path / "missing.md", {}, lambda rid: None) is None
