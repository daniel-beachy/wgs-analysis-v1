"""`wgs export-site` (ADR-028): the public web copy of the demo release."""

from __future__ import annotations

import json

import duckdb
import pandas as pd
import pytest

from wgs import site


def _releases(tmp_path, truth=True, statement="About c1."):
    rd = tmp_path / "releases"
    (rd / "objects").mkdir(parents=True)
    rows = [{"chrom": c, "chrom_order": i, "pos": p, "rsid": f"rs{i}{p}", "af_1kg": 0.25, "gnomad_popmax_af": 0.125,
             "gnomad_nhomalt": 3.0, "revel": 0.5, "am_score": 0.5}
            for i, c in enumerate(["1", "2", "3", "X"], start=1) for p in range(1000, 1300)]
    pd.DataFrame(rows).sample(frac=1, random_state=1).to_parquet(rd / "objects" / "annotations-a.parquet")
    pd.DataFrame([{"claim_id": "c1", "statement": statement}]).to_parquet(rd / "objects" / "claims-a.parquet")
    (rd / "objects" / "qc-a.json").write_text(json.dumps({"sections": {"reference": {"fasta": "/Volumes/X/ref.fa"}}}))
    tables = {n: {"path": f"objects/{n}-a.parquet", "bytes": (rd / "objects" / f"{n}-a.parquet").stat().st_size,
                  "sha256": "x", "rows": r} for n, r in [("annotations", len(rows)), ("claims", 1)]}
    manifest = {"schema": 1, "id": "r1", "sample": "HG002", "previous": None, "tables": tables,
                "documents": {"qc": {"path": "objects/qc-a.json", "bytes": 1, "sha256": "x"}}}
    if truth:
        manifest["truth_set"] = {"label": "Genome in a Bottle HG002"}
    (rd / "r1").mkdir()
    (rd / "r1" / "manifest.json").write_text(json.dumps(manifest))
    old = {"id": "r0", "sample": "HG002", "manifest": "r0/manifest.json", "previous": None}
    entry = {"id": "r1", "sample": "HG002", "manifest": "r1/manifest.json", "previous": "r0"}
    (rd / "index.json").write_text(json.dumps({"schema": 1, "releases": [entry, old], "latest": "r1"}))
    return rd


def _manifest(out):
    idx = json.loads((out / "data" / "index.json").read_text())
    return idx, json.loads((out / "data" / idx["releases"][0]["manifest"]).read_text())


def test_a_personal_release_is_refused(tmp_path):
    with pytest.raises(ValueError, match="public truth set"):
        site.export(_releases(tmp_path, truth=False), tmp_path / "site")
    with pytest.raises(ValueError, match="private ID"):
        site.export(_releases(tmp_path / "b"), tmp_path / "site", forbidden=["HG002"])
    assert not (tmp_path / "site" / "data").exists()


def test_export_leaves_out_columns_and_local_paths(tmp_path):
    out = tmp_path / "site"
    res = site.export(_releases(tmp_path), out)
    idx, m = _manifest(out)
    assert [r["id"] for r in idx["releases"]] == ["r1"] and res["split"] == {}
    con = duckdb.connect()
    ann = out / "data" / m["tables"]["annotations"]["path"]
    cols = {c[0]: c[1] for c in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{ann.as_posix()}')").fetchall()}
    assert not {"revel", "am_score", "gnomad_nhomalt"} & set(cols)
    assert cols["af_1kg"] == cols["gnomad_popmax_af"] == "FLOAT"
    assert set(m["export"]["dropped_columns"]["annotations"]) == set(site.DROP["annotations"])
    qc = json.loads((out / "data" / m["documents"]["qc"]["path"]).read_text())
    assert qc["sections"]["reference"]["fasta"] == "ref.fa"
    assert (out / "data" / m["tables"]["claims"]["path"]).exists()


def test_large_tables_are_split_by_chromosome(tmp_path, monkeypatch):
    monkeypatch.setattr(site, "PART_BYTES", 1)  # every chromosome gets its own part
    out = tmp_path / "site"
    res = site.export(_releases(tmp_path), out)
    _, m = _manifest(out)
    parts = m["tables"]["annotations"]["parts"]
    assert res["split"] == {"annotations": 4} and len(parts) == 4
    assert m["tables"]["annotations"]["path"] == parts[0]["path"]
    con = duckdb.connect()
    files = ", ".join(f"'{(out / 'data' / p['path']).as_posix()}'" for p in parts)
    rows = con.execute(f"SELECT chrom_order, pos FROM read_parquet([{files}])").fetchall()
    assert len(rows) == m["tables"]["annotations"]["rows"] == 1200 and rows == sorted(rows)
    spans = [con.execute(f"SELECT min(chrom_order), max(chrom_order) FROM read_parquet({f})").fetchone()
             for f in files.split(", ")]
    assert all(a[1] < b[0] for a, b in zip(spans, spans[1:], strict=False))


def test_private_strings_stop_the_export(tmp_path):
    out = tmp_path / "site"
    with pytest.raises(ValueError, match="private strings"):
        site.export(_releases(tmp_path, statement="Kit KIT123 says hi."), out, forbidden=["KIT123"])
    assert not (out / "data").exists()
