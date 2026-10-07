import json
from types import SimpleNamespace

import pyarrow as pa
import pyarrow.parquet as pq
from typer.testing import CliRunner

from wgs import cli


def test_query_skips_non_json_documents(tmp_path, monkeypatch):
    rel = tmp_path / "r1"
    (rel / "objects").mkdir(parents=True)
    pq.write_table(pa.table({"x": [1, 2]}), rel / "objects" / "t.parquet")
    (rel / "objects" / "report.html").write_text("<html></html>")
    (rel / "objects" / "checks.json").write_text(json.dumps({"checks": []}))
    (rel / "manifest.json").write_text(json.dumps({
        "sample": "S",
        "tables": {"t": {"path": "r1/objects/t.parquet", "rows": 2}},
        "documents": {"pgx_report": {"path": "r1/objects/report.html"},
                      "checks": {"path": "r1/objects/checks.json"}},
    }))
    (tmp_path / "index.json").write_text(json.dumps({"releases": [{"id": "r1", "manifest": "r1/manifest.json"}]}))
    monkeypatch.setattr(cli, "_cfg", lambda: SimpleNamespace(releases_dir=tmp_path))
    res = CliRunner().invoke(cli.app, ["query", "select sum(x) as s from t"])
    assert res.exit_code == 0, res.output
    assert "3" in res.output
