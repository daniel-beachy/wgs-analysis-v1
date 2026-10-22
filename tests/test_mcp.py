"""`wgs mcp` (Step 12): a read-only MCP server over the newest release."""

from __future__ import annotations

import io
import json

import pandas as pd

from wgs.mcp_server import handle, serve


def _releases(tmp_path):
    rd = tmp_path / "releases"
    (rd / "objects").mkdir(parents=True)
    pd.DataFrame([{"claim_id": "c1", "gene": "TPMT", "overall": "Strong", "statement": "About c1."},
                  {"claim_id": "c2", "gene": "GENEB", "overall": "Limited", "statement": "About c2."}]
                 ).to_parquet(rd / "objects" / "claims.parquet")
    pd.DataFrame([{"rank": 1, "kind_label": "Medicine gene", "section": "pgx", "title": "TPMT", "statement": "Slower.",
                   "notable_because": "Intermediate", "overall": "Strong", "evidence_level": "Strong",
                   "call_level": "High", "sensitive": False, "ref_table": "claims", "ref_key": "c1",
                   "highlight_id": "h1"}]).to_parquet(rd / "objects" / "highlights.parquet")
    (rd / "objects" / "checks.json").write_text(json.dumps({"counts": {"pass": 2, "fail": 0}, "checks": []}))
    manifest = {"id": "r2", "created": "2026-10-22", "pipeline_version": "0.1.0", "sample": "S",
                "knowledge": {"versions": {"clinvar": "2026-10-19"}},
                "tables": {"claims": {"path": "objects/claims.parquet", "rows": 2, "description": "Claims"},
                           "highlights": {"path": "objects/highlights.parquet", "rows": 1, "description": "H"}},
                "documents": {"checks": {"path": "objects/checks.json"}}}
    (rd / "r2.json").write_text(json.dumps(manifest))
    (rd / "index.json").write_text(json.dumps({"releases": [{"id": "r2", "manifest": "r2.json"}], "latest": "r2"}))
    return rd


def _call(rd, name, args=None, mid=1):
    reply = handle(rd, {"jsonrpc": "2.0", "id": mid, "method": "tools/call",
                        "params": {"name": name, "arguments": args or {}}})
    text = reply["result"]["content"][0]["text"]
    return reply["result"]["isError"], text


def test_handshake_and_tool_list(tmp_path):
    rd = _releases(tmp_path)
    init = handle(rd, {"jsonrpc": "2.0", "id": 0, "method": "initialize", "params": {"protocolVersion": "2025-03-26"}})
    assert init["result"]["protocolVersion"] == "2025-03-26" and "guide" in init["result"]["instructions"]
    unknown = handle(rd, {"jsonrpc": "2.0", "id": 0, "method": "initialize", "params": {"protocolVersion": "1999"}})
    assert unknown["result"]["protocolVersion"] == "2025-06-18"
    assert handle(rd, {"jsonrpc": "2.0", "method": "notifications/initialized"}) is None
    names = [t["name"] for t in handle(rd, {"jsonrpc": "2.0", "id": 1, "method": "tools/list"})["result"]["tools"]]
    assert names == ["guide", "overview", "highlights", "result", "query"]
    assert handle(rd, {"jsonrpc": "2.0", "id": 2, "method": "resources/list"})["error"]["code"] == -32601


def test_tools_read_the_release(tmp_path):
    rd = _releases(tmp_path)
    err, text = _call(rd, "overview")
    ov = json.loads(text)
    assert not err and ov["release"] == "r2" and ov["self_checks"] == {"pass": 2, "fail": 0}
    assert "sample" not in ov and ov["documents"] == ["doc_checks"]
    err, text = _call(rd, "highlights")
    assert not err and json.loads(text)[0]["ref_key"] == "c1"
    err, text = _call(rd, "result", {"table": "claims", "key": "c2"})
    assert not err and json.loads(text)[0]["gene"] == "GENEB"
    err, text = _call(rd, "query", {"sql": "SELECT gene FROM claims ORDER BY gene", "limit": 1})
    assert not err and json.loads(text) == [{"gene": "GENEB"}]
    err, text = _call(rd, "guide")
    assert not err and "Guide for AI assistants" in text


def test_nothing_can_be_written_or_read_outside_the_release(tmp_path):
    rd = _releases(tmp_path)
    (tmp_path / "secret.csv").write_text("a\n1\n")
    for sql in ["DELETE FROM claims", "COPY claims TO 'x.csv'", "SELECT 1; SELECT 2",
                f"SELECT * FROM read_csv('{(tmp_path / 'secret.csv').as_posix()}')"]:
        err, text = _call(rd, "query", {"sql": sql})
        assert err, sql
    err, text = _call(rd, "result", {"table": "claims", "key": "nope"})
    assert err and "no claims row" in text
    err, text = _call(rd, "result", {"table": "variants", "key": "1"})
    assert err


def test_serve_speaks_one_json_message_per_line(tmp_path):
    rd = _releases(tmp_path)
    lines = [{"jsonrpc": "2.0", "id": 1, "method": "ping"}, {"jsonrpc": "2.0", "method": "notifications/initialized"}]
    stdin = io.StringIO("\n".join(json.dumps(x) for x in lines) + "\nnot json\n")
    stdout = io.StringIO()
    serve(rd, stdin, stdout)
    replies = [json.loads(x) for x in stdout.getvalue().splitlines()]
    assert replies == [{"jsonrpc": "2.0", "id": 1, "result": {}},
                       {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "parse error"}}]
