"""`wgs mcp`: a local, read-only Model Context Protocol server over the newest release (Step 12).

An AI assistant that speaks MCP (Claude Code, Claude Desktop, …) starts this as a subprocess and talks JSON-RPC 2.0
over stdin/stdout, one message per line. It can read the release the same way `wgs query` does (ADR-027): only a
single SELECT over the release's own tables and documents, nothing written, nothing outside the releases folder. The
genome itself (FASTQ, CRAM, VCF) is never exposed; the assistant sees only rows it asks for.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import IO

from . import __version__
from .readonly import open_release, select
from .take import CITABLE

PROTOCOLS = ("2025-06-18", "2025-03-26", "2024-11-05")
MAX_ROWS = 1000
GUIDE = Path(__file__).resolve().parents[2] / "AGENTS.md"
INSTRUCTIONS = ("Read-only access to one person's whole-genome analysis (a dated release of graded results). "
                "Call `guide` before discussing results and follow it: give each result's grade and reasons, "
                "separate the science from the genotype call, ask before sensitive topics, never present anything "
                "as a diagnosis, and quote only the rows you need. Start with `highlights`.")

TOOLS = [
    {"name": "guide", "description": "How to read and talk about these results: every table and column, the "
     "grading rules, and how to explain results responsibly. Read this first.",
     "inputSchema": {"type": "object", "properties": {}}},
    {"name": "overview", "description": "The release: its id and date, knowledge versions, self-check results, "
     "and every table (row count, description) and JSON document (as a doc_<name> view).",
     "inputSchema": {"type": "object", "properties": {}}},
    {"name": "highlights", "description": "The notable results from every section, strongest evidence first, "
     "each with why it was picked and the table and key of the row it came from.",
     "inputSchema": {"type": "object", "properties": {}}},
    {"name": "result", "description": "Every column of one result, by table and key (e.g. claims + claim_id, "
     "pgx_genes + gene, pgs_traits + trait_id, highlights + highlight_id).",
     "inputSchema": {"type": "object", "required": ["table", "key"], "properties": {
         "table": {"type": "string", "enum": sorted(CITABLE)}, "key": {"type": "string"}}}},
    {"name": "query", "description": "Run one read-only SQL SELECT (DuckDB dialect) over the release's tables and "
     "doc_<name> views. Prefer narrow queries over dumping whole tables.",
     "inputSchema": {"type": "object", "required": ["sql"], "properties": {
         "sql": {"type": "string"},
         "limit": {"type": "integer", "minimum": 1, "maximum": MAX_ROWS, "default": 200}}}},
]


def _rows(cols: list[str], rows: list[tuple]) -> list[dict]:
    return [dict(zip(cols, r, strict=True)) for r in rows]


def call(releases_dir: Path, name: str, args: dict) -> object:
    if name == "guide":
        return GUIDE.read_text() if GUIDE.exists() else INSTRUCTIONS
    r = open_release(releases_dir)
    try:
        if name == "overview":
            m = r.manifest
            checks = (_rows(*select(r.con, "SELECT counts FROM doc_checks", 1))[0]["counts"]
                      if "checks" in r.docs else None)
            return {"release": r.id, "created": m.get("created"), "pipeline_version": m.get("pipeline_version"),
                    "truth_set": m.get("truth_set"), "knowledge": m.get("knowledge"), "self_checks": checks,
                    "tables": {k: {"rows": t.get("rows"), "description": t.get("description")}
                               for k, t in m["tables"].items()},
                    "documents": [f"doc_{d}" for d in r.docs]}
        if name == "highlights":
            if "highlights" not in r.manifest["tables"]:
                return []
            return _rows(*select(r.con, "SELECT rank, kind_label, section, title, statement, notable_because, "
                                        "overall, evidence_level, call_level, sensitive, ref_table, ref_key "
                                        "FROM highlights ORDER BY rank", MAX_ROWS))
        if name == "result":
            table, key = args["table"], str(args["key"])
            if table not in CITABLE or table not in r.manifest["tables"]:
                raise ValueError(f"no table {table} in this release")
            col = CITABLE[table][0]
            rows = r.con.execute(f"SELECT * FROM {table} WHERE CAST({col} AS VARCHAR) = ?", [key])
            found = _rows([d[0] for d in rows.description], rows.fetchall())
            if not found:
                raise ValueError(f"no {table} row with {col} = {key}")
            return found
        if name == "query":
            limit = max(1, min(int(args.get("limit", 200)), MAX_ROWS))
            return _rows(*select(r.con, args["sql"], limit))
    finally:
        r.con.close()
    raise ValueError(f"unknown tool {name}")


def handle(releases_dir: Path, msg: dict) -> dict | None:
    method, mid = msg.get("method"), msg.get("id")
    if mid is None:  # a notification (initialized, cancelled, …): no reply
        return None
    if method == "initialize":
        asked = msg.get("params", {}).get("protocolVersion")
        result = {"protocolVersion": asked if asked in PROTOCOLS else PROTOCOLS[0],
                  "capabilities": {"tools": {}}, "serverInfo": {"name": "wgs", "version": __version__},
                  "instructions": INSTRUCTIONS}
    elif method == "ping":
        result = {}
    elif method == "tools/list":
        result = {"tools": TOOLS}
    elif method == "tools/call":
        p = msg.get("params", {})
        try:
            out = call(releases_dir, p.get("name", ""), p.get("arguments") or {})
            text = out if isinstance(out, str) else json.dumps(out, default=str, indent=1)
            result = {"content": [{"type": "text", "text": text}], "isError": False}
        except Exception as e:  # a bad query is the assistant's to fix: report it as a tool error
            result = {"content": [{"type": "text", "text": f"{type(e).__name__}: {e}"}], "isError": True}
    else:
        return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32601, "message": f"unknown method {method}"}}
    return {"jsonrpc": "2.0", "id": mid, "result": result}


def serve(releases_dir: Path, stdin: IO[str] = sys.stdin, stdout: IO[str] = sys.stdout) -> None:
    for line in stdin:
        if not line.strip():
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            reply = {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "parse error"}}
        else:
            reply = handle(releases_dir, msg) if isinstance(msg, dict) else {
                "jsonrpc": "2.0", "id": None, "error": {"code": -32600, "message": "batches are not supported"}}
        if reply is not None:
            stdout.write(json.dumps(reply, default=str) + "\n")
            stdout.flush()
