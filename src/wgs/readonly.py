"""Read-only SQL over one release, shared by `wgs query` and `wgs mcp`.

Every table is a view over the release's Parquet files and every JSON document a `doc_<name>` view. The connection
can read only the releases folder and nothing can be written; only a single SELECT (or WITH … SELECT) is accepted.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import duckdb


@dataclass
class Release:
    con: duckdb.DuckDBPyConnection
    id: str
    manifest: dict
    docs: list[str]


def open_release(releases_dir: Path, release: str | None = None) -> Release:
    index = json.loads((releases_dir / "index.json").read_text())
    rel = next((r for r in index["releases"] if r["id"] == release), None) if release else index["releases"][0]
    if rel is None:
        raise ValueError(f"no release {release}")
    m = json.loads((releases_dir / rel["manifest"]).read_text())
    con = duckdb.connect()
    for name, t in m["tables"].items():
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{(releases_dir / t['path']).as_posix()}')")
    docs = {name: releases_dir / (d["path"] if isinstance(d, dict) else d)
            for name, d in m.get("documents", {}).items()}
    docs = {name: path for name, path in docs.items() if path.suffix == ".json"}
    for name, path in docs.items():
        con.execute(f"CREATE VIEW doc_{name} AS SELECT * FROM read_json_auto('{path.as_posix()}')")
    # Read-only sandbox: only the releases folder is readable, nothing can be written.
    con.execute(f"SET allowed_directories = ['{releases_dir.as_posix()}/']")
    con.execute("SET enable_external_access = false")
    con.execute("SET lock_configuration = true")
    return Release(con, rel["id"], m, list(docs))


def select(con: duckdb.DuckDBPyConnection, sql: str, limit: int) -> tuple[list[str], list[tuple]]:
    stmts = con.extract_statements(sql)
    if len(stmts) != 1 or stmts[0].type != duckdb.StatementType.SELECT:
        raise ValueError("only a single SELECT (or WITH … SELECT) query is allowed")
    rel = con.sql(sql)
    return rel.columns, rel.fetchmany(limit)
