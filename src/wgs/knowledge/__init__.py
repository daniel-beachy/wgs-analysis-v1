"""Knowledge layer: versioned, local copies of public genomic knowledge (ADR-012).

Privacy rule: knowledge comes to the data, never the other way round. Sources are downloaded in
full and joined locally; no genotype or variant list is ever sent to a web service.

Layout (inside the workspace)::

    knowledge/index.json                    which versions exist and which one is current, per source
    knowledge/<source>/<version>/*.parquet  normalised tables, ready for DuckDB
    knowledge/<source>/<version>/source.json  provenance: URLs, upstream dates, sizes, sha256, licence

A refresh probes each source upstream (HTTP Last-Modified / ETag / dated filenames), downloads only
what changed, builds a new version directory and makes it current. Old versions are kept (a few,
configurable) so releases built from them stay explainable and diffs can name the version change.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from ..config import Config
from ..pipeline import console, now, write_json


@dataclass
class Upstream:
    url: str
    last_modified: str | None = None
    etag: str | None = None
    size: int | None = None

    def fingerprint(self) -> str:
        return f"{self.url}|{self.last_modified}|{self.etag}|{self.size}"


@dataclass
class Built:
    version: str                       # human-meaningful, e.g. "2026-09-28" or "r75"
    tables: dict[str, Path]            # name -> parquet written in the version dir
    upstream: list[Upstream] = field(default_factory=list)
    notes: dict = field(default_factory=dict)


class Source:
    """A public knowledge source. Subclasses implement probe() and build()."""

    id: str = ""
    title: str = ""
    homepage: str = ""
    licence: str = ""
    cadence: str = ""                  # how often upstream changes, for display
    description: str = ""
    schema: int = 1                    # bump when build() output changes, so stored versions are rebuilt
    workspace: Path = Path(".")        # set by refresh(); lets a build reuse files in <workspace>/cache

    def probe(self, pin: str | None = None) -> list[Upstream]:
        raise NotImplementedError

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        raise NotImplementedError


class Store:
    def __init__(self, cfg: Config):
        self.root = cfg.workspace / "knowledge"
        self.root.mkdir(parents=True, exist_ok=True)
        self.index_file = self.root / "index.json"

    def index(self) -> dict:
        if self.index_file.exists():
            return json.loads(self.index_file.read_text())
        return {"schema": 1, "sources": {}}

    def current(self, source_id: str) -> dict | None:
        s = self.index()["sources"].get(source_id)
        if not s or not s.get("current"):
            return None
        return next((v for v in s["versions"] if v["version"] == s["current"]), None)

    def table(self, source_id: str, name: str) -> Path | None:
        cur = self.current(source_id)
        if not cur or name not in cur["tables"]:
            return None
        p = self.root / cur["tables"][name]["path"]
        return p if p.exists() else None

    def file(self, source_id: str, name: str) -> Path | None:
        """A non-table file stored alongside a source's current version (e.g. PharmCAT's jar)."""
        cur = self.current(source_id)
        if not cur:
            return None
        p = self.root / source_id / cur["version"] / name
        return p if p.exists() else None

    def versions(self) -> dict[str, str]:
        """{source: current version} — recorded in every release so results are traceable."""
        return {k: v["current"] for k, v in self.index()["sources"].items() if v.get("current")}

    def record(self, src: Source, built: Built, vdir: Path, keep: int) -> None:
        idx = self.index()
        entry = idx["sources"].setdefault(src.id, {"versions": []})
        entry.update({"title": src.title, "homepage": src.homepage, "licence": src.licence, "cadence": src.cadence,
                      "description": src.description})
        tables = {}
        for name, p in built.tables.items():
            import duckdb

            rows = duckdb.connect().execute(f"SELECT count(*) FROM '{p}'").fetchone()[0]
            tables[name] = {"path": str(p.relative_to(self.root)), "rows": rows, "bytes": p.stat().st_size,
                            "sha256": _sha256(p)}
        version = {"version": built.version, "built": now(), "schema": src.schema, "tables": tables,
                   "upstream": [u.__dict__ for u in built.upstream], "notes": built.notes}
        write_json(vdir / "source.json", {"source": src.id, **version})
        entry["versions"] = [v for v in entry["versions"] if v["version"] != built.version] + [version]
        entry["versions"].sort(key=lambda v: v["version"], reverse=True)
        entry["current"] = built.version
        entry["checked"] = now()
        cur = next(v for v in entry["versions"] if v["version"] == built.version)
        kept = [cur] + [v for v in entry["versions"] if v is not cur][:max(keep - 1, 0)]
        for old in entry["versions"]:
            if old not in kept:
                shutil.rmtree(self.root / src.id / old["version"], ignore_errors=True)
        entry["versions"] = sorted(kept, key=lambda v: v["version"], reverse=True)
        write_json(self.index_file, idx)

    def kept_matching(self, source_id: str, fp: list[str], schema: int = 1) -> dict | None:
        """A stored version whose upstream matches `fp` and whose files are intact (lets un-pinning skip downloads)."""
        for v in self.index()["sources"].get(source_id, {}).get("versions", []):
            if _fingerprint(v) == fp and v.get("schema", 1) == schema and \
                    all((self.root / t["path"]).exists() for t in v["tables"].values()):
                return v
        return None

    def activate(self, source_id: str, version: str) -> None:
        idx = self.index()
        idx["sources"][source_id]["current"] = version
        idx["sources"][source_id]["checked"] = now()
        write_json(self.index_file, idx)

    def mark_checked(self, source_id: str) -> None:
        idx = self.index()
        if source_id in idx["sources"]:
            idx["sources"][source_id]["checked"] = now()
            write_json(self.index_file, idx)


def _fingerprint(version: dict) -> list[str]:
    return [f"{u['url']}|{u.get('last_modified')}|{u.get('etag')}|{u.get('size')}" for u in version["upstream"]]


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def refresh(cfg: Config, scratch: Path, only: list[str] | None = None, pin: dict[str, str] | None = None,
            force: bool = False, keep: int = 4) -> dict[str, str]:
    """Bring sources up to date. Returns {source: 'updated <v>' | 'up to date <v>' | 'failed: ...'}."""
    from .registry import SOURCES

    store = Store(cfg)
    out: dict[str, str] = {}
    for src in SOURCES:
        if only and src.id not in only:
            continue
        p = (pin or {}).get(src.id)
        src.workspace = cfg.workspace
        try:
            up = src.probe(p)
            cur = store.current(src.id)
            fp = [u.fingerprint() for u in up]
            same = cur and _fingerprint(cur) == fp and cur.get("schema", 1) == src.schema
            kept = None if same or force else store.kept_matching(src.id, fp, src.schema)
            if kept:
                store.activate(src.id, kept["version"])
                out[src.id] = f"updated ({kept['version']})"
                console.print(f"[green]✓[/] {src.title}: now {kept['version']} (already stored, no download)")
                continue
            if same and not force:
                store.mark_checked(src.id)
                out[src.id] = f"up to date ({cur['version']})"
                console.print(f"[green]✓[/] {src.title}: up to date ({cur['version']})")
                continue
            console.print(f"[cyan]↓[/] {src.title}: fetching {', '.join(u.url.rsplit('/', 1)[-1] for u in up)}")
            work = scratch / "knowledge" / src.id
            shutil.rmtree(work, ignore_errors=True)
            work.mkdir(parents=True)
            staging = store.root / src.id / ".building"
            shutil.rmtree(staging, ignore_errors=True)
            staging.mkdir(parents=True)
            built = src.build(up, work, staging)
            vdir = store.root / src.id / built.version
            shutil.rmtree(vdir, ignore_errors=True)
            staging.rename(vdir)
            built.tables = {k: vdir / v.name for k, v in built.tables.items()}
            store.record(src, built, vdir, keep)
            shutil.rmtree(work, ignore_errors=True)
            out[src.id] = f"updated ({built.version})"
            console.print(f"[green]✓[/] {src.title}: now {built.version}")
        except Exception as exc:  # one failing source must not block the others
            out[src.id] = f"failed: {exc}"
            console.print(f"[red]✗ {src.title}: {exc}[/]")
    return out
