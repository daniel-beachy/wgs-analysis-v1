"""Export the public demo release as static files for a web host (ADR-028).

Only a release built from a public truth set (its manifest has `truth_set`: the HG002 demo) can be exported; a
personal release never leaves the machine. The export copies the newest release into `<out>/data/`, the folder the
viewer reads, with these changes for the web:

* `annotations` loses the columns no page reads and the bulk predictor scores whose licences allow only
  non-commercial use (the few graded findings that quote a score keep it, with its source);
* tables larger than PART_BYTES are rewritten with Parquet V2 encodings, in position order (`sorted_by`), and split by
  chromosome into parts, so no file passes MAX_FILE_BYTES (Cloudflare Pages allows 25 MiB per file, GitHub 100 MB);
* `index.json` lists only that release, and its manifest records what changed under `export`;
* absolute paths on this machine inside JSON documents are cut down to the file name.

Every output file is then checked for the private strings of the personal config (sample and kit IDs).
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path

import duckdb

PART_BYTES = 20_000_000
MAX_FILE_BYTES = 25_000_000
LOCAL_PATH = re.compile(r"^(/Volumes/|/Users/|/home/|/private/|/var/folders/|/tmp/|[A-Za-z]:\\)")

DROP = {
    "annotations": {
        "revel": "REVEL is free for non-commercial use only, so its scores are not republished in bulk",
        "am_score": "AlphaMissense is CC BY-NC-SA 4.0, so its scores are not republished in bulk",
        "gnomad_popmax_af_controls": "not shown on any page",
        "gnomad_nhomalt": "not shown on any page",
        "gnomad_nhomalt_controls": "not shown on any page",
    },
}


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _float_column(name: str) -> bool:
    """Allele frequencies need no more than single precision, which halves their size."""
    return name.startswith("af_") or name.endswith("_af")


def _rewrite(src: Path, name: str, drop: dict[str, str], objects: Path) -> tuple[list[dict], bool]:
    """Write one table's parts (Parquet V2, sorted by position, split by chromosome); returns their refs and
    whether the rows are in position order (the viewer then skips ORDER BY and stops at the first page)."""
    con = duckdb.connect()
    cols = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{src.as_posix()}')").fetchall()
    names = [c[0] for c in cols]
    sel = ", ".join(f'CAST("{n}" AS FLOAT) AS "{n}"' if t == "DOUBLE" and _float_column(n) else f'"{n}"'
                    for n, t, *_ in cols if n not in drop)
    positional = "chrom_order" in names and "pos" in names
    order = " ORDER BY chrom_order, pos" if positional else ""
    opts = "FORMAT parquet, COMPRESSION zstd, COMPRESSION_LEVEL 9, ROW_GROUP_SIZE 65536, PARQUET_VERSION V2"

    def write(where: str, dest: Path) -> int:
        con.execute(f"COPY (SELECT {sel} FROM read_parquet('{src.as_posix()}') {where}{order}) TO "
                    f"'{dest.as_posix()}' ({opts})")
        return con.execute(f"SELECT count(*) FROM read_parquet('{dest.as_posix()}')").fetchone()[0]

    whole = objects / f"{name}.tmp.parquet"
    rows = write("", whole)
    size = whole.stat().st_size
    groups: list[tuple[str, int]] = [("", rows)]
    if size > PART_BYTES and positional:
        # Greedy: consecutive chromosomes share a part while the estimated size stays under PART_BYTES.
        per_row = size / max(rows, 1)
        counts = con.execute(f"SELECT chrom_order, count(*) FROM read_parquet('{whole.as_posix()}') "
                             "GROUP BY 1 ORDER BY 1").fetchall()
        groups, lo, hi, n = [], None, None, 0
        for co, c in counts:
            if lo is not None and (n + c) * per_row > PART_BYTES:
                groups.append((f"WHERE chrom_order BETWEEN {lo} AND {hi}", n))
                lo, n = None, 0
            lo = co if lo is None else lo
            hi, n = co, n + c
        groups.append((f"WHERE chrom_order BETWEEN {lo} AND {hi}", n))
    parts = []
    for i, (where, _) in enumerate(groups):
        tmp = objects / f"{name}.part{i}.tmp.parquet"
        if where:
            n = write(where, tmp)
        else:
            whole.replace(tmp)
            n = rows
        sha = _sha(tmp)
        dest = objects / (f"{name}-{i + 1}-{sha[:12]}.parquet" if len(groups) > 1 else f"{name}-{sha[:12]}.parquet")
        tmp.replace(dest)
        parts.append({"path": f"objects/{dest.name}", "bytes": dest.stat().st_size, "sha256": sha, "rows": n})
    whole.unlink(missing_ok=True)
    return parts, positional


def _strip_paths(o):
    """Cut absolute paths on this machine (in QC documents, for example) down to the file name."""
    if isinstance(o, dict):
        return {k: _strip_paths(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_strip_paths(v) for v in o]
    if isinstance(o, str) and LOCAL_PATH.match(o):
        return re.split(r"[/\\]", o)[-1]
    return o


def _forbidden_hits(out: Path, forbidden: list[str]) -> list[str]:
    """Files whose text (JSON) or string values (Parquet) contain a private string."""
    words = [w for w in forbidden if w]
    if not words:
        return []
    hits = []
    con = duckdb.connect()
    for f in sorted(out.rglob("*")):
        if not f.is_file():
            continue
        if f.suffix == ".parquet":
            cols = [c[0] for c in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{f.as_posix()}')").fetchall()]
            cond = " OR ".join(f'contains(CAST("{c}" AS VARCHAR), $w)' for c in cols)
            if any(con.execute(f"SELECT count(*) FROM read_parquet('{f.as_posix()}') WHERE {cond}", {"w": w})
                   .fetchone()[0] for w in words):
                hits.append(str(f.relative_to(out)))
        elif any(w.encode() in f.read_bytes() for w in words):
            hits.append(str(f.relative_to(out)))
    return hits


def export(releases: Path, out: Path, forbidden: list[str] | None = None) -> dict:
    """Export the newest release in `releases` to `out/data`; returns a summary. Raises ValueError when refused."""
    idx = json.loads((releases / "index.json").read_text())
    entry = next(r for r in idx["releases"] if r["id"] == idx["latest"])
    manifest = json.loads((releases / entry["manifest"]).read_text())
    if not manifest.get("truth_set"):
        raise ValueError(f"Release {entry['id']} was not built from a public truth set; only the demo release "
                         "(config/demo.toml) can be exported.")
    for w in forbidden or []:
        if w and w == manifest.get("sample"):
            raise ValueError("The release's sample matches a private ID from the personal config.")

    data = out / "data"
    if data.exists():
        shutil.rmtree(data)
    objects = data / "objects"
    objects.mkdir(parents=True)

    split: dict[str, int] = {}
    for name, ref in manifest["tables"].items():
        src = releases / ref["path"]
        if name in DROP or ref["bytes"] > PART_BYTES:
            parts, positional = _rewrite(src, name, DROP.get(name, {}), objects)
            new = {**ref, "path": parts[0]["path"], "bytes": sum(p["bytes"] for p in parts),
                   "rows": sum(p["rows"] for p in parts)}
            if positional:
                new["sorted_by"] = ["chrom_order", "pos"]
            new["sha256"] = parts[0]["sha256"] if len(parts) == 1 else ""
            if len(parts) > 1:
                new["parts"] = parts
                split[name] = len(parts)
            manifest["tables"][name] = new
        else:
            shutil.copyfile(src, data / ref["path"])
    for ref in manifest["documents"].values():
        src, dest = releases / ref["path"], data / ref["path"]
        if src.suffix == ".json":
            doc = json.loads(src.read_text())
            if (clean := _strip_paths(doc)) != doc:
                dest.write_text(json.dumps(clean))
                continue
        shutil.copyfile(src, dest)

    manifest["export"] = {
        "created": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "from_release": entry["id"],
        "dropped_columns": {t: cols for t, cols in DROP.items() if t in manifest["tables"]},
        "split": split,
        "note": "Public copy of the demo release for a static web host (ADR-028). Large tables are split by "
                "chromosome; the columns listed were left out.",
    }
    rel_dir = data / entry["id"]
    rel_dir.mkdir()
    (rel_dir / "manifest.json").write_text(json.dumps(manifest, indent=1))
    (data / "index.json").write_text(json.dumps({**idx, "releases": [entry], "latest": entry["id"]}, indent=1))

    too_big = [str(f.relative_to(out)) for f in data.rglob("*") if f.is_file() and f.stat().st_size > MAX_FILE_BYTES]
    hits = _forbidden_hits(data, forbidden or [])
    if too_big or hits:
        shutil.rmtree(data)
        raise ValueError(f"Export removed: files over {MAX_FILE_BYTES} bytes {too_big}; private strings in {hits}")
    files = [f for f in data.rglob("*") if f.is_file()]
    return {"release": entry["id"], "files": len(files), "bytes": sum(f.stat().st_size for f in files),
            "largest": max(f.stat().st_size for f in files), "split": split}
