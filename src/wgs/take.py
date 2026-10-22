"""Claude's take (Step 11): an AI-written reading of the highlights, published as a clearly labelled opinion.

The take is the one place in a release where the words are not produced by a rule (ADR-026). It is written in a
Claude Code session into the private workspace (`<workspace>/notes/take.md`), never into this public repo, and the
pipeline holds it to three rules before it is shown:

- every paragraph and list cites at least one result, as `[[table:key]]` (e.g. `[[highlights:<id>]]`,
  `[[claims:<claim_id>]]`, `[[pgx_genes:<gene>]]`);
- every citation resolves to a row of this release; the dashboard shows each one with that row's current grade;
- the front matter names the release it was written for (`written_for`); if a cited row's grade or statement has
  changed since, or the row is gone, the take is marked stale.

The markdown is parsed here into blocks and inline segments, so the dashboard never renders raw HTML.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from pathlib import Path

import duckdb

CITE = re.compile(r"\[\[([a-z_]+):([^\]]+)\]\]")
INLINE = re.compile(r"\[\[([a-z_]+):([^\]]+)\]\]|\*\*(.+?)\*\*|\*(.+?)\*")

# table -> (key column, title expression, columns whose change makes a take stale, page that shows it)
CITABLE = {
    "highlights": ("highlight_id", "title", ["overall", "statement"], "section"),
    "claims": ("claim_id", "coalesce(gene, subject)", ["overall", "statement"], "section"),
    "pgx_genes": ("gene", "gene || ' ' || coalesce(diplotype, '')", ["phenotype", "overall"], "'pgx'"),
    "pgs_traits": ("trait_id", "label", ["CAST(round(percentile) AS INTEGER)", "overall"], "section"),
    "trait_snps": ("rsid", "label", ["genotype", "overall"], "'traits'"),
    "hla_associations": ("association_id", "condition", ["notable", "overall"], "'health'"),
    "special_genes": ("gene", "gene", ["status", "overall"], "'carrier'"),
    "repeats": ("locus_id", "coalesce(gene, locus_id)", ["status", "overall"], "'health'"),
    "cnv_events": ("event_id", "array_to_string(genes, ', ')", ["status", "overall"], "'health'"),
    "ancestry_components": ("superpop", "name", ["detected"], "'ancestry'"),
    "haplogroups": ("line", "haplogroup", ["haplogroup"], "'ancestry'"),
}


def path(workspace: Path) -> Path:
    return workspace / "notes" / "take.md"


def _front(text: str) -> tuple[dict, str]:
    if not text.startswith("---\n"):
        return {}, text
    head, _, body = text[4:].partition("\n---\n")
    meta = {}
    for line in head.splitlines():
        k, sep, v = line.partition(":")
        if sep:
            meta[k.strip()] = v.strip()
    return meta, body


def inline(text: str) -> list[dict]:
    """Text, bold, italic and citation segments: the only inline formatting a take may use."""
    out, at = [], 0
    for m in INLINE.finditer(text):
        if m.start() > at:
            out.append({"t": "text", "v": text[at:m.start()]})
        if m.group(1):
            out.append({"t": "cite", "ref": f"{m.group(1)}:{m.group(2)}"})
        else:
            out.append({"t": "b" if m.group(3) else "i", "v": m.group(3) or m.group(4)})
        at = m.end()
    if at < len(text):
        out.append({"t": "text", "v": text[at:]})
    return out


def blocks(body: str) -> list[dict]:
    """Headings (`#`, `##`, `###`), bullet lists (`- `) and paragraphs separated by blank lines."""
    out = []
    for chunk in re.split(r"\n\s*\n", body.strip()):
        lines = [x.strip() for x in chunk.strip().splitlines() if x.strip()]
        if not lines:
            continue
        if lines[0].startswith("#"):
            out.append({"type": "h", "text": lines[0].lstrip("#").strip()})
            lines = lines[1:]
            if not lines:
                continue
        if all(x.startswith("- ") for x in lines):
            out.append({"type": "ul", "raw": "\n".join(lines), "items": [inline(x[2:]) for x in lines]})
        else:
            raw = " ".join(lines)
            out.append({"type": "p", "raw": raw, "segments": inline(raw)})
    return out


def _rows(con, table: Path, spec: tuple, keys: list[str]) -> dict[str, dict]:
    key, title, fields, page = spec
    cols = ", ".join(f"CAST({f} AS VARCHAR)" for f in fields)
    rows = con.execute(f"SELECT CAST({key} AS VARCHAR), {title}, {page}, list_value({cols}), "
                       f"{'overall' if 'overall' in fields else 'NULL'} FROM '{table}' "
                       f"WHERE CAST({key} AS VARCHAR) IN (SELECT unnest(?))", [keys]).fetchall()
    return {r[0]: {"title": r[1], "page": r[2], "fingerprint": r[3], "overall": r[4]} for r in rows}


def build(file: Path, tables: dict[str, Path],
          tables_of: Callable[[str | None], dict[str, Path] | None]) -> dict | None:
    """The take document for a release, or None when no take has been written. `tables_of(release_id)` gives
    the tables of the release the take was written for (None if that release is unknown)."""
    if not file.exists():
        return None
    meta, body = _front(file.read_text())
    written_for_tables = tables_of(meta.get("written_for"))
    bl = blocks(body)
    problems = [f"Cites nothing: “{b['raw'][:80]}…”" for b in bl if b["type"] != "h" and not CITE.search(b["raw"])]
    refs = sorted({f"{t}:{k}" for b in bl if b["type"] != "h" for t, k in CITE.findall(b["raw"])})
    con = duckdb.connect()
    found: dict[str, dict] = {}
    before: dict[str, dict] = {}
    for t in sorted({r.split(":", 1)[0] for r in refs}):
        keys = [r.split(":", 1)[1] for r in refs if r.split(":", 1)[0] == t]
        if t not in CITABLE:
            problems += [f"Can't cite table {t}: {t}:{k}" for k in keys]
            continue
        if t in tables:
            found |= {f"{t}:{k}": v for k, v in _rows(con, tables[t], CITABLE[t], keys).items()}
        if written_for_tables and t in written_for_tables:
            before |= {f"{t}:{k}": v for k, v in _rows(con, written_for_tables[t], CITABLE[t], keys).items()}
    con.close()
    citations = []
    for r in refs:
        t, k = r.split(":", 1)
        row = found.get(r)
        if row is None and t in CITABLE:
            problems.append(f"Not in this release: {r}")
        changed = (written_for_tables is not None and t in CITABLE
                   and (before.get(r) or {}).get("fingerprint") != (row or {}).get("fingerprint"))
        citations.append({"ref": r, "table": t, "key": k, "found": row is not None, "changed": changed,
                          **({k2: row[k2] for k2 in ("title", "page", "overall")} if row else {})})
    for b in bl:
        b.pop("raw", None)
    return {"author": meta.get("author", "Claude"), "written": meta.get("written"),
            "written_for": meta.get("written_for"), "blocks": bl, "citations": citations, "problems": problems,
            "stale": any(c["changed"] for c in citations) or written_for_tables is None}


def to_json(doc: dict) -> str:
    return json.dumps(doc, indent=1, ensure_ascii=False)
