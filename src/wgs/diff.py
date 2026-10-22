"""What changed between two releases — the heart of "this is not a snapshot" (ADR-012, ADR-024).

Compares the previous release's claims and annotations with the new ones and explains each
difference in terms of *why*: a knowledge source moved to a new version, the evidence model
changed, or your own data/pipeline changed. Output is a small JSON document stored in the release.
"""

from __future__ import annotations

import json
from pathlib import Path

import duckdb

CLAIM_FIELDS = ["category", "clinvar_stars", "evidence_level", "call_level", "overall", "gene_validity",
                "section", "statement"]
LIMIT = 500

# Results that are not graded claims but can still move when knowledge or the pipeline changes (Step 10). For each:
# the key that identifies a row across releases, the column that names it, the fields compared (with the smallest
# change that counts for numbers, so float noise is not news) and the knowledge sources the table is built from.
RESULTS = [
    {"table": "pgx_genes", "label": "Medicine gene", "key": "gene", "name": "gene",
     "fields": {"diplotype": None, "phenotype": None, "activity_score": 0.01}, "sources": ["pharmcat"]},
    {"table": "pgs_traits", "label": "Polygenic score", "key": "trait_id", "name": "label",
     "fields": {"percentile": 1.0, "agreement": None, "overall": None},
     "sources": ["pgs_catalog", "pgs_reference", "1000genomes"]},
    {"table": "ancestry_components", "label": "Ancestry", "key": "superpop", "name": "name",
     "fields": {"share": 0.01, "detected": None}, "sources": ["1000genomes", "igsr_populations"]},
    {"table": "haplogroups", "label": "Haplogroup", "key": "line", "name": "line", "fields": {"haplogroup": None},
     "sources": []},
    {"table": "special_genes", "label": "Hard-to-read gene", "key": "gene", "name": "gene",
     "fields": {"status": None, "overall": None}, "sources": []},
    {"table": "repeats", "label": "Repeat", "key": "locus_id", "name": "gene",
     "fields": {"status": None, "overall": None}, "sources": ["strchive"]},
]


def _cols(con, p: Path) -> set[str]:
    return {r[0] for r in con.execute(f"DESCRIBE SELECT * FROM '{p}'").fetchall()}


def diff(prev_manifest: dict | None, prev_root: Path, new_tables: dict[str, Path], new_knowledge: dict,
         new_model: str | None) -> dict:
    out: dict = {"schema": 2, "previous": None, "knowledge": {}, "evidence_model": None,
                 "claims": {"added": [], "removed": [], "changed": [], "counts": {}},
                 "results": {"items": [], "count": 0},
                 "reclassified": {"items": [], "count": 0}, "first": prev_manifest is None}
    if prev_manifest is None:
        return out
    out["previous"] = prev_manifest.get("id")
    old_k = prev_manifest.get("knowledge", {}) or {}
    old_versions = old_k.get("versions", old_k) if isinstance(old_k, dict) else {}
    for src in sorted(set(old_versions) | set(new_knowledge)):
        a, b = old_versions.get(src), new_knowledge.get(src)
        if a != b:
            out["knowledge"][src] = [a, b]
    old_model = old_k.get("evidence_model") if isinstance(old_k, dict) else None
    if old_model != new_model:
        out["evidence_model"] = [old_model, new_model]

    con = duckdb.connect()
    old_t = prev_manifest.get("tables", {})
    out["tables"] = {"added": sorted(set(new_tables) - set(old_t)), "removed": sorted(set(old_t) - set(new_tables))}
    oc = prev_root / old_t["claims"]["path"] if "claims" in old_t else None
    nc = new_tables.get("claims")
    if nc and nc.exists() and "claim_id" in _cols(con, nc):
        if oc and oc.exists() and "claim_id" in _cols(con, oc):
            fields = [f for f in CLAIM_FIELDS if f in _cols(con, oc) and f in _cols(con, nc)]
            sens = "coalesce(sensitive, false)" if "sensitive" in _cols(con, nc) else "false"
            src = lambda p: "CAST(sources AS VARCHAR)" if "sources" in _cols(con, p) else "NULL"  # noqa: E731
            base = (f"claim_id, section, subject, gene, rsid, category_label, overall, statement, {sens} AS sensitive, "
                    f"{src(nc)} AS _sources")
            old_base = base if "sensitive" in _cols(con, oc) else base.replace(sens, "false")
            old_base = old_base.replace(f"{src(nc)} AS _sources", f"{src(oc)} AS _sources")
            out["claims"]["added"] = _rows(con, f"""
              SELECT {base} FROM '{nc}' WHERE claim_id NOT IN (SELECT claim_id FROM '{oc}')
              ORDER BY section, claim_id LIMIT {LIMIT}""")
            out["claims"]["removed"] = _rows(con, f"""
              SELECT {old_base} FROM '{oc}' WHERE claim_id NOT IN (SELECT claim_id FROM '{nc}')
              ORDER BY section, claim_id LIMIT {LIMIT}""")
            diffs = " OR ".join(f"o.{f} IS DISTINCT FROM n.{f}" for f in fields)
            # "regraded": grade, stars, level or section moved; "reworded": only the explanatory text changed.
            graded = " OR ".join(f"o.{f} IS DISTINCT FROM n.{f}" for f in fields if f != "statement") or "false"
            rows = _rows(con, f"""
              SELECT n.claim_id, n.section, n.subject, n.gene, n.rsid, n.category_label, n.overall,
                     {"coalesce(n.sensitive, false)" if "sensitive" in _cols(con, nc) else "false"} AS sensitive,
                     CASE WHEN {graded} THEN 'regraded' ELSE 'reworded' END AS change_kind,
                     {"CAST(n.sources AS VARCHAR)" if "sources" in _cols(con, nc) else "NULL"} AS _sources,
                     {', '.join(f'o.{f} AS old_{f}, n.{f} AS new_{f}' for f in fields)}
              FROM '{oc}' o JOIN '{nc}' n USING (claim_id) WHERE {diffs}
              QUALIFY row_number() OVER (PARTITION BY change_kind ORDER BY n.section, n.claim_id) <= {LIMIT}
              ORDER BY change_kind, n.section, n.claim_id""")
            for r in rows:
                r["fields"] = {f: [r.pop(f"old_{f}"), r.pop(f"new_{f}")] for f in fields
                               if r[f"old_{f}"] != r[f"new_{f}"]}
                r.pop("statement", None)
                r["kind"] = r.pop("change_kind")
            out["claims"]["changed"] = rows
            for r in out["claims"]["added"] + out["claims"]["removed"] + rows:
                r["why"] = why(_source_ids(r.pop("_sources", None)), out)
            q = {"added": f"FROM '{nc}' WHERE claim_id NOT IN (SELECT claim_id FROM '{oc}')",
                 "removed": f"FROM '{oc}' WHERE claim_id NOT IN (SELECT claim_id FROM '{nc}')",
                 "changed": f"FROM '{oc}' o JOIN '{nc}' n USING (claim_id) WHERE {diffs}",
                 "regraded": f"FROM '{oc}' o JOIN '{nc}' n USING (claim_id) WHERE {graded}"}
            for kind, frm in q.items():
                out["claims"]["counts"][kind] = con.execute(f"SELECT count(*) {frm}").fetchone()[0]
        else:
            out["claims"]["counts"]["added"] = con.execute(f"SELECT count(*) FROM '{nc}'").fetchone()[0]
            out["claims"]["note"] = "Previous release had no graded claims; all current claims are new."

    out["results"] = results(con, old_t, prev_root, new_tables, out)

    oa = prev_root / old_t["annotations"]["path"] if "annotations" in old_t else None
    na = new_tables.get("annotations")
    if oa and oa.exists() and na and na.exists():
        joined = f"""
          FROM (SELECT * FROM '{oa}' WHERE clinvar_id IS NOT NULL) o
          FULL JOIN (SELECT * FROM '{na}' WHERE clinvar_id IS NOT NULL) n USING (chrom, pos, ref, alt)"""
        # Only a change of class is a reclassification; same class with a new review status or wording is
        # counted separately so routine ClinVar churn does not drown out the changes that matter.
        sql = f"{joined} WHERE o.clinvar_class IS DISTINCT FROM n.clinvar_class"
        out["reclassified"]["count"] = con.execute(f"SELECT count(*) {sql}").fetchone()[0]
        out["reclassified"]["review_changed"] = con.execute(f"""SELECT count(*) {joined}
          WHERE o.clinvar_class = n.clinvar_class AND o.clinvar_stars IS DISTINCT FROM n.clinvar_stars""").fetchone()[0]
        out["reclassified"]["detail_changed"] = con.execute(f"""SELECT count(*) {joined}
          WHERE o.clinvar_class = n.clinvar_class AND o.clinvar_stars IS NOT DISTINCT FROM n.clinvar_stars
            AND o.clinvar_significance IS DISTINCT FROM n.clinvar_significance""").fetchone()[0]
        out["reclassified"]["by_change"] = _rows(con, f"""
          SELECT coalesce(o.clinvar_class, 'not in ClinVar') AS old, coalesce(n.clinvar_class, 'not in ClinVar') AS new,
                 count(*) AS n {sql} GROUP BY ALL ORDER BY n DESC""")
        out["reclassified"]["items"] = _rows(con, f"""
          SELECT chrom, pos, ref, alt, coalesce(n.rsid, o.rsid) AS rsid, coalesce(n.gene, o.gene) AS gene,
                 o.clinvar_class AS old_class, n.clinvar_class AS new_class,
                 o.clinvar_stars AS old_stars, n.clinvar_stars AS new_stars,
                 coalesce(n.clinvar_id, o.clinvar_id) AS clinvar_id
          {sql}
          ORDER BY (coalesce(n.clinvar_class, '') IN ('pathogenic', 'likely_pathogenic')
                    OR coalesce(o.clinvar_class, '') IN ('pathogenic', 'likely_pathogenic')) DESC, chrom, pos
          LIMIT {LIMIT}""")
    con.close()
    return out


def _source_ids(text: str | None) -> set[str]:
    """Knowledge-source ids named in a claim's `sources` JSON (curated quotes name papers, which never match)."""
    try:
        return {s.get("source") for s in json.loads(text or "[]") if isinstance(s, dict)} - {None}
    except (TypeError, ValueError):
        return set()


def why(sources: set[str] | list[str], ch: dict) -> list[str]:
    """Why a result may have moved: the knowledge sources it rests on that changed version, plus the evidence model
    if the grading rules changed. Empty means neither, so the analysis itself changed (code or settings)."""
    out = sorted(s for s in sources if s in ch["knowledge"])
    if ch.get("evidence_model"):
        out.append("evidence_model")
    return out


def results(con, old_t: dict, prev_root: Path, new_tables: dict[str, Path], ch: dict) -> dict:
    """Added, removed or changed rows in the result tables that are not claims (RESULTS)."""
    items, counts = [], {}
    for spec in RESULTS:
        t = spec["table"]
        o = prev_root / old_t[t]["path"] if t in old_t else None
        n = new_tables.get(t)
        if not (o and o.exists() and n and n.exists()):
            continue
        oc, nc = _cols(con, o), _cols(con, n)
        key, name = spec["key"], spec["name"]
        if key not in oc or key not in nc:
            continue
        fields = [f for f in spec["fields"] if f in oc and f in nc]
        nm = name if name in oc and name in nc else key
        sens = "coalesce(n.sensitive, false)" if "sensitive" in nc else "false"

        def moved(f: str, tol: float | None) -> str:
            if tol is None:
                return f"o.{f} IS DISTINCT FROM n.{f}"
            return (f"(o.{f} IS NULL) <> (n.{f} IS NULL) OR abs(TRY_CAST(o.{f} AS DOUBLE) - TRY_CAST(n.{f} AS DOUBLE)) "
                    f">= {tol}")
        cond = " OR ".join(f"({moved(f, spec['fields'][f])})" for f in fields) or "false"
        rows = _rows(con, f"""
          SELECT coalesce(n.{key}, o.{key}) AS key, coalesce(n.{nm}, o.{nm}) AS name, {sens} AS sensitive,
                 CASE WHEN o.{key} IS NULL THEN 'added' WHEN n.{key} IS NULL THEN 'removed' ELSE 'changed' END AS kind,
                 {', '.join(f'o.{f} AS old_{f}, n.{f} AS new_{f}' for f in fields) or 'NULL AS _none'}
          FROM '{o}' o FULL JOIN '{n}' n USING ({key})
          WHERE o.{key} IS NULL OR n.{key} IS NULL OR {cond}
          ORDER BY kind, key LIMIT {LIMIT}""")
        for r in rows:
            r["table"], r["label"] = t, spec["label"]
            r["fields"] = {f: [r.pop(f"old_{f}"), r.pop(f"new_{f}")] for f in fields}
            if r["kind"] == "changed":
                r["fields"] = {f: v for f, v in r["fields"].items() if v[0] != v[1]}
            r.pop("_none", None)
            r["why"] = why(spec["sources"], ch)
        counts[t] = len(rows)
        items += rows
    return {"items": items, "count": len(items), "counts": counts}


def _rows(con, sql: str) -> list[dict]:
    cur = con.execute(sql)
    names = [d[0] for d in cur.description]
    return [dict(zip(names, (_plain(v) for v in row), strict=True)) for row in cur.fetchall()]


def _plain(v):
    if isinstance(v, float) and v != v:
        return None
    if hasattr(v, "item"):
        return v.item()
    return v


def summarise(ch: dict) -> str:
    if ch.get("first"):
        return "First release."
    parts = []
    if ch["knowledge"]:
        parts.append(", ".join(f"{k} added ({b})" if a is None else f"{k} removed" if b is None else f"{k} {a} → {b}"
                               for k, (a, b) in ch["knowledge"].items()))
    c = ch["claims"]["counts"]
    if any(c.get(k) for k in ("added", "removed", "changed")):
        txt = f"findings +{c.get('added', 0)} / −{c.get('removed', 0)}"
        if not c.get("removed") and not c.get("changed"):
            txt = f"findings +{c.get('added', 0)}"
        elif "regraded" in c:
            txt += f", {c['regraded']} regraded, {c.get('changed', 0) - c['regraded']} reworded"
        else:
            txt += f" / ~{c.get('changed', 0)}"
        parts.append(txt)
    if ch.get("results", {}).get("count"):
        parts.append(f"{ch['results']['count']} other results changed")
    if ch["reclassified"]["count"]:
        parts.append(f"{ch['reclassified']['count']} of your ClinVar variants reclassified")
    t = ch.get("tables", {})
    if t.get("added"):
        parts.append("new: " + ", ".join(t["added"]))
    if t.get("removed"):
        parts.append("no longer included: " + ", ".join(t["removed"]))
    return "; ".join(parts) or "No changes in findings."
