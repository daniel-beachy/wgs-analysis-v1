"""Self-checks run on every release, before it is published.

Each check is a SQL query over the release's own tables that returns the rows that break a rule; an empty result is
a pass. They catch the kind of bug a human would only notice by reading every card: a guideline that matched but has
no text, a drug that falls through every bucket, a claim with no source. Results are written to the release as
`checks.json` and shown in the QC tab. A failing check never blocks a release (the data is still useful); it is shown
loudly instead. Add a check whenever a bug is found, so it can't come back silently.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import duckdb


@dataclass(frozen=True)
class Check:
    id: str
    title: str
    why: str
    tables: tuple[str, ...]
    sql: str  # returns one row per violation; first column is a short label for that row


CHECKS: list[Check] = [
    Check("claims.sourced", "Every finding names its source, evidence grade and call confidence",
          "A claim without a source or grade can't be judged or refreshed.", ("claims",),
          """SELECT claim_id FROM claims WHERE sources IS NULL OR sources IN ('', '[]') OR evidence_level IS NULL
             OR call_level IS NULL OR overall IS NULL"""),
    Check("claims.unique", "Every finding has a unique id", "Ids are how 'What changed' matches findings across "
          "releases; duplicates would hide changes.", ("claims",),
          "SELECT claim_id FROM claims GROUP BY claim_id HAVING count(*) > 1"),
    Check("claims.drug_response_not_sensitive", "Drug-response findings are never hidden as 'sensitive'",
          "Sensitive is for disease risk you may not want to see; medicine responses should always be visible.",
          ("claims",), "SELECT claim_id FROM claims WHERE \"group\" = 'drug_response' AND sensitive"),
    Check("claims.common_is_limited", "Disease findings for common variants are graded Limited",
          "A variant carried by more than 5% of some population can't on its own cause a rare disease (ACMG BA1).",
          ("claims",),
          """SELECT coalesce(gene, rsid, claim_id) FROM claims
             WHERE category IN ('pathogenic', 'likely_pathogenic', 'conflicting', 'predicted')
               AND greatest(coalesce(af_1kg, 0), coalesce(gnomad_popmax_af, 0)) >= 0.05
               AND evidence_level <> 'Limited'"""),
    Check("pgx.matched_has_text", "Every guideline that matches you has recommendation text",
          "A matched guideline with no text shows as an empty card.", ("pgx_drugs",),
          """SELECT drug || ' (' || source || ')' FROM pgx_drugs
             WHERE matched AND source IN ('CPIC', 'DPWG') AND coalesce(trim(recommendation), '') = ''"""),
    Check("pgx.tier_valid", "Every guideline row has a valid tier",
          "Buckets on the Medicines page come only from the tier set by the pipeline.", ("pgx_drugs",),
          """SELECT drug || ' (' || source || ')' FROM pgx_drugs
             WHERE tier IS NULL OR tier NOT IN ('change', 'note', 'standard', 'none')
                OR (NOT matched AND tier NOT IN ('note', 'none')) OR (matched AND tier = 'none')"""),
    Check("pgx.tier_from_flags", "'Guidance differs' is set only by the curated dose/alternative-drug flags",
          "The bucket must follow ClinPGx's structured flags, never the wording of the recommendation.",
          ("pgx_drugs",),
          """SELECT drug || ' (' || source || ')' FROM pgx_drugs WHERE matched AND source IN ('CPIC', 'DPWG') AND
             ((tier = 'change') <> (list_contains(flags, 'dosingInformation')
                                    OR list_contains(flags, 'alternateDrugAvailable')))
             AND NOT (tier = 'standard' AND (classification ILIKE 'no recommendation' OR genes = ['CFTR']))"""),
    Check("pgx.notes_are_guidance", "Drug notes without a recommendation come only from CPIC/DPWG guidelines",
          "PharmCAT also attaches data caveats to FDA label rows; those are not prescribing advice.", ("pgx_drugs",),
          """SELECT drug || ' (' || source || ')' FROM pgx_drugs
             WHERE NOT matched AND tier = 'note' AND source NOT IN ('CPIC', 'DPWG')"""),
    Check("pgx.gene_described", "Every medicine gene has a sourced description",
          "Gene descriptions come from MedlinePlus or NCBI Gene, not from text written into the app.",
          ("pgx_genes", "gene_about"),
          """SELECT gene FROM pgx_genes WHERE gene NOT IN
             (SELECT gene FROM gene_about WHERE medlineplus_text IS NOT NULL OR ncbi_summary IS NOT NULL)"""),
    Check("pgx.gene_called_has_level", "Every medicine gene states how sure the call is",
          "A diplotype without a confidence level reads as certain.", ("pgx_genes",),
          "SELECT gene FROM pgx_genes WHERE call_level IS NULL OR evidence_level IS NULL"),
]


def run(tables: dict[str, Path]) -> dict:
    con = duckdb.connect()
    for name, path in tables.items():
        con.execute(f"CREATE VIEW \"{name}\" AS SELECT * FROM read_parquet('{path}')")
    results = []
    for c in CHECKS:
        r = {"id": c.id, "title": c.title, "why": c.why}
        if not all(t in tables for t in c.tables):
            results.append({**r, "status": "skipped", "detail": f"needs {', '.join(c.tables)}"})
            continue
        try:
            bad = [str(row[0]) for row in con.execute(c.sql).fetchall()]
        except duckdb.Error as e:
            results.append({**r, "status": "error", "detail": str(e).splitlines()[0]})
            continue
        results.append({**r, "status": "fail" if bad else "pass", "count": len(bad), "examples": bad[:10]})
    con.close()
    counts = {s: sum(1 for r in results if r["status"] == s) for s in ("pass", "fail", "error", "skipped")}
    return {"counts": counts, "checks": results}
