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
    Check("traits.lift_sane", "Trait variants sit at their GRCh37 positions",
          "A wrong genome-build lift would read a different base; rs12913832 (eye colour) is a fixed landmark.",
          ("trait_snps",),
          """SELECT rsid || ' at ' || chrom || ':' || pos FROM trait_snps
             WHERE rsid = 'rs12913832' AND NOT (chrom = '15' AND pos = 28365618)
             UNION ALL SELECT 'rs12913832 missing' WHERE NOT EXISTS
                 (SELECT 1 FROM trait_snps WHERE rsid = 'rs12913832')"""),
    Check("traits.graded", "Every trait variant has an evidence grade and call confidence",
          "Ungraded rows read as certain.", ("trait_snps",),
          "SELECT rsid FROM trait_snps WHERE evidence_level IS NULL OR call_confidence IS NULL OR overall IS NULL"),
    Check("traits.direction_unknown_limited", "A variant with no known effect allele is never graded above Limited",
          "Without the effect allele your genotype cannot be read as higher or lower.", ("trait_snps",),
          "SELECT rsid FROM trait_snps WHERE NOT direction_known AND evidence_level <> 'Limited'"),
    Check("pgs.match_rate", "Every scored polygenic score matched at least 75% of its variants",
          "Below this, the score is mostly filled-in averages and its percentile is meaningless (pgsc_calc rule).",
          ("pgs_scores",),
          """SELECT pgs_id || ' (' || coalesce(round(100 * match_rate, 1)::VARCHAR, '?') || '%)' FROM pgs_scores
             WHERE call_confidence <> 'Not callable' AND (match_rate IS NULL OR match_rate < 0.75)"""),
    Check("pgs.percentile_range", "Every percentile lies between 0 and 100",
          "A percentile outside 0–100 means the adjustment step misread its input.", ("pgs_scores",),
          "SELECT pgs_id FROM pgs_scores WHERE percentile IS NOT NULL AND (percentile < 0 OR percentile > 100)"),
    Check("pgs.featured_scored", "Every featured polygenic score produced a percentile",
          "A featured trait shown without a result would be a silent gap.", ("pgs_scores",),
          "SELECT pgs_id || ' ' || label FROM pgs_scores WHERE featured AND percentile IS NULL"),
    Check("pgs.genotype_calibration", "You are rarely called 'reference' where almost everyone carries the other "
          "allele", "At sites where >99% of the panel carries the alternate allele, a reference call should be rare "
          "(<2%); more means calls are being made wrongly — an indel-spelling mismatch once did this to 18% of "
          "such indels, biasing every score.", ("pgs_genotype_qc",),
          """SELECT kind || ': ' || round(hom_ref_pct, 1) || '% called reference' FROM pgs_genotype_qc
             WHERE band = 6 AND sites >= 100 AND hom_ref_pct > 2"""),
    Check("pgs.score_spread", "Across all scores, your z-scores look like one person's",
          "One genome's z-scores over many unrelated scores should average near 0 with spread near 1; a shifted "
          "or squeezed set points to a systematic genotyping or scaling error rather than biology.", ("pgs_scores",),
          """SELECT 'mean z ' || round(avg(z), 2) || ', SD ' || round(stddev(z), 2) FROM pgs_scores
             WHERE z IS NOT NULL HAVING count(*) >= 20 AND (abs(avg(z)) > 0.5 OR stddev(z) NOT BETWEEN 0.6 AND 1.5)"""),
    Check("pgs.never_strong", "Polygenic scores are never graded Strong",
          "A score shifts odds modestly; it is capped below Strong by the evidence model.", ("pgs_scores",),
          "SELECT pgs_id FROM pgs_scores WHERE evidence_level = 'Strong' OR overall = 'Strong'"),
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
