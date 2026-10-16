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
    Check("claims.exact_record", "Every finding links to its exact source record, not just a website",
          "A link to a database's home page can't be checked; each claim must point at the record it came from.",
          ("claims",),
          """SELECT claim_id FROM claims WHERE NOT list_bool_or(list_transform(
               coalesce(json_extract_string(sources, '$[*].url'), []),
               u -> coalesce(regexp_matches(u, '://[^/]+/[^?#]*[^/?#]/[^/?#]+|[?][^#]*='
                                                   || '|pubmed[.]ncbi[.]nlm[.]nih[.]gov/[0-9]+'), false)))"""),
    Check("claims.ba1_role", "Common variants never get a carrier or affected label",
          "An allele carried by 5% or more of a population can't cause a rare disease on its own (ACMG/AMP BA1); "
          "only expert-panel (3★+) ClinVar records, which already weighed frequency, are exempt.", ("claims",),
          """SELECT coalesce(gene, rsid, claim_id) FROM claims
             WHERE role IN ('carrier', 'affected', 'possible') AND coalesce(clinvar_stars, 0) < 3
               AND greatest(coalesce(af_1kg, 0), coalesce(gnomad_popmax_af, 0)) >= 0.05"""),
    Check("claims.no_rarer_allele_text", "No finding blames 'the rarer allele' when yours is the common one",
          "This wording was wrong for variants where the reference genome carries the rare version.", ("claims",),
          """SELECT claim_id FROM claims WHERE statement ILIKE '%rarer allele%' OR role_reason ILIKE '%rarer allele%'
             OR CAST(evidence_reasons AS VARCHAR) ILIKE '%rarer allele%'"""),
    Check("claims.protein_hgvs", "Protein changes use standard HGVS names (p.Arg402Gln)",
          "Readers and search engines recognise the standard form; anything else is a formatting bug.", ("claims",),
          """SELECT coalesce(gene, claim_id) || ' ' || aa_change FROM claims WHERE aa_change IS NOT NULL
             AND (NOT regexp_matches(aa_change, '^p[.]') OR length(aa_change) > 40)"""),
    Check("claims.brain_agree", "Findings and the Brain & mind list agree on your genotype",
          "The same variant must not read as one copy on one page and two on another.",
          ("claims", "brain_variants"),
          """SELECT b.gene || ' ' || b.chrom || ':' || b.pos FROM brain_variants b JOIN claims c
               ON c.chrom = b.chrom AND c.pos = b.pos AND c.alt = b.alt
             WHERE c.zygosity IS NOT NULL AND b.zygosity IS NOT NULL AND c.zygosity <> b.zygosity"""),
    Check("brain.sourced", "Every Brain & mind variant names its sources", "Each row needs links to the records "
          "behind it (SFARI Gene, gnomAD, ClinVar when classified).", ("brain_variants",),
          """SELECT gene || ' ' || chrom || ':' || pos FROM brain_variants
             WHERE sources IS NULL OR sources IN ('', '[]')"""),
    Check("traits.sourced", "Every trait variant names its GWAS study and paper",
          "A trait association must point at the study (GCST) and publication it came from.", ("trait_snps",),
          "SELECT rsid FROM trait_snps WHERE study IS NULL OR pmid IS NULL"),
    Check("claims.penetrance_single_gene", "Every 'how often it leads to disease' figure was rated for that gene alone",
          "ClinGen sometimes rates several genes together; that likelihood is not a figure for any one of them.",
          ("claims",),
          """SELECT gene FROM claims WHERE json_type(json_extract(actionability, '$.penetrance')) = 'OBJECT'
             AND coalesce(json_array_length(json_extract(actionability, '$.penetrance.genes')), 0) <> 1"""),
    Check("pgs.abs_risk_sourced", "Every 'your risk in numbers' result cites both of its numbers",
          "An absolute risk is a typical risk times the score's effect; both must link to their source.",
          ("pgs_scores",),
          """SELECT pgs_id FROM pgs_scores WHERE json_extract_string(abs_risk, '$.status') = 'ok'
             AND (json_extract_string(abs_risk, '$.baseline.url') IS NULL
                  OR json_extract_string(abs_risk, '$.effect.url') IS NULL
                  OR json_extract_string(abs_risk, '$.effect.ppm_id') IS NULL
                  OR (json_extract_string(abs_risk, '$.baseline.kind') <> 'cohort'
                      AND json_extract_string(abs_risk, '$.baseline.quote') IS NULL)
                  OR CAST(json_extract(abs_risk, '$.you') AS DOUBLE) NOT BETWEEN 0 AND 1)"""),
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
    Check("pgs.featured_scored", "Every featured trait has a result from at least one score",
          "A featured trait shown without a result would be a silent gap.", ("pgs_traits",),
          "SELECT trait_id || ' ' || label FROM pgs_traits WHERE featured AND percentile IS NULL"),
    Check("pgs.trait_headline", "Each trait's headline is the median across studies of its scores",
          "Recomputed from the score rows, so the headline can't drift from the scores shown beneath it, and several "
          "scores from one paper count as one study (ADR-019).", ("pgs_scores", "pgs_traits"),
          """WITH ok AS (SELECT * FROM pgs_scores WHERE passed AND percentile IS NOT NULL AND z IS NOT NULL),
               s AS (SELECT trait_group, coalesce(pgp_id, pgs_id) AS study, median(percentile) AS p
                     FROM ok GROUP BY ALL),
               t AS (SELECT trait_group, median(p) AS p, count(*) AS n FROM s GROUP BY ALL),
               r AS (SELECT trait_group, min(percentile) AS lo, max(percentile) AS hi, count(*) AS n
                     FROM ok GROUP BY ALL)
             SELECT x.trait_id FROM pgs_traits x LEFT JOIN t ON t.trait_group = x.trait_id
               LEFT JOIN r ON r.trait_group = x.trait_id
             WHERE abs(coalesce(x.percentile, -1) - coalesce(t.p, -1)) > 1e-9 OR x.n_studies <> coalesce(t.n, 0)
                OR x.n_scored <> coalesce(r.n, 0) OR abs(coalesce(x.low, -1) - coalesce(r.lo, -1)) > 1e-9
                OR abs(coalesce(x.high, -1) - coalesce(r.hi, -1)) > 1e-9
                OR x.percentile NOT BETWEEN x.low AND x.high"""),
    Check("pgs.trait_risk_named",
          "A trait's risk in numbers comes from the most-tested score that has one, and names it",
          "The absolute risk uses one score's published effect, so it must say which; picking any other score "
          "would let a weakly tested outlier set the number (ADR-019).", ("pgs_scores", "pgs_traits"),
          """WITH ok AS (SELECT *, row_number() OVER (PARTITION BY trait_group ORDER BY
                   (json_extract_string(abs_risk, '$.status') = 'ok') DESC, independent_pubs DESC, independent_n DESC,
                   evaluated_n DESC) AS k
                 FROM pgs_scores WHERE passed AND percentile IS NOT NULL AND z IS NOT NULL AND abs_risk IS NOT NULL)
             SELECT x.trait_id FROM pgs_traits x LEFT JOIN ok ON ok.trait_group = x.trait_id AND ok.k = 1
             WHERE coalesce(x.risk_pgs_id, '') <> coalesce(ok.pgs_id, '')
                OR coalesce(x.abs_risk, '') <> coalesce(ok.abs_risk, '')"""),
    Check("pgs.sourced", "Every computed polygenic score links to its paper",
          "Each percentile on the page must lead to the study that built the score.", ("pgs_scores",),
          "SELECT pgs_id FROM pgs_scores WHERE percentile IS NOT NULL AND pmid IS NULL AND doi IS NULL"),
    Check("ancestry.shares_valid", "Continental shares are between 0 and 1 and add up to 100%",
          "Shares that do not add up mean the fit or its inputs went wrong.", ("ancestry_components",),
          """SELECT 'sum ' || round(sum(share), 4) FROM ancestry_components HAVING count(*) > 0
             AND (abs(sum(share) - 1) > 1e-6 OR min(share) < 0 OR max(share) > 1)"""),
    Check("ancestry.noise_labelled", "Shares below their noise floor are never presented as found",
          "Reference people with none of an ancestry still get a small share by chance; a share below that floor "
          "must be graded Limited and marked not detected (ADR-020).", ("ancestry_components",),
          """SELECT superpop FROM ancestry_components WHERE (share <= noise_floor AND (detected OR evidence_level <>
             'Limited')) OR (share > noise_floor AND NOT detected)"""),
    Check("ancestry.neighbours_counted", "Every nearest reference person is counted once",
          "The nearest-population vote must account for exactly the 10 neighbours used (ancestry.K_NEIGHBOURS).",
          ("ancestry_populations",),
          "SELECT 'sum ' || sum(neighbours) FROM ancestry_populations HAVING count(*) > 0 AND sum(neighbours) <> 10"),
    Check("haplogroups.sourced", "Every haplogroup names its tool, tree, quality and a link to the tree",
          "A lineage label is only checkable with the tree it was placed on.", ("haplogroups",),
          """SELECT line FROM haplogroups WHERE tool IS NULL OR tree IS NULL OR url IS NULL OR quality IS NULL
             OR quality NOT BETWEEN 0 AND 1 OR call_level IS NULL"""),
    Check("haplogroups.no_placeholders", "Haplogroup mutation lists contain only mutations",
          "An empty haplogrep field was once read as the text 'nan' and counted as a missing mutation.",
          ("haplogroups",),
          """SELECT line FROM haplogroups WHERE regexp_matches(coalesce(found, '') || coalesce(missing, '')
             || coalesce(extra, ''), '"(nan|none|null)"', 'i')"""),
    Check("special.genes_sourced", "Every hard-to-read gene result names its tool, grades and sources",
          "SMN1 and GBA1 are re-read by specialised callers; each result must say which and how sure.",
          ("special_genes",),
          """SELECT gene FROM special_genes WHERE tool IS NULL OR sources IS NULL OR sources IN ('', '[]')
             OR call_level IS NULL OR overall IS NULL"""),
    Check("special.undetermined_not_reassuring", "An undetermined result is never shown as typical",
          "When a caller can't decide (e.g. SMN1 copy number between 1 and 2), the page must say 'not determined' "
          "and grade it Not callable, never fall back to 'not a carrier'.", ("special_genes",),
          """SELECT gene FROM special_genes WHERE (status = 'not_determined') <> (call_level = 'Not callable')"""),
    Check("repeats.unreadable_not_typical", "Repeats that couldn't be read are never called typical",
          "No reads across a repeat (e.g. FMR1, which is very GC-rich) means unknown, not normal.", ("repeats",),
          """SELECT locus_id FROM repeats WHERE (status = 'not_readable') <> (call_level = 'Not callable')
             OR (status = 'typical' AND (alleles IS NULL OR alleles = ''))"""),
    Check("repeats.out_of_range_is_a_claim", "Every repeat outside the typical range is also a graded finding",
          "A repeat in the intermediate or disease range must appear on the Health or Carrier page, not only in "
          "the table.", ("repeats", "claims"),
          """SELECT r.locus_id FROM repeats r WHERE r.status IN ('intermediate', 'disease_range')
             AND NOT EXISTS (SELECT 1 FROM claims c WHERE c.claim_id = 'special:repeat:' || r.locus_id)"""),
    Check("hla.quoted", "Every HLA association quotes its source", "The HLA statements are only as good as the "
          "quoted study; each must carry at least one quote with a PubMed link (ADR-021).", ("hla_associations",),
          """SELECT association_id FROM hla_associations WHERE coalesce(json_array_length(quotes), 0) = 0
             OR NOT list_bool_and(list_transform(CAST(json_extract(quotes, '$[*].url') AS VARCHAR[]),
                                                 u -> regexp_matches(u, 'pubmed[.]ncbi[.]nlm[.]nih[.]gov/[0-9]+')))"""),
    Check("blood_group.sourced", "Every blood-group result quotes its source and was read, not assumed",
          "The O position matches the reference genome, so 'no variant' means O; that is only safe where the "
          "position was readable, and the result must cite the studies behind the rule.", ("blood_group",),
          """SELECT system FROM blood_group WHERE coalesce(json_array_length(quotes), 0) = 0
             OR (result IS NULL) <> (call_level = 'Not callable')"""),
    Check("cnv.two_signals", "A deletion or duplication is only reported as possible when two methods agree",
          "Read depth alone dips at GC-rich, GC-poor and hard-to-map exons; a reported copy-number change needs "
          "breakpoint reads (Delly) and a GC-corrected depth drop, and must not be common in gnomAD-SV (ADR-021).",
          ("cnv_events",),
          """SELECT event_id FROM cnv_events WHERE status IN ('possible_exon_deletion', 'possible_exon_duplication')
             AND NOT (delly AND depth_support)"""),
    Check("cnv.findings_are_claims", "Every possible deletion or duplication is also a graded finding",
          "A copy-number change in an ACMG gene must appear on the Health page, not only in the scan table.",
          ("cnv_events", "claims"),
          """SELECT e.event_id FROM cnv_events e WHERE e.status IN ('possible_exon_deletion',
             'possible_exon_duplication') AND NOT EXISTS (SELECT 1 FROM claims c
             WHERE c.claim_id = 'special:cnv:' || e.event_id)"""),
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
