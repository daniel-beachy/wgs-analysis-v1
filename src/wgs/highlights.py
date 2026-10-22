"""Highlights (Step 11): every notable way the genome differs from a typical one, ranked Strong to Limited.

Nothing here is new interpretation. Each highlight is a row that another stage already graded and explained, picked
by a fixed rule (ADR-025): single-gene findings and carrier status, computer predictions (as hypotheses), risk factors
graded at least Moderate, medicine genes whose effect isn't normal, notable HLA associations, polygenic percentiles in
the top or bottom 5% of a featured trait whose studies don't disagree, uncommon genotypes for well-supported single-
variant traits, and anything the hard-to-read-gene, repeat and deletion/duplication scans flagged. A highlight keeps
its source row's grade, statement and sources; `notable_because` names the rule that picked it.
"""

from __future__ import annotations

from pathlib import Path

INPUTS = ["claims", "pgx_genes", "pgx_drugs", "hla_associations", "pgs_traits", "trait_snps", "special_genes",
          "repeats", "cnv_events"]
EXTREME = 5          # polygenic percentile at or below this, or at or above 100 minus it, is a highlight
UNCOMMON = 0.25      # a single-variant trait genotype shared by fewer than this share of people is a highlight
NOT_NORMAL = ("poor", "intermediate", "decreased", "increased", "rapid", "ultrarapid", "positive")

COLUMNS = ["kind", "kind_label", "section", "title", "statement", "notable_because", "overall", "evidence_level",
           "call_level", "ref_table", "ref_key", "sensitive", "extra", "sources"]
KIND_ORDER = ["monogenic", "special", "cnv", "repeat", "carrier", "pgx", "hla", "risk", "pgs", "trait", "predicted"]


def _percentile(col: str) -> str:
    """'below the 1st', 'the 3rd', 'the 97th', 'above the 99th': the dashboard's wording for a percentile."""
    n = f"CAST(round({col}) AS INTEGER)"
    return (f"CASE WHEN {col} < 0.5 THEN 'below the 1st' WHEN {col} > 99.5 THEN 'above the 99th' ELSE 'the ' || {n} "
            f"|| CASE WHEN {n} % 100 BETWEEN 11 AND 13 THEN 'th' WHEN {n} % 10 = 1 THEN 'st' WHEN {n} % 10 = 2 "
            f"THEN 'nd' WHEN {n} % 10 = 3 THEN 'rd' ELSE 'th' END END || ' percentile'")


def _claims(p: Path, kind: str, label: str, section: str, where: str) -> str:
    return f"""SELECT '{kind}', '{label}', '{section}', coalesce(gene, subject), statement,
        coalesce(evidence_reasons[1], category_label), overall, evidence_level, call_level, 'claims', claim_id,
        coalesce(sensitive, false), NULL, sources
      FROM '{p}' WHERE {where}"""


def branches(t: dict[str, Path]) -> list[str]:
    """One SELECT per rule, over whichever input tables exist."""
    out = []
    if c := t.get("claims"):
        out += [
            _claims(c, "monogenic", "Single-gene finding", "health", "section = 'health' AND \"group\" = 'monogenic'"),
            _claims(c, "carrier", "Carrier status", "carrier", "section = 'carrier'"),
            _claims(c, "risk", "Risk factor", "health",
                    "section = 'health' AND \"group\" = 'risk' AND overall IN ('Strong', 'Moderate')"),
            _claims(c, "predicted", "Computer prediction (a hypothesis)", "health",
                    "section = 'health' AND \"group\" = 'predicted'"),
        ]
    if g := t.get("pgx_genes"):
        drugs = t.get("pgx_drugs")
        # as on the Medicines page: only a CPIC or DPWG guideline that matched you puts a medicine in "change"
        changed = (f"(SELECT count(DISTINCT drug) FROM '{drugs}' d WHERE d.tier = 'change' AND d.matched "
                   f"AND d.source IN ('CPIC', 'DPWG') AND list_contains(d.genes, g.gene))") if drugs else "0"
        def claim(col: str) -> str:  # the gene's pgx claim carries the statement and the ClinPGx link
            return (f"(SELECT any_value(c.{col}) FROM '{t['claims']}' c WHERE c.\"group\" = 'pgx_gene' "
                    f"AND c.gene = g.gene)") if t.get("claims") else "NULL"
        out.append(f"""SELECT 'pgx', 'Medicine gene', 'pgx', g.gene,
            coalesce({claim('statement')}, g.phenotype),
            g.phenotype || CASE WHEN {changed} > 0 THEN '; the CPIC/DPWG advice differs for ' || {changed}
              || CASE WHEN {changed} = 1 THEN ' medicine' ELSE ' medicines' END ELSE '' END,
            g.overall, g.evidence_level, g.call_level, 'pgx_genes', g.gene, false, NULL, {claim('sources')}
          FROM '{g}' g WHERE g.phenotype_class IN {NOT_NORMAL}""")
    if h := t.get("hla_associations"):
        out.append(f"""SELECT 'hla', 'Immune-system (HLA) type', 'health', condition, statement,
            CASE kind WHEN 'rule_out' THEN 'Rule-out test: you lack ' || markers
              ELSE 'You carry the HLA type behind a ' || direction || ' association' END,
            overall, evidence_level, call_level,
            'hla_associations', association_id, false, NULL, quotes
          FROM '{h}' WHERE notable""")
    if p := t.get("pgs_traits"):
        out.append(f"""SELECT 'pgs', 'Polygenic score', section, label, agreement_text,
            'In the ' || CASE WHEN percentile >= 50 THEN 'top' ELSE 'bottom' END || ' {EXTREME}%: '
              || {_percentile('percentile')} || ' (median across studies)',
            overall, evidence_level, call_confidence, 'pgs_traits', trait_id, false, abs_risk,
            '[{{"source": "pgs_catalog", "record": "' || coalesce(risk_pgs_id, best_pgs_id) || '", "url": '
              || '"https://www.pgscatalog.org/score/' || coalesce(risk_pgs_id, best_pgs_id) || '/"}}]'
          FROM '{p}' WHERE coalesce(section, '') <> '' AND (percentile <= {EXTREME} OR percentile >= {100 - EXTREME})
            AND coalesce(agreement, '') <> 'disagree'""")
    if s := t.get("trait_snps"):
        out.append(f"""SELECT 'trait', 'Trait', 'traits', label, effect_text,
            'About ' || CAST(round(genotype_share * 100) AS INTEGER) || ' in 100 people share your genotype',
            overall, evidence_level, call_confidence, 'trait_snps', rsid, false, NULL,
            '[{{"source": "gwas_catalog", "record": "PMID ' || pmid || '", "url": '
              || '"https://pubmed.ncbi.nlm.nih.gov/' || pmid || '/"}}]'
          FROM '{s}' WHERE genotype_share < {UNCOMMON} AND overall IN ('Strong', 'Moderate')""")
    if s := t.get("special_genes"):
        out.append(f"""SELECT 'special', 'Hard-to-read gene', 'carrier', gene, statement,
            'Re-read with ' || coalesce(tool, 'a specialised caller'),
            overall, evidence_level, call_level, 'special_genes', gene, false, NULL, sources
          FROM '{s}' WHERE status NOT IN ('typical', 'not_determined')""")
    if r := t.get("repeats"):
        out.append(f"""SELECT 'repeat', 'Repeat expansion', 'health', coalesce(gene, locus_id), disease,
            'Repeat length in the ' || replace(status, '_', ' '), overall, evidence_level, call_level, 'repeats',
            locus_id, coalesce(sensitive, false), NULL,
            '[{{"source": "strchive", "record": "' || locus_id || '", "url": "' || coalesce(url, '') || '"}}]'
          FROM '{r}' WHERE status IN ('intermediate', 'disease_range')""")
    if v := t.get("cnv_events"):
        out.append(f"""SELECT 'cnv', 'Deletion or duplication', 'health', array_to_string(genes, ', '),
            replace(status, '_', ' '), 'Seen by both breakpoint reads and depth; rare in gnomAD-SV',
            overall, evidence_level, call_level, 'cnv_events', event_id, true, NULL, NULL
          FROM '{v}' WHERE status LIKE 'possible_%'""")
    return out


def sql(t: dict[str, Path]) -> str | None:
    """The highlights table, ranked by overall grade, then kind, then title. None when no input exists."""
    parts = branches(t)
    if not parts:
        return None
    cols = ", ".join(COLUMNS)
    kinds = ", ".join(f"'{k}'" for k in KIND_ORDER)
    union = "\nUNION ALL\n".join(f"SELECT * FROM ({b}) AS b({cols})" for b in parts)
    return f"""SELECT row_number() OVER (ORDER BY CASE overall WHEN 'Strong' THEN 0 WHEN 'Moderate' THEN 1
                 WHEN 'Limited' THEN 2 ELSE 3 END, list_position([{kinds}], kind), title) AS rank,
             md5(ref_table || ':' || ref_key) AS highlight_id, *
      FROM ({union}) ORDER BY rank"""
