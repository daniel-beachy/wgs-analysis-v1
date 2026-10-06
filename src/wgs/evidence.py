"""Evidence model (ADR-013): every claim carries two independent grades, and the overall grade is the weaker.

* **Evidence strength** — how good is the science linking this genotype to this outcome?
  Strong / Moderate / Limited. Derived from the curating body's review level (e.g. ClinVar stars,
  CPIC level, ClinGen gene validity) and then *downgraded* for known red flags: a "disease" variant
  that is common in healthy people (ACMG/AMP BA1/BS1), a gene whose disease link is disputed, or
  conflicting submissions. Nothing is ever upgraded beyond what the curators themselves concluded.

* **Call confidence** — how sure are we that *you* really have this genotype?
  High / Medium / Low / Not callable. From the variant caller's filter, genotype quality (GQ),
  read depth (DP), allele balance (VAF) and agreement with the provider's independent genotype file.

Every grade comes with plain-language reasons so the dashboard can show *why*, not just *what*.
The rules live in this one module so they can be read, tested and versioned (MODEL_VERSION is
recorded in each release; changing a rule shows up in "What changed").
"""

from __future__ import annotations

MODEL_VERSION = "1"

EVIDENCE = ["Limited", "Moderate", "Strong"]
CALL = ["Not callable", "Low", "Medium", "High"]
OVERALL = ["Not callable", "Limited", "Moderate", "Strong"]

# ClinVar review status → stars → base evidence strength
STARS_TO_EVIDENCE = {4: "Strong", 3: "Strong", 2: "Moderate", 1: "Limited", 0: "Limited"}
STAR_TEXT = {
    4: "ClinVar 4★: part of a professional practice guideline",
    3: "ClinVar 3★: reviewed by a ClinGen expert panel",
    2: "ClinVar 2★: two or more labs agree, with their criteria published",
    1: "ClinVar 1★: a single lab's assessment (or labs disagree)",
    0: "ClinVar 0★: submitted without assessment criteria",
}

# ClinGen gene–disease validity; disputed/refuted links undercut any variant-level "pathogenic"
GENE_VALIDITY_RANK = {"Definitive": 5, "Strong": 4, "Moderate": 3, "Limited": 2, "Disputed": 1, "Refuted": 0,
                      "No Known Disease Relationship": 0}

# Allele-frequency red flags for pathogenic claims (ACMG/AMP BA1 = 5 %, BS1 ≈ "greater than expected";
# 1 % is a conservative generic BS1 threshold for rare Mendelian disease)
BA1_AF = 0.05
BS1_AF = 0.01

DISEASE_CLASSES = {"pathogenic", "likely_pathogenic"}


def _down(level: str, steps: int = 1) -> str:
    return EVIDENCE[max(0, EVIDENCE.index(level) - steps)]


def clinvar_evidence(stars: int, sig_class: str, *, popmax_af: float | None = None, kg_af: float | None = None,
                     gene_validity: str | None = None, low_penetrance: bool = False) -> tuple[str, list[str]]:
    """Evidence strength for a ClinVar assertion, with reasons (positive reasons first)."""
    stars = int(stars or 0)
    level = STARS_TO_EVIDENCE.get(stars, "Limited")
    reasons = [STAR_TEXT.get(stars, f"ClinVar {stars}★")]
    if sig_class == "drug_response" and level == "Strong":
        level = "Moderate"
        reasons.append("ClinVar drug-response entries summarise PharmGKB annotations of mixed strength; dosing "
                       "guidance needs a CPIC/DPWG guideline (see Pharmacogenomics)")
    if sig_class == "conflicting":
        level = "Limited"
        reasons.append("Labs disagree about this variant (conflicting classifications)")
    if sig_class in DISEASE_CLASSES:
        af = max([a for a in (popmax_af, kg_af) if a is not None], default=None)
        if af is not None and af >= BA1_AF:
            level = "Limited"
            reasons.append(f"Common in healthy people (up to {af:.1%}) — too common to cause a rare disease on its "
                           "own (ACMG BA1)")
        elif af is not None and af >= BS1_AF:
            level = _down(level)
            reasons.append(f"More common than expected for a rare disease variant ({af:.1%}; ACMG BS1)")
        elif popmax_af is not None:
            reasons.append(f"Rare in gnomAD (max {popmax_af:.3%} in any population), consistent with disease")
        elif kg_af is None:
            reasons.append("Not seen in 1000 Genomes or gnomAD")
        if gene_validity:
            rank = GENE_VALIDITY_RANK.get(gene_validity)
            if rank is not None and rank <= 1:
                level = "Limited"
                reasons.append(f"ClinGen rates this gene's disease link as {gene_validity}")
            elif rank == 2:
                level = _down(level)
                reasons.append("ClinGen rates this gene's disease link as only Limited")
            elif rank is not None:
                reasons.append(f"ClinGen: gene–disease link is {gene_validity}")
        if low_penetrance:
            reasons.append("Classified as low penetrance: many carriers never develop the condition")
    return level, reasons


def call_confidence(*, filter: str | None, gq: int | None, dp: int | None, vaf: float | None,
                    zygosity: str | None, provider: str | None = None) -> tuple[str, list[str]]:
    """How sure we are of the genotype. `provider` is 'agree', 'disagree' or None (not on their file)."""
    reasons: list[str] = []
    if filter not in (None, "PASS", "."):
        return "Low", [f"Variant caller flagged this site ({filter}) — not a confident call"]
    gq, dp = gq or 0, dp or 0
    if gq >= 30 and dp >= 15:
        level = "High"
        reasons.append(f"Strong read support: {dp} reads, genotype quality {gq}")
    elif gq >= 20 and dp >= 10:
        level = "Medium"
        reasons.append(f"Moderate read support: {dp} reads, genotype quality {gq}")
    else:
        level = "Low"
        reasons.append(f"Thin read support: {dp} reads, genotype quality {gq}")
    if vaf is not None and zygosity:
        if zygosity == "het" and not 0.2 <= vaf <= 0.8:
            level = CALL[max(1, CALL.index(level) - 1)]
            reasons.append(f"Unbalanced reads for a heterozygote ({vaf:.0%} show the variant)")
        elif zygosity == "hom" and vaf < 0.8:
            level = CALL[max(1, CALL.index(level) - 1)]
            reasons.append(f"Only {vaf:.0%} of reads show the variant for a two-copy call")
    if provider == "agree":
        reasons.append("tellmeGen's own genotype file agrees")
    elif provider == "disagree":
        level = "Low"
        reasons.append("tellmeGen's genotype file disagrees with this call")
    return level, reasons


def overall(evidence: str, call: str) -> str:
    """The weaker of the two grades, on the shared Strong / Moderate / Limited / Not callable scale."""
    if call == "Not callable":
        return "Not callable"
    call_cap = {"High": "Strong", "Medium": "Moderate", "Low": "Limited"}[call]
    return OVERALL[min(OVERALL.index(evidence), OVERALL.index(call_cap))]
