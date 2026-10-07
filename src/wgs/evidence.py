"""Evidence model (ADR-013): every claim carries two independent grades, and the overall grade is the weaker.

* **Evidence strength** — how good is the science linking this genotype to this outcome?
  Strong / Moderate / Limited. Derived from the curating body's review level (e.g. ClinVar stars,
  CPIC level, ClinGen gene validity) and then *downgraded* for known red flags: a "disease" variant
  that is common in healthy people (ACMG/AMP BA1/BS1), a gene whose disease link is disputed, or
  conflicting submissions. Nothing is ever upgraded beyond what the curators themselves concluded.

* **Call confidence** — how sure are we that *you* really have this genotype?
  High / Medium / Low / Not callable. From the variant caller's filter, genotype quality (GQ),
  read depth (DP), allele balance (VAF) and agreement with the provider's independent genotype file.

Model 2 adds:

* **Inheritance-aware interpretation** — what a variant means depends on how the condition is
  inherited (see ``wgs.inheritance``). One copy in a recessive condition → carrier; one copy in a
  dominant condition → may matter for you. Inheritance is looked up per *condition* (ClinVar's
  condition IDs → HPO / Orphanet via Mondo) before falling back to everything known about the gene.
* **Medically actionable genes** — ACMG SF v3.3 reporting rules and ClinGen actionability scores.
* **Predicted evidence** — AlphaMissense, REVEL and gene constraint for rare variants no expert has
  classified. Always labelled "predicted" and capped at Limited.

Every grade comes with plain-language reasons so the dashboard can show *why*, not just *what*.
The rules live in this one module so they can be read, tested and versioned (MODEL_VERSION is
recorded in each release; changing a rule shows up in "What changed").
"""

from __future__ import annotations

import re

MODEL_VERSION = "4"

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


def too_common(*, popmax_af: float | None, kg_af: float | None, stars: int | None) -> str | None:
    """ACMG/AMP BA1: an allele at ≥5% in a population is too common to cause a rare Mendelian disease on its own.

    Returns the reason when BA1 applies, so carrier/affected labels are withheld. Expert-panel (3★+) ClinVar
    classifications already weighed frequency (ClinGen's BA1 exceptions such as HFE C282Y), so they are kept."""
    if int(stars or 0) >= 3:
        return None
    af = max([a for a in (popmax_af, kg_af) if a is not None], default=None)
    if af is None or af < BA1_AF:
        return None
    share = "over 99%" if af >= 0.995 else f"{af:.0%}"
    return (f"Carried on up to {share} of chromosomes in some populations: too common to cause a rare inherited "
            "condition on its own (ACMG/AMP rule BA1), so carrier/affected labels are not applied")


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


# ---------------------------------------------------------------------------------------------------------------
# Model 2: inheritance, actionability, predictions
# ---------------------------------------------------------------------------------------------------------------

LOF = {"stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost"}
REVEL_SUPPORTING, REVEL_MODERATE, REVEL_STRONG = 0.644, 0.773, 0.932   # Pejaver et al. 2022 (ClinGen PP3)
AM_LIKELY_PATHOGENIC = 0.564                                           # Cheng et al. 2023
LOEUF_CONSTRAINED = 0.6                                                # gnomAD v4 guidance
PREDICTED_MAX_AF = 0.001                                               # rare enough to be worth a look
PREDICTED_CAP = "Limited"


def is_lof(consequence: str | None) -> bool:
    return bool(set((consequence or "").split("&")) & LOF)


def role(modes: set[str], zygosity: str | None, chrom: str | None) -> tuple[str, str]:
    """What carrying this genotype means given the inheritance.

    Returns ('carrier' | 'affected' | 'possible' | 'unknown', why)."""
    from . import inheritance as inh

    core = modes - {"DG", "MF", "SP"}
    if not core:
        return "unknown", "How this condition is inherited is not recorded"
    if zygosity == "hom" and core & (inh.RECESSIVE | inh.X_LINKED):
        return "affected", "Two copies of a variant in a recessive gene: both copies of the gene are affected"
    if chrom == "X" and zygosity == "hemi" and core & inh.X_LINKED:
        return "affected", "X-linked gene and a single X chromosome: the one copy you have carries the variant"
    if zygosity == "het" and core <= inh.RECESSIVE:
        return "carrier", ("Recessive condition: one copy usually makes you a healthy carrier who could pass it on; "
                           "the condition needs two non-working copies")
    if zygosity == "het" and chrom == "X" and core <= {"XLR", "XL"}:
        return "carrier", "X-linked recessive: with two X chromosomes, one copy usually makes you a carrier"
    if core & inh.DOMINANT:
        if core & inh.RECESSIVE:
            return "possible", ("This gene causes both dominant and recessive conditions: one copy may matter for the "
                                "dominant one and makes you a carrier for the recessive one")
        return "possible", "Dominant condition: one copy can be enough, though many carriers never develop it"
    return "unknown", "Inheritance pattern does not settle what one copy means"


def acmg_reportable(rule: str | None, *, copies: int, zygosity: str | None, consequence: str | None,
                    aa_change: str | None, compound: bool = False) -> bool:
    """Would a clinical lab report this P/LP genotype as an ACMG secondary finding?"""
    if not rule:
        return False
    r = rule.lower()
    if "c282y" in r:
        return zygosity == "hom" and bool(re.search(r"C282Y|282C>282Y|Cys282Tyr", aa_change or ""))
    if "biallelic" in r:
        return copies == 2 or compound
    if "truncating" in r:
        return is_lof(consequence)
    return True


def predicted_evidence(*, consequence: str | None, am_score: float | None, revel: float | None,
                       loeuf: float | None, popmax_af: float | None, kg_af: float | None,
                       recessive_gene: bool, clinvar: str | None = None) -> tuple[str | None, list[str]]:
    """Grade a computational prediction for a rare, unclassified variant. Returns (None, []) if not flagged.

    Only flags variants where predictors agree or are calibrated to strong evidence, so the list stays short
    enough to be worth reading. Never above Limited: predictions are not clinical classifications."""
    afs = [a for a in (popmax_af, kg_af) if a is not None]
    if afs and max(afs) >= PREDICTED_MAX_AF:
        return None, []
    reasons: list[str] = []
    flagged = False
    if is_lof(consequence):
        if loeuf is not None and loeuf < LOEUF_CONSTRAINED:
            flagged = True
            reasons.append(f"Predicted to break the gene ({consequence.replace('_', ' ')}), and the gene is intolerant "
                           f"of broken copies (gnomAD LOEUF {loeuf:.2f} < {LOEUF_CONSTRAINED})")
        elif recessive_gene:
            flagged = True
            reasons.append(f"Predicted to break the gene ({consequence.replace('_', ' ')}) in a gene where broken "
                           "copies cause recessive disease")
        if flagged:
            reasons.append("Not every predicted break is real (e.g. near the end of a gene the protein may survive)")
    elif "missense" in (consequence or ""):
        am_lp = am_score is not None and am_score > AM_LIKELY_PATHOGENIC
        if revel is not None and revel >= REVEL_STRONG:
            flagged = True
            reasons.append(f"REVEL {revel:.3f} ≥ {REVEL_STRONG}: ClinGen-calibrated strong computational evidence")
        elif revel is not None and revel >= REVEL_SUPPORTING and am_lp:
            flagged = True
            level = "moderate" if revel >= REVEL_MODERATE else "supporting"
            reasons.append(f"REVEL {revel:.3f} ({level} computational evidence) and AlphaMissense {am_score:.2f} "
                           "(likely pathogenic) agree")
        elif revel is None and am_score is not None and am_score >= 0.9:
            flagged = True
            reasons.append(f"AlphaMissense {am_score:.2f} (strongly likely pathogenic); no REVEL score")
    if not flagged:
        return None, []
    reasons.append("Rare: " + (f"at most {max(afs):.3%} in any population" if afs else "not seen in gnomAD or 1000 "
                                                                                        "Genomes"))
    if clinvar:
        reasons.append(f"ClinVar: {clinvar}, i.e. labs have not settled whether it causes disease; the "
                       "computer prediction is capped at Limited evidence")
    else:
        reasons.append("Computer prediction only: no lab or expert has classified this variant, so it is capped "
                       "at Limited evidence")
    return PREDICTED_CAP, reasons


# ---------------------------------------------------------------------------------------------------------------
# Model 3: common-variant associations and polygenic scores (ADR-017)
# ---------------------------------------------------------------------------------------------------------------

GENOME_WIDE_MLOG10P = 7.3          # p < 5e-8, the conventional genome-wide significance threshold
PGS_MIN_MATCH = 0.75               # pgsc_calc default: below this the score is not computed
PGS_HIGH_MATCH = 0.9
PGS_MODERATE_IND_PUBS, PGS_MODERATE_IND_N = 2, 20_000


def gwas_evidence(*, publications: int | None, mlog10p: float | None,
                  direction_known: bool) -> tuple[str, list[str]]:
    """How solid is a single-variant GWAS association? Replication across publications is what counts."""
    pubs = int(publications or 0)
    reasons = [f"Reported by {pubs} publication{'s' if pubs != 1 else ''} in the GWAS Catalog"]
    level = "Strong" if pubs >= 3 else "Moderate" if pubs == 2 else "Limited"
    if mlog10p is not None and mlog10p < GENOME_WIDE_MLOG10P:
        level = "Limited"
        reasons.append("Strongest result is below genome-wide significance (p ≥ 5×10⁻⁸)")
    elif mlog10p is not None:
        reasons.append(f"Strongest result p ≈ 10^-{mlog10p:.0f} (genome-wide significant)")
    if not direction_known:
        level = "Limited"
        reasons.append("The catalog does not record which allele has the effect, so your genotype cannot be read "
                       "as higher or lower")
    reasons.append("A common variant nudges a trait; it does not decide it")
    return level, reasons


def pgs_evidence(*, independent_pubs: int | None, independent_n: int | None,
                 evaluated_n: int | None) -> tuple[str, list[str]]:
    """How well tested is a polygenic score? Never Strong: a score shifts odds, it never diagnoses."""
    ip, inn = int(independent_pubs or 0), int(independent_n or 0)
    if ip >= PGS_MODERATE_IND_PUBS and inn >= PGS_MODERATE_IND_N:
        level = "Moderate"
        reasons = [f"Independently tested by {ip} other research groups in {inn:,} people"]
    elif ip >= 1:
        level = "Limited"
        reasons = [f"Independently tested by {ip} other research group{'s' if ip != 1 else ''} in {inn:,} people"]
    else:
        level = "Limited"
        reasons = ["Only tested by the team that built it" + (f" ({int(evaluated_n):,} people)" if evaluated_n
                                                               else "")]
    reasons.append("Polygenic scores shift the odds by modest amounts and are capped below Strong")
    return level, reasons


def pgs_call(*, match_rate: float | None, passed: bool | None) -> tuple[str, list[str]]:
    """How completely we could compute the score from your genome."""
    if not passed or match_rate is None or match_rate < PGS_MIN_MATCH:
        pct = f"{match_rate:.0%}" if match_rate is not None else "too few"
        return "Not callable", [f"Only {pct} of the score's variants could be matched (need ≥ {PGS_MIN_MATCH:.0%})"]
    level = "High" if match_rate >= PGS_HIGH_MATCH else "Medium"
    return level, [f"{match_rate:.1%} of the score's variants matched your genome and the reference panel"]
