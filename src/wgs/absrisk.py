"""Absolute risk: "about N in 100 people like you" next to "N in 100 people overall" (ADR-018).

A percentile alone hides how much a score matters. This module turns a polygenic-score position into an
absolute risk, but only when every input is published and fit for the purpose; otherwise it returns a
"not established" reason. Nothing here is hand-tuned per disease.

Inputs and their sources
* z: your score in standard deviations of the ancestry-matched reference panel (pgsc_calc, 1000 Genomes).
* Discrimination: an AUROC from a PGS Catalog evaluation of this exact score, adjusted for nothing beyond
  genetic ancestry/technical terms (AUROC with age, sex or clinical factors mixed in overstates the score), OR
  an odds/hazard ratio that a cited paper states is per 1 SD (vendored table `knowledge/data/baseline_risk.toml`).
* Baseline K: a cited population lifetime risk (same table), else the share of participants diagnosed in that
  same evaluation cohort, only when the cohort is not case-enriched.

Models (one per input type, all calibrated so the population average equals K)
* AUROC → equal-variance binormal model: cases ~ N(d, 1), non-cases ~ N(0, 1), AUROC = Φ(d/√2). This is the
  model behind Wray et al. 2010 (PLoS Genet, doi:10.1371/journal.pgen.1000864) and Pain et al. 2022 (Eur J
  Hum Genet, doi:10.1038/s41431-021-01028-z); risk(z) follows from Bayes' rule.
* OR per SD → logistic model, logit risk = a + ln(OR)·z, with a solved so the population mean risk is K.
* HR per SD → proportional hazards: risk = 1 − S^exp(ln(HR)·z), with S solved so the mean risk is K.
"""

from __future__ import annotations

import math
import re
import tomllib
from functools import cache
from pathlib import Path
from statistics import NormalDist

N = NormalDist()
TABLE = Path(__file__).parent / "knowledge" / "data" / "baseline_risk.toml"

# Covariates that only control for ancestry or technical batch; anything else (age, sex, BMI, family history,
# clinical risk factors) inflates an AUROC beyond what the score alone achieves.
TECHNICAL = {"pc", "pcs", "principal", "component", "components", "ancestry", "genetic", "genotyping", "genotype",
             "array", "arrays", "chip", "platform", "batch", "study", "studies", "centre", "center", "site",
             "assessment",
             "cohort", "of", "and", "the", "first", "top", "to", "for", "pca", "ancestral", "population", "structure",
             "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "twenty"}
MAX_CASE_SHARE = 0.20   # above this a cohort was enriched for cases, so its case share is not a population risk
MIN_CASES = 100
MIN_N = 1000
LIFETIME_PREFIX = re.compile(r"^lifetime,?\s*", re.I)
MEASURE = {"lifetime_risk": "lifetime risk", "prevalence": "share of people who have or had it (prevalence)"}
# Words that qualify how cases were counted without changing which disease it is.
FILLER = {"incident", "prevalent", "primary", "personal", "history", "risk", "total", "any", "overall", "the", "of",
          "onset", "diagnosed", "diagnosis", "all", "and"}


def _words(s: str | None) -> set[str]:
    return set(re.findall(r"[a-z0-9]+", (s or "").lower().replace("coeliac", "celiac"))) - FILLER


def same_outcome(eval_trait: str | None, score_traits: list[str | None]) -> bool:
    """The evaluation measured the score's own disease, not a subtype, special population or other outcome."""
    w = _words(eval_trait)
    return bool(w) and any(w <= _words(t) for t in score_traits if t)


def _gh(n: int = 40) -> list[tuple[float, float]]:
    """Gauss–Hermite nodes/weights for E[f(Z)], Z ~ N(0,1)."""
    import numpy as np
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return list(zip(x.tolist(), (w / w.sum()).tolist(), strict=True))


GH = _gh()


def _solve(f, lo: float, hi: float, target: float) -> float:
    for _ in range(100):
        mid = (lo + hi) / 2
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def from_auroc(z: float, auc: float, k: float) -> float:
    d = math.sqrt(2) * N.inv_cdf(auc)
    # z is measured against the whole population (a K:(1−K) mixture), so map it onto the non-case scale first.
    s = k * d + z * math.sqrt(1 + k * (1 - k) * d * d)
    a, b = k * N.pdf(s - d), (1 - k) * N.pdf(s)
    return a / (a + b)


def _expit(x: float) -> float:
    return 1 / (1 + math.exp(-x)) if x > -700 else 0.0


def from_or_per_sd(z: float, or_sd: float, k: float) -> float:
    beta = math.log(or_sd)
    a = _solve(lambda a: sum(w * _expit(a + beta * x) for x, w in GH), -40, 40, k)
    return _expit(a + beta * z)


def from_hr_per_sd(z: float, hr_sd: float, k: float) -> float:
    beta = math.log(hr_sd)
    # cumulative hazard H; risk = 1 − exp(−H·e^{βz}); solve H so the mean risk is K
    h = _solve(lambda h: sum(w * (1 - math.exp(-h * math.exp(beta * x))) for x, w in GH), 0, 50, k)
    return 1 - math.exp(-h * math.exp(beta * z))


def _num(s) -> tuple[float | None, float | None, float | None]:
    """'0.633 [0.624,0.641]' → (0.633, 0.624, 0.641)."""
    xs = [float(x) for x in re.findall(r"\d*\.\d+|\d+", str(s or ""))]
    return (xs[0] if xs else None, xs[1] if len(xs) >= 3 else None, xs[2] if len(xs) >= 3 else None)


def technical_only(covariates: str | None) -> bool:
    return all(w in TECHNICAL for w in re.findall(r"[a-z]+", (covariates or "").lower()))


@cache
def table(path: Path = TABLE) -> dict:
    return tomllib.loads(path.read_text()) if path.exists() else {}


def _norm(s: str | None) -> str:
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def baseline(trait_id: str | None, labels: list[str | None], sex: str | None, path: Path = TABLE) -> dict | None:
    """A cited lifetime risk for this disease and sex (falling back to both sexes)."""
    rows = [r for r in table(path).get("baseline", [])
            if r.get("measure") == "lifetime_risk" and risk_stated(r["value"], r.get("quote"))]
    names = {_norm(x) for x in labels if x}
    hits = [r for r in rows if (trait_id and r.get("efo") == trait_id) or _norm(r["disease"]) in names]
    sx = (sex or "").upper()
    want = "male" if sx.startswith("XY") else "female" if sx.startswith("XX") else None
    for s in (want, "both"):
        for r in hits:
            if r.get("sex") == s:
                return r
    return None


def _quoted_numbers(quote: str | None) -> list[float]:
    q = (quote or "").replace("\u00b7", ".")
    return [float(x.replace(",", "")) for x in re.findall(r"\d[\d,]*(?:\.\d+)?", q)]


def risk_stated(value: float, quote: str | None) -> bool:
    """The quote states this risk itself: as a percent, per 1,000 / per 100,000, or '1 in N' (no derived values)."""
    for n in _quoted_numbers(quote):
        if any(math.isclose(n, value * k, rel_tol=1e-6) for k in (100, 1000, 100000)):
            return True
        if n > 1 and math.isclose(1 / n, value, rel_tol=0.02):
            return True
    return False


def effect_stated(estimate: float, quote: str | None) -> bool:
    """The quote states this ratio: '1.47' or as a percent increase per SD ('47% increased risk')."""
    return any(math.isclose(n, estimate, abs_tol=1e-9) or math.isclose(n, round((estimate - 1) * 100), abs_tol=1e-9)
               and f"{n:g}%" in (quote or "").replace(" %", "%") for n in _quoted_numbers(quote))


def catalog_effect(r: dict, evals: list[dict]) -> dict | None:
    """The PGS Catalog record this per-SD row annotates, if the catalog reports the same ratio."""
    for e in evals:
        if e.get("ppm_id") == r.get("ppm_id"):
            est = _num(e.get("hr") if r.get("metric") == "HR" else e.get("or"))[0]
            return e if est is not None and math.isclose(est, r["estimate"], abs_tol=0.005) else None
    return None


def per_sd_problem(r: dict, evals: list[dict], compared_with: str | None) -> str | None:
    """Why this per-SD row can't be used (None = usable): the PGS Catalog record must exist and report the same
    ratio, its cohort ancestry must match yours, the model may add nothing beyond age, sex and technical
    covariates, and the paper must be quoted saying the ratio is per standard deviation."""
    e = catalog_effect(r, evals)
    if not e:
        return "its PGS Catalog record does not report the same ratio"
    if not _ancestry_ok(e, compared_with):
        return f"it was measured in a cohort of different ancestry ({e.get('ancestry')})"
    if not (technical_only_or_age_sex(e.get("covariates")) and technical_only_or_age_sex(r.get("adjusted_for"))):
        return "it was measured with clinical or lifestyle factors in the model"
    if not re.search(r"per (1[- ])?(one )?(SD|standard deviation)|standard deviation increase|each standard deviation",
                     r.get("quote") or "", re.I):
        return "the quoted paper does not say the ratio is per standard deviation"
    return None


def confirmed_per_sd(r: dict, evals: list[dict], compared_with: str | None) -> bool:
    return per_sd_problem(r, evals, compared_with) is None


def per_sd(pgs_id: str, path: Path = TABLE) -> list[dict]:
    return [r for r in table(path).get("per_sd", []) if r.get("pgs_id") == pgs_id]


def _ancestry_ok(e: dict, compared_with: str | None) -> bool:
    anc = (e.get("ancestry") or "").lower()
    if (compared_with or "").upper() == "EUR":
        return "european" in anc and "," not in anc
    return True


def eligible_auroc(evals: list[dict], compared_with: str | None, traits: list[str | None]) -> list[dict]:
    out = []
    for e in evals:
        auc, lo, hi = _num(e.get("auroc"))
        if auc is None or not 0.5 < auc < 1 or not technical_only(e.get("covariates")):
            continue
        if (e.get("n") or 0) < MIN_N or (e.get("cases") or 0) < MIN_CASES:
            continue
        if not same_outcome(e.get("trait_reported"), traits):
            continue
        if not _ancestry_ok(e, compared_with):
            continue
        out.append({**e, "_auc": auc, "_lo": lo, "_hi": hi})
    # independent first, then the largest sample
    return sorted(out, key=lambda e: (bool(e.get("independent")), e.get("n") or 0), reverse=True)


def cohort_baseline(e: dict) -> float | None:
    n, cases = e.get("n"), e.get("cases")
    if not n or not cases or cases < MIN_CASES:
        return None
    k = cases / n
    return k if k <= MAX_CASE_SHARE else None


def _cite(e: dict) -> str:
    return f"{e.get('eval_author') or 'PGS Catalog'} {str(e.get('eval_published') or '')[:4]}".strip()


def _ppm_url(e: dict) -> str | None:
    # PGS Catalog has no per-PPM page; the score page lists every performance record by PPM ID.
    return f"https://www.pgscatalog.org/score/{e['pgs_id']}/" if e.get("pgs_id") else None


def per_100(p: float) -> str:
    v = 100 * p
    return f"{v:.0f}" if v >= 9.5 else f"{v:.1f}".rstrip("0").rstrip(".") if v >= 0.1 else "<0.1"


def estimate(*, pgs_id: str, trait_id: str | None, labels: list[str | None], z: float | None,
             compared_with: str | None, evals: list[dict], sex: str | None, path: Path = TABLE) -> dict:
    """Absolute risk for one score, or {'status': 'not_established', 'reason': ...}."""
    if z is None:
        return {"status": "not_established", "reason": "no ancestry-adjusted score position for you"}
    base = baseline(trait_id, labels, sex, path)
    aucs = eligible_auroc(evals, compared_with, labels)
    all_sds = per_sd(pgs_id, path)
    sds = [r for r in all_sds if confirmed_per_sd(r, evals, compared_with)]

    def kinfo(e: dict | None) -> tuple[float, float | None, float | None, dict] | None:
        if base:
            return (base["value"], base.get("ci_low"), base.get("ci_high"),
                    {"kind": base.get("measure", "lifetime_risk"),
                     "text": f"{MEASURE.get(base.get('measure'), 'lifetime risk')} "
                             f"({LIFETIME_PREFIX.sub('', base['age_window'])}), {base['population']}",
                     "source": base["source"], "url": base["url"], "quote": base.get("quote")})
        k = cohort_baseline(e) if e else None
        if k is None:
            return None
        return (k, None, None, {"kind": "cohort", "text": f"share of the {e.get('cohorts') or 'evaluation'} cohort "
                                f"diagnosed during the study ({e['cases']:,} of {e['n']:,})",
                                "source": _cite(e), "url": _ppm_url(e)})

    for e in aucs:
        k = kinfo(e)
        if not k:
            continue
        kv = k[0]
        p = from_auroc(z, e["_auc"], kv)
        rng = sorted(from_auroc(z, a, kv) for a in (e["_lo"], e["_hi"]) if a and 0.5 < a < 1)
        return {"status": "ok", "method": "auroc", "you": p, "typical": kv,
                "you_low": rng[0] if rng else None, "you_high": rng[-1] if rng else None,
                "you_per_100": per_100(p), "typical_per_100": per_100(kv), "ratio": p / kv,
                "baseline": k[3], "effect": {"metric": "AUROC", "estimate": e["_auc"], "ci": [e["_lo"], e["_hi"]],
                                             "ppm_id": e.get("ppm_id"), "covariates": e.get("covariates"),
                                             "cohort": e.get("cohorts"), "n": e.get("n"), "cases": e.get("cases"),
                                             "url": _ppm_url(e), "source": _cite(e), "pmid": e.get("eval_pmid"),
                                             "trait": e.get("trait_reported")}}
    for r in sds:
        if not base:
            break
        kv = base["value"]
        f = from_or_per_sd if r["metric"] == "OR" else from_hr_per_sd
        p = f(z, r["estimate"], kv)
        rng = sorted(f(z, x, kv) for x in (r.get("ci_low"), r.get("ci_high")) if x)
        k = kinfo(None)
        return {"status": "ok", "method": f"{r['metric'].lower()}_per_sd", "you": p, "typical": kv,
                "you_low": rng[0] if rng else None, "you_high": rng[-1] if rng else None,
                "you_per_100": per_100(p), "typical_per_100": per_100(kv), "ratio": p / kv, "baseline": k[3],
                "effect": {"metric": f"{r['metric']} per SD", "estimate": r["estimate"],
                           "ci": [r.get("ci_low"), r.get("ci_high")], "ppm_id": r.get("ppm_id"),
                           "covariates": r.get("adjusted_for"), "cohort": r.get("cohort"), "url": r["url"],
                           "source": r["source"], "pmid": r.get("pmid"), "quote": r.get("quote")}}
    if not aucs and all_sds and not sds:
        why = ("no published AUROC fits, and the per-standard-deviation effect on file can't be used because "
               + "; ".join(sorted({per_sd_problem(r, evals, compared_with) or "" for r in all_sds})))
    elif not aucs and not sds:
        why = ("no published evaluation of this score reports an AUROC without age, sex or clinical factors mixed in, "
               "or an odds/hazard ratio confirmed to be per standard deviation")
    elif not base:
        why = ("the published evaluations either enrolled extra cases (so their case share is not a population "
               "risk) or are too small, and no cited lifetime risk for this disease is on file")
    else:
        why = "no usable combination of effect size and baseline risk"
    return {"status": "not_established", "reason": why}


def technical_only_or_age_sex(adjusted: str | None) -> bool:
    """Per-SD ratios adjusted for age and sex are still the score's own effect (they remove confounding, not
    add predictors); anything clinical is excluded."""
    s = re.sub(r"\bage at [a-z ]+?(?=[,;(]|$)|\b(age|sex|gender)\b(\s*\(where relevant\))?", "", adjusted or "",
               flags=re.I)
    return technical_only(s)
