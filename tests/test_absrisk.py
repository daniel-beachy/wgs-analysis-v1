import math
import re

import pytest

from wgs import absrisk as ar


@pytest.mark.parametrize("f,arg", [(ar.from_auroc, 0.7), (ar.from_or_per_sd, 1.6), (ar.from_hr_per_sd, 1.6)])
def test_models_are_calibrated_to_the_population_risk(f, arg):
    for k in (0.01, 0.1, 0.3):
        mean = sum(w * f(x, arg, k) for x, w in ar.GH)
        assert math.isclose(mean, k, rel_tol=0.01)
        assert f(2, arg, k) > k > f(-2, arg, k)


def test_auroc_half_means_no_information():
    assert math.isclose(ar.from_auroc(2.5, 0.5, 0.07), 0.07)


def test_covariate_gate():
    assert ar.technical_only(None) and ar.technical_only("Study, PCs(1-10)")
    assert ar.technical_only("first four genetic PCs, genotyping array")
    assert not ar.technical_only("age, PC 1-9")
    assert not ar.technical_only("family history of cancer (in first-degree relatives), genotyping array")
    assert ar.technical_only_or_age_sex("age, sex, 4 PCs") and not ar.technical_only_or_age_sex("age, sex, BMI")


def test_case_enriched_cohort_is_not_a_baseline():
    assert ar.cohort_baseline({"n": 400812, "cases": 4340}) == pytest.approx(4340 / 400812)
    assert ar.cohort_baseline({"n": 1950, "cases": 974}) is None
    assert ar.cohort_baseline({"n": 5000, "cases": 20}) is None


def test_estimate_uses_eligible_auroc_and_cohort_baseline(tmp_path):
    evals = [{"ppm_id": "PPM1", "auroc": "0.7 [0.69,0.71]", "covariates": "age, PC 1-9", "n": 9000, "cases": 900,
              "ancestry": "European", "independent": True},
             {"ppm_id": "PPM2", "auroc": "0.63 [0.62,0.64]", "covariates": "Genotyping array", "n": 400000,
              "cases": 4000, "ancestry": "European", "independent": True, "cohorts": "UKB",
              "trait_reported": "Incident X disease"}]
    r = ar.estimate(pgs_id="PGSX", trait_id="EFO_X", labels=["x disease"], z=2.0, compared_with="EUR", evals=evals,
                    sex="XY", path=tmp_path / "none.toml")
    assert r["status"] == "ok" and r["effect"]["ppm_id"] == "PPM2" and r["baseline"]["kind"] == "cohort"
    assert r["you_low"] < r["you"] < r["you_high"] and r["you"] > r["typical"]
    r = ar.estimate(pgs_id="PGSX", trait_id="EFO_X", labels=["x disease"], z=2.0, compared_with="EUR", evals=evals[:1],
                    sex="XY", path=tmp_path / "none.toml")
    assert r["status"] == "not_established" and "AUROC" in r["reason"]


def test_outcome_must_be_the_scores_disease():
    t = ["breast carcinoma", "Breast cancer"]
    assert ar.same_outcome("Incident primary breast cancer", t)
    assert not ar.same_outcome("Breast cancer intrinsic-like subtype (triple negative)", t)
    assert not ar.same_outcome("Coronary artery disease in childhood cancer survivors", ["Coronary artery disease"])
    assert ar.same_outcome("Coeliac disease", ["celiac disease"])


def test_cited_lifetime_risk_preferred(tmp_path):
    p = tmp_path / "t.toml"
    p.write_text('[[baseline]]\ndisease = "x disease"\nefo = "EFO_X"\nsex = "male"\nmeasure = "lifetime_risk"\n'
                 'value = 0.2\nage_window = "to 85"\npopulation = "P"\nsource = "S"\nurl = "https://u"\n'
                 'quote = "lifetime risk was 20%"\n'
                 '[[per_sd]]\npgs_id = "PGSX"\nppm_id = "PPMX"\nquote = "HR per standard deviation"\n'
                 'metric = "HR"\nestimate = 1.5\nci_low = 1.4\nci_high = 1.6\n'
                 'adjusted_for = "age, sex, 10 PCs"\nsource = "S2"\nurl = "https://v"\n')
    ar.table.cache_clear()
    r = ar.estimate(pgs_id="PGSX", trait_id="EFO_X", labels=[], z=1.0, compared_with="EUR",
                    evals=[{"ppm_id": "PPMX", "hr": "1.5 [1.4,1.6]", "ancestry": "European", "covariates": "age"}],
                    sex="XY", path=p)
    assert r["method"] == "hr_per_sd" and r["typical"] == 0.2 and r["baseline"]["kind"] == "lifetime_risk"
    assert ar.baseline("EFO_X", [], "XX", p) is None


def test_curated_quote_numbers_must_appear_at_source():
    from wgs.audit import numbers_present
    quote = "risk was 12.5% (95% CI 11.0-14.1) in 1,234 men"
    assert numbers_present(quote, "12.5 % (95 % CI 11.0 to 14.1) among 1 234") == []
    assert numbers_present("lifetime risk 1 in 8", "lifetime risk 1 in 9") == ["8"]


def test_vendored_table_rows_state_their_numbers():
    ar.table.cache_clear()
    t = ar.table()
    for r in t["baseline"]:
        assert r.get("url") and r.get("quote") and r.get("source"), r["disease"]
        if r["measure"] == "lifetime_risk":
            assert ar.risk_stated(r["value"], r["quote"]), r["disease"]
    for r in t["per_sd"]:
        assert r.get("ppm_id") and r.get("quote") and r.get("url"), r["pgs_id"]
        assert re.search(r"per (1[- ])?(one )?(SD|standard deviation)|standard deviation increase|each standard "
                         r"deviation", r["quote"], re.I), r["pgs_id"]


def test_per_sd_needs_matching_catalog_record():
    r = {"pgs_id": "P", "ppm_id": "PPM1", "metric": "HR", "estimate": 1.43, "adjusted_for": "age, sex",
         "quote": "HR per one standard deviation (SD) increase"}
    e = {"ppm_id": "PPM1", "hr": "1.43 [1.36,1.49]", "ancestry": "European", "covariates": "age, sex, PCs(1-10)"}
    assert ar.confirmed_per_sd(r, [e], "EUR")
    assert not ar.confirmed_per_sd(r, [{**e, "hr": "1.53 [1.4,1.6]"}], "EUR")
    assert not ar.confirmed_per_sd(r, [{**e, "ancestry": "African unspecified"}], "EUR")
    assert not ar.confirmed_per_sd(r, [{**e, "covariates": "age, sex, BMI, smoking"}], "EUR")
    assert not ar.confirmed_per_sd(r, [], "EUR")
    assert ar.risk_stated(0.402, "was 40·2% (95% CI 39·2-41·3)") and ar.risk_stated(0.0072, "7.2 per 1,000")
    assert not ar.risk_stated(0.044, "1.0% for BP-I, 1.1% for BP-II, and 2.4% for subthreshold")
