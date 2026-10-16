import numpy as np

from wgs.stages import ancestry


def test_mix_recovers_known_weights_and_stays_on_the_simplex():
    centres = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0], [10.0, 10.0], [5.0, -8.0]])
    x = 0.7 * centres[0] + 0.3 * centres[1]
    w = ancestry.mix(centres, x)
    assert abs(w.sum() - 1) < 1e-9 and (w >= 0).all()
    assert np.allclose(centres.T @ w, x, atol=1e-9)
    # A point outside the hull is projected onto it, never given negative weights.
    w = ancestry.mix(centres, np.array([-5.0, -5.0]))
    assert (w >= 0).all() and abs(w.sum() - 1) < 1e-9


def _panel(rng, groups):
    pcs, sp, pop = [], [], []
    for (s, p), centre in groups.items():
        pcs.append(rng.normal(centre, 0.3, size=(40, 10)))
        sp += [s] * 40
        pop += [p] * 40
    return {"pcs": np.vstack(pcs), "superpop": np.array(sp), "population": np.array(pop)}


def test_components_noise_floor_and_nearest_group():
    rng = np.random.default_rng(1)
    z = np.zeros(10)
    groups = {(s, s + "1"): z + 10 * np.eye(10)[i] for i, s in enumerate(ancestry.SOURCES)}
    # Two European populations that overlap completely, one that is distinct.
    groups[("EUR", "EUR2")] = groups[("EUR", "EUR1")]
    groups[("EUR", "EUR3")] = groups[("EUR", "EUR1")] + 3 * np.eye(10)[7]
    ref = _panel(rng, groups)
    rows, cal = ancestry.components(ref, groups[("EUR", "EUR1")])
    eur = next(r for r in rows if r["superpop"] == "EUR")
    assert eur["share"] > 0.95 and eur["detected"]
    assert all(not r["detected"] for r in rows if r["superpop"] != "EUR")
    pops, near = ancestry.nearest(ref, groups[("EUR", "EUR1")], "EUR")
    assert set(near["group"]) == {"EUR1", "EUR2"}  # indistinguishable pair reported together
    assert sum(r["neighbours"] for r in pops) == ancestry.K_NEIGHBOURS
    assert next(r for r in pops if r["code"] == "EUR3")["loo_recall"] > 0.9


def test_grades():
    rows = ancestry.grade_components([
        {"superpop": "EUR", "share": 0.9, "noise_floor": 0.05, "detected": True, "reference_own_share": 0.98},
        {"superpop": "SAS", "share": 0.08, "noise_floor": 0.03, "detected": True, "reference_own_share": 0.94},
        {"superpop": "AFR", "share": 0.02, "noise_floor": 0.03, "detected": False, "reference_own_share": 0.99}],
        pca_variants=60_000)
    assert [r["evidence_level"] for r in rows] == ["Strong", "Moderate", "Limited"]
    assert ancestry._level(0.91, ancestry.MT_HIGH, ancestry.MT_MEDIUM) == "High"
    assert ancestry._level(0.85, ancestry.MT_HIGH, ancestry.MT_MEDIUM) == "Medium"
    assert ancestry._level(None, 1, 1) == "Not callable"
