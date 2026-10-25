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


def _genotypes(rng, freqs, n):
    return rng.binomial(2, np.clip(freqs, 0.02, 0.98), size=(n, len(freqs))).astype(float)


def test_finer_comparison_pools_panels_and_flags_a_batch_effect():
    rng = np.random.default_rng(3)
    base = rng.uniform(0.1, 0.9, 3000)
    tsi = np.clip(base + rng.normal(0, 0.12, 3000), 0, 1)
    fin = np.clip(base + rng.normal(0, 0.12, 3000), 0, 1)

    def panel(tuscan_hgdp):
        x = np.vstack([_genotypes(rng, tsi, 60), _genotypes(rng, tuscan_hgdp, 15), _genotypes(rng, fin, 60)])
        labels = np.array(["TSI"] * 60 + ["TuscanHGDP"] * 15 + ["FIN"] * 60)
        return x, labels

    x, labels = panel(tsi)  # the same population sampled twice
    x[5, 7] = np.nan
    reg = ancestry.regional(x, labels, _genotypes(rng, fin, 1)[0], "EUR")
    assert reg["near"]["nearest"] == "FIN" and reg["near"]["positions"] > 2500
    assert reg["near"]["pcs"] in ancestry.REGION_PCS and len(reg["my_pcs"]) == max(ancestry.REGION_PCS)
    assert reg["near"]["twins"] == {"TuscanHGDP": {"twin": "TSI", "confused": True}}
    assert sum(r["neighbours"] for r in reg["rows"]) == ancestry.K_NEIGHBOURS
    # A technical shift (like low- vs high-depth genotypes) separates the twins, and the check must say so.
    shifted = np.clip(tsi + np.where(rng.random(3000) < 0.3, 0.35, 0), 0, 1)
    x, labels = panel(shifted)
    reg = ancestry.regional(x, labels, _genotypes(rng, fin, 1)[0], "EUR")
    assert reg["near"]["twins"]["TuscanHGDP"]["confused"] is False


def test_panel_alleles_are_aligned_to_the_grch37_site(tmp_path):
    from wgs.knowledge import ancestry_panels as ap

    assert ap.align(("A", "G", "+"), "G", "A") == 1       # a1 is the alt allele
    assert ap.align(("A", "G", "+"), "A", "G") == -1      # a1 is the ref allele
    assert ap.align(("A", "G", "-"), "C", "T") == 1       # lifted to the minus strand: complement first
    assert ap.align(("A", "T", "+"), "T", "A") == 0       # A/T can't be oriented across builds
    assert ap.align(("A", "G", "+"), "A", "C") == 0
    g = np.array([2, 1, 0, ap.MISSING], dtype=np.uint8)
    assert ap.alt_copies(g, -1).tolist() == [0, 1, 2, ap.MISSING]
    # Two sites x three people; PLINK 1 packs four people per byte, low bits first.
    (tmp_path / "x.bim").write_text("1\t1:10:A:G\t0\t10\tG\tA\n1\t1:20:C:T\t0\t20\tT\tC\n")
    (tmp_path / "x.fam").write_text("a a 0 0 1 -9\nb b 0 0 2 -9\nc c 0 0 1 -9\n")
    (tmp_path / "x.bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0b00011000, 0b00001111]))
    bim, fam, geno = ap.read_bed(tmp_path / "x")
    assert fam == ["a", "b", "c"] and len(bim) == 2
    assert geno.tolist() == [[2, 1, ap.MISSING], [0, 0, 2]]


def test_igsr_rows_include_hgdp_populations():
    from wgs.knowledge.ancestry import rows

    recs = [{"code": "YRI", "name": "Yoruba", "dataCollections": [{"title": "1000 Genomes phase 3 release"}],
             "superpopulation": {"code": "AFR"}},
            {"code": None, "elasticId": "YorubaHGDP", "name": "Yoruba", "description": "Yoruba in Nigeria (HGDP)",
             "dataCollections": [{"title": "Human Genome Diversity Project"}], "superpopulation": {"code": None}},
            {"code": None, "elasticId": "SomeSGDP", "dataCollections": [{"title": "Simons Genome Diversity Project"}]}]
    out = {r["code"]: r for r in rows(recs)}
    assert set(out) == {"YRI", "YorubaHGDP"}
    assert out["YorubaHGDP"]["panel"] == "hgdp" and out["YorubaHGDP"]["url"].endswith("/YorubaHGDP")
    assert ancestry.code({"population": "Yoruba", "panel": "hgdp"}) == "YorubaHGDP"
    assert ancestry.code({"population": "YRI", "panel": "1kg30x"}) == "YRI"


def test_populations_confused_through_a_chain_are_reported_together():
    # Han is mixed up with CHB, CHB with CHS, CHS with Dai; JPT with nobody.
    c = {"HAN": {"HAN": 0.5, "CHB": 0.3}, "CHB": {"CHB": 0.7, "CHS": 0.25}, "CHS": {"CHS": 0.4, "DAI": 0.1},
         "DAI": {"DAI": 0.1, "CHS": 0.7}, "JPT": {"JPT": 1.0}}
    conf = {a: {b: c[a].get(b, 0.0) for b in c} for a in c}
    assert ancestry.confusable(conf, "HAN") == {"HAN", "CHB", "CHS", "DAI"}
    assert ancestry.confusable(conf, "JPT") == {"JPT"}
