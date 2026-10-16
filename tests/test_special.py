from wgs.stages import special

T = {"benign_max": 26, "intermediate_min": 27, "pathogenic_min": 36}


def test_pure_run_handles_interruptions_and_iupac():
    ar = "TG" + "CAG" * 22 + "CAAGAGACTAGCCCCAGG" + "CAG" * 5
    assert special.pure_run(ar, "CAG") == 22
    assert special.pure_run("GCGGCAGCC" * 2, "GCN") == 6


def test_clinical_scale_is_relative_to_the_reference():
    # AR: the catalogue region counts 35 units for a reference with 22 pure CAG; +2 units → 24 clinical.
    assert special.to_clinical(37, 35, 22) == 24
    # A loose motif can over-count the reference; it is capped at the top of the typical range.
    assert special.to_clinical(15, 15, 15, benign_max=14) == 14
    # Insertion-type loci are counted directly.
    assert special.to_clinical(400, 1, 0) == 400


def test_classify_repeat():
    assert special.classify_repeat([17, 18], [17, 18], 17, 17, T) == ("typical", "High")
    assert special.classify_repeat([17, 40], [17, 44], 17, 17, T) == ("disease_range", "Medium")
    assert special.classify_repeat([17, 30], [17, 31], 17, 17, T) == ("intermediate", "Medium")
    # Most likely typical, but the interval reaches the disease range: Low call, never a silent all-clear.
    assert special.classify_repeat([17, 18], [17, 60], 17, 17, T) == ("typical", "Low")
    assert special.classify_repeat([0, 0], [714, 714], 17, 17, T) == ("not_readable", "Not callable")
    # Inconsistent source ranges are not used.
    assert special.classify_repeat([5, 5], [5, 5], 5, 4, {"benign_max": 4, "pathogenic_min": 3})[0] == "no_thresholds"
    # A second caller (the VCF) agreeing on reference length, with ExpansionHunter one unit off: typical, High.
    assert special.classify_repeat([16, 16], [16, 16], 15, 15, {"benign_max": 14, "pathogenic_min": 15},
                                   vcf_same=True) == ("typical", "High")
    # ...but not when ExpansionHunter sees two different alleles further apart (HTT 15/18 vs reference 17).
    assert special.classify_repeat([15, 18], [15, 18], 17, 17, T, vcf_same=True) == ("typical", "High")


def test_smn_rows():
    amb = special.smn_row({"SMN1": None, "SMN2": None, "SMN1_CN_raw": [1.4, 1.5], "Info": "Ambiguous"}, uneven=True)
    assert amb["status"] == "not_determined" and amb["call_level"] == "Not callable" and amb["next_step"]
    assert special.smn_row({"SMN1": 1, "SMN2": 2}, uneven=False)["status"] == "carrier"
    assert special.smn_row({"SMN1": 2, "SMN2": 2, "g27134TG_CN": 1}, uneven=False)["status"] == \
        "possible_silent_carrier"
    assert special.smn_row({"SMN1": 2, "SMN2": 1}, uneven=False)["call_level"] == "High"
    claims = special.gene_claims([special.smn_row({"SMN1": 1, "SMN2": 2}, uneven=False)])
    assert claims and claims[0]["section"] == "carrier" and claims[0]["role"] == "carrier"


def test_gba_rows():
    clean = special.gba_row({"deletion_CN": 2, "variant_raw_count": {"c.1226A>G": "0,30"}}, uneven=False)
    assert clean["status"] == "typical" and clean["call_level"] == "High"
    hit = special.gba_row({"deletion_CN": 2, "GBAP1_like_variant_exon9_11": "", "other_variants": "c.1226A>G",
                           "variant_raw_count": {"c.1226A>G": "14,16"}}, uneven=False)
    assert hit["status"] == "carrier"


def test_hla_matching_rules(monkeypatch):
    table = [
        {"id": "ro", "condition": "Coeliac disease", "kind": "rule_out", "alleles": ["DQB1*02", "DQB1*03:02"],
         "markers": "DQ2 or DQ8", "quotes": [], "caveats": [{"if_present": ["DQB1*03:01"], "quote": "dq7"}]},
        {"id": "hap", "condition": "Type 1 diabetes", "kind": "haplotype", "direction": "protective",
         "alleles": ["DRB1*15:01", "DQB1*06:02"], "markers": "DR15-DQ6", "quotes": []},
        {"id": "risk", "condition": "Narcolepsy", "kind": "risk", "alleles": ["DQB1*06:02"], "markers": "x",
         "quotes": []},
    ]
    monkeypatch.setattr(special, "hla_table", lambda: table)
    typed = {"DQB1": ["06:02", "03:01"], "DRB1": ["15:01", "04:01"]}
    rows = {r["association_id"]: r for r in special.hla_rows(typed, {"DQB1": [36, 43], "DRB1": [36, 29]})}
    assert rows["ro"]["notable"] and "don't carry" in rows["ro"]["statement"]
    assert "dq7" in rows["ro"]["quotes"]  # caveat shown because DQB1*03:01 is present
    assert rows["hap"]["notable"] and rows["hap"]["direction"] == "protective"
    assert rows["risk"]["notable"] and "narcolepsy" in rows["risk"]["statement"]
    # Field-boundary matching: DQB1*02 must not match DQB1*020x-style names or other loci.
    assert special.carries({"DQB1": ["02:01", "03:02"]}, "DQB1*02") == ["DQB1*02:01"]
    assert special.carries({"DQB1": ["03:01"]}, "DQB1*03:02") == []
    # A rule needing a gene that wasn't typed says nothing.
    assert not [r for r in special.hla_rows({"DRB1": ["15:01"]}, {}) if r["association_id"] == "ro"]


def test_curated_first_author_check():
    from wgs.audit import first_author_matches
    assert first_author_matches("Kallweit U, et al. J Clin Sleep Med. 2022", {"authors": [{"name": "Kallweit U"}]})
    assert not first_author_matches("Ziemann LC, et al. J Clin Sleep Med. 2022", {"authors": [{"name": "Kallweit U"}]})


def test_abo_rules():
    assert special.abo_type(0, 0, 0)[0] == "O"
    assert special.abo_type(1, 0, 0)[0] == "A"
    assert special.abo_type(1, 1, 0)[0] == "B"
    assert special.abo_type(2, 1, 0)[0] == "AB"
    assert special.abo_type(2, 0, 1)[1].endswith("(A2 subgroup)")
    assert special.abo_type(None, 0, 0)[0] is None          # unreadable position: no guess
    assert special.abo_type(0, 1, 0)[0] is None             # B markers on an O backbone: needs a lab


def test_rhd_rules():
    t = {"zero_below": 0.15, "one_between": [0.3, 0.7], "two_above": 0.8}
    assert special.rhd_type(0.02, t)[0] == 0
    assert special.rhd_type(0.5, t)[0] == 1
    assert special.rhd_type(1.0, t)[0] == 2
    assert special.rhd_type(0.2, t)[0] is None and special.rhd_type(None, t)[0] is None
