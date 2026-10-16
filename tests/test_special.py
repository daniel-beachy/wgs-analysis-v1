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


def _exons(ratios, gcs=None, gene="G"):
    gcs = gcs or [0.45] * len(ratios)
    return [{"gene": gene, "chrom": "1", "start": i * 1000, "end": i * 1000 + 100, "exon": i + 1, "ratio": r,
             "ratio_all": r, "gc": g} for i, (r, g) in enumerate(zip(ratios, gcs, strict=True))]


def test_gc_correction_removes_gc_driven_dips():
    # 40 normal mid-GC exons plus 40 high-GC exons that all read at ~0.6: a GC effect, not a deletion.
    ex = _exons([1.0 + (i % 5) * 0.02 for i in range(40)] + [0.6 + (i % 5) * 0.02 for i in range(40)],
                [0.45] * 40 + [0.70] * 40)
    special.gc_correct(ex)
    assert not special.depth_runs(ex)


def test_real_deletion_survives_gc_correction():
    ratios = [1.0 + (i % 7) * 0.01 for i in range(60)]
    ratios[20:23] = [0.5, 0.48, 0.52]          # one copy lost over three neighbouring exons
    ex = _exons(ratios)
    special.gc_correct(ex)
    runs = special.depth_runs(ex)
    assert len(runs) == 1 and runs[0]["svtype"] == "DEL" and runs[0]["exons"] == [21, 22, 23]


def test_extreme_gc_and_unmappable_exons_give_no_signal():
    e = {"gc": 0.80, "gc_ratio": 0.2, "gc_ratio_all": 0.2, "gc_z": -9}
    assert special.depth_kind(e) is None
    # Dip only among uniquely placed reads (all reads normal): mappability, not a deletion.
    assert special.depth_kind({"gc": 0.45, "gc_ratio": 0.3, "gc_ratio_all": 1.0, "gc_z": -9}) is None


def test_cnv_classification():
    base = {"svtype": "DEL", "exons": ["G exon 2"], "gnomad": None}
    assert special.classify_cnv({**base, "delly": True, "depth_support": True})[0] == "possible_exon_deletion"
    assert special.classify_cnv({**base, "delly": False, "depth_support": True})[0] == "unconfirmed_signal"
    assert special.classify_cnv({**base, "delly": True, "depth_support": True,
                                 "gnomad": {"af": 0.2, "popmax_af": 0.3}})[0] == "common"
    assert special.classify_cnv({**base, "exons": [], "delly": True, "depth_support": False})[0] == "no_exon"
    many = {**base, "exons": [f"G exon {i}" for i in range(6)], "delly": False, "depth_support": True}
    assert special.classify_cnv(many)[0] == "long_depth_signal"


def test_adult_onset_uses_typical_onset():
    # Huntington's: rare juvenile cases (earliest onset 2) but typically 30+; must count as adult-onset.
    assert special.adult_onset({"typ_age_onset_min": 30, "age_onset_min": 2})
    assert not special.adult_onset({"typ_age_onset_min": 1, "age_onset_min": 0})
    assert special.adult_onset({"typ_age_onset_min": None, "age_onset_min": 40})
    assert not special.adult_onset({})


def test_provider_copies_reads_indel_letters():
    """tellmeGen writes the ABO O deletion as D / I; GRCh37's reference lacks the base, so I = non-O."""
    non_o = {"ref": "-", "alt": "C"}
    assert special.provider_copies("DD", non_o) == 0
    assert special.provider_copies("DI", non_o) == 1
    b = {"ref": "G", "alt": "T"}
    assert special.provider_copies("GG", b) == 0 and special.provider_copies("GT", b) == 1
    assert special.provider_copies("--", b) is None and special.provider_copies(None, b) is None


def test_panel_exons_breaks_ties_the_same_way_every_time(tmp_path, monkeypatch):
    """TNNT2 has no MANE Select match in Ensembl r87 and several transcripts with the same number of exons; the pick
    was arbitrary, so the exon list (and the HG002 regions built from it) changed between runs."""
    import gzip

    import pandas as pd

    genes = tmp_path / "genes.parquet"
    pd.DataFrame([{"symbol": "TNNT2", "gene_id": "ENSG1", "chrom": "1", "start": 100, "end": 900}]).to_parquet(genes)
    mane = tmp_path / "mane.parquet"
    pd.DataFrame([{"transcript": "ENST99", "status": "MANE Select"}]).to_parquet(mane)
    gff = tmp_path / "a.gff3.gz"
    lines = []
    for tx, starts in (("ENST00000509001", (100, 300)), ("ENST00000367318", (110, 310))):
        lines.append(f"1\t.\tmRNA\t100\t900\t.\t+\t.\tID=transcript:{tx};Parent=gene:ENSG1;biotype=protein_coding")
        lines += [f"1\t.\texon\t{s}\t{s + 50}\t.\t+\t.\tParent=transcript:{tx};rank={i}"
                  for i, s in enumerate(starts, 1)]
    with gzip.open(gff, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    acmg = tmp_path / "acmg.csv"
    acmg.write_text("gene,inheritance,condition\nTNNT2,AD,Cardiomyopathy\n")
    monkeypatch.setattr(special, "ACMG_LIST", acmg)

    class St:
        def table(self, src, name):
            return {"genes": genes, "transcripts": mane}[name]

        def file(self, src, name):
            return gff

    picks = {e["transcript"] for _ in range(5) for e in special.panel_exons(St())}
    assert picks == {"ENST00000367318"}
