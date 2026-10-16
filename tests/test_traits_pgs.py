import csv
import json
from pathlib import Path

from wgs import evidence as ev
from wgs.knowledge.pgs import select
from wgs.stages import pgs
from wgs.stages.pgs import _pvar_header, best_evaluation, effect_text
from wgs.stages.traits import effect_copies, genotype_share, gwas_effect


def test_effect_copies_reads_forward_then_complement():
    assert effect_copies(["A", "G"], "A", {"A", "G"}) == (1, "forward")
    assert effect_copies(["T", "T"], "A", {"T", "C"}) == (2, "reverse")
    assert effect_copies(["A", "G"], None, {"A", "G"}) == (None, None)
    assert effect_copies(None, "A", {"A"}) == (None, None)


def test_genotype_share_is_hardy_weinberg():
    assert genotype_share(0, 0.1) == (0.9 ** 2)
    assert abs(genotype_share(1, 0.1) - 0.18) < 1e-12
    assert abs(genotype_share(2, 0.1) - 0.01) < 1e-12
    assert genotype_share(None, 0.1) is None


def test_gwas_effect_reads_catalog_convention():
    beta = gwas_effect(0.23, "[0.21-0.25] unit decrease", "T", "Hair shape")
    assert beta["effect_kind"] == "beta" and beta["effect_direction"] == -1 and "down" in beta["effect_text"]
    odds = gwas_effect(1.744, "NR", "T", "Red hair")
    assert odds["effect_kind"] == "or" and odds["effect_direction"] == 1 and "raises the odds" in odds["effect_text"]
    assert gwas_effect(0.6, "[0.5-0.7]", "G", "x")["effect_direction"] == -1
    assert gwas_effect(float("nan"), "NR", "G", "x")["effect_text"] is None
    assert gwas_effect(1.2, "NR", None, "x")["effect_text"] is None


def test_gwas_evidence_needs_replication_significance_and_direction():
    assert ev.gwas_evidence(publications=3, mlog10p=50, direction_known=True)[0] == "Strong"
    assert ev.gwas_evidence(publications=2, mlog10p=50, direction_known=True)[0] == "Moderate"
    assert ev.gwas_evidence(publications=5, mlog10p=6, direction_known=True)[0] == "Limited"
    assert ev.gwas_evidence(publications=5, mlog10p=50, direction_known=False)[0] == "Limited"


def test_pgs_grades_never_strong_and_need_match():
    assert ev.pgs_evidence(independent_pubs=9, independent_n=10**7, evaluated_n=10**7)[0] == "Moderate"
    assert ev.pgs_evidence(independent_pubs=1, independent_n=10**6, evaluated_n=0)[0] == "Limited"
    assert ev.pgs_evidence(independent_pubs=0, independent_n=0, evaluated_n=5000)[0] == "Limited"
    assert ev.pgs_call(match_rate=0.95, passed=True)[0] == "High"
    assert ev.pgs_call(match_rate=0.8, passed=True)[0] == "Medium"
    assert ev.pgs_call(match_rate=0.6, passed=False)[0] == "Not callable"
    assert ev.overall("Moderate", "Not callable") == "Not callable"


def test_best_evaluation_prefers_independent_then_effect_then_size():
    evals = [{"independent": False, "or": 1.5, "n": 10**6},
             {"independent": True, "n": 10**6, "other": "R2 0.1"},
             {"independent": True, "auroc": 0.62, "n": 5000}]
    b = best_evaluation(evals)
    assert b["auroc"] == 0.62 and effect_text(b) == "AUROC 0.62"
    assert effect_text(evals[1]) == "R2 0.1"
    assert best_evaluation([]) is None and effect_text(None) is None


def test_pvar_header_counts_meta_lines(tmp_path):
    p = tmp_path / "x.pvar"
    p.write_text("##fileformat=PVARv1.0\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\n1\t10\trs1\tA\tG\n")
    assert _pvar_header(p) == (3, ["CHROM", "POS", "ID", "REF", "ALT"])


def _write(p, rows):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def test_select_prefers_independently_tested_and_honours_pins(tmp_path):
    def score(k, trait, pgp, date="2020-01-01"):
        return {"Polygenic Score (PGS) ID": k, "Mapped Trait(s) (EFO ID)": trait, "Mapped Trait(s) (EFO label)": trait,
                "PGS Publication (PGP) ID": pgp, "Number of Interaction Terms": "0", "Release Date": date,
                "Number of Variants": "100"}
    _write(tmp_path / "pgs_all_metadata_scores.csv", [
        score("PGS1", "T1", "P1"), score("PGS2", "T1", "P2"), score("PGS3", "T2", "P3"),
        score("PGS4", "T1,T2", "P4")])
    _write(tmp_path / "pgs_all_metadata_evaluation_sample_sets.csv",
           [{"PGS Sample Set (PSS)": "S1", "Number of Individuals": "1000"}])
    _write(tmp_path / "pgs_all_metadata_performance_metrics.csv", [
        {"Evaluated Score": "PGS2", "PGS Sample Set (PSS)": "S1", "PGS Publication (PGP) ID": "P9"},
        {"Evaluated Score": "PGS3", "PGS Sample Set (PSS)": "S1", "PGS Publication (PGP) ID": "P3"}])
    rows = select(tmp_path, [{"trait_id": "T1", "label": "Trait one", "section": "traits"}])
    # T2's only score lacks independent tests: not browsable. PGS1 rides along as T1's stand-in.
    assert [(r["pgs_id"], r["alternate_for"]) for r in rows] == [("PGS2", ""), ("PGS1", "T1")]
    assert rows[0]["independent_pubs"] == 1 and rows[0]["featured"] and not rows[1]["featured"]
    pinned = select(tmp_path, [{"trait_id": "T1", "pin": "PGS1"}])
    assert pinned[0]["pgs_id"] == "PGS1" and pinned[0]["pinned"] and pinned[0]["rule_pick"] == "PGS2"
    assert pinned[1]["pgs_id"] == "PGS2" and pinned[1]["alternate_for"] == "T1"


def _score(pid, group, pgp, pct, z, featured=False, pubs=0, risk=None, passed=True):
    return {"pgs_id": pid, "trait_group": group, "pgp_id": pgp, "percentile": pct, "z": z, "passed": passed,
            "featured": featured, "label": group, "section": "health", "independent_pubs": pubs,
            "independent_n": 0, "evaluated_n": 0, "abs_risk": json.dumps(risk) if risk else None,
            "overall": "Limited", "evidence_level": "Limited", "call_confidence": "High"}


def test_trait_headline_is_median_over_studies_not_scores():
    # Two scores from one paper (P1) count once: median of P1's 10/30 is 20, then median of (20, 90) = 55.
    rows = [_score("A", "T", "P1", 10.0, -1.3, featured=True), _score("B", "T", "P1", 30.0, -0.5),
            _score("C", "T", "P2", 90.0, 1.3), _score("D", "T", "P3", None, None, passed=False)]
    (t,) = pgs.summarise(rows)
    assert t["percentile"] == 55.0 and t["n_studies"] == 2 and t["n_scored"] == 3 and t["n_scores"] == 4
    assert (t["low"], t["high"]) == (10.0, 90.0) and t["featured"]
    assert t["agreement"] == "disagree"


def test_trait_risk_from_most_tested_score_with_a_number():
    ok = {"status": "ok", "you": 0.1}
    rows = [_score("A", "T", "P1", 50.0, 0.0, pubs=9, risk={"status": "not_established", "reason": "x"}),
            _score("B", "T", "P2", 60.0, 0.2, pubs=3, risk=ok), _score("C", "T", "P3", 70.0, 0.5, pubs=1, risk=ok)]
    (t,) = pgs.summarise(rows)
    assert t["risk_pgs_id"] == "B" and t["best_pgs_id"] == "A"


def test_agreement_categories():
    assert pgs.agreement([], 0)[0] == "none"
    assert pgs.agreement([2.0], 1)[0] == "single"
    assert pgs.agreement([0.2, 1.1], 2)[0] == "consistent"
    assert pgs.agreement([0.4, 2.4, 1.4], 3)[0] == "same_direction"
    assert pgs.agreement([-3.3, -1.2], 2)[0] == "same_direction"
    assert pgs.agreement([-0.1, 1.6], 2)[0] == "disagree"
    # Two methods from one paper that disagree are not hidden; the text says they are one study.
    kind, text = pgs.agreement([-0.7, 1.2], 1)
    assert kind == "disagree" and "one study" in text


def test_files_skips_appledouble(tmp_path):
    for n in ["self_ALL_additive_0.scorefile.gz", "._self_ALL_additive_0.scorefile.gz", "b.scorefile.gz"]:
        (tmp_path / n).write_bytes(b"")
    assert [p.name for p in pgs._files(tmp_path, "*.scorefile.gz")] == ["b.scorefile.gz",
                                                                         "self_ALL_additive_0.scorefile.gz"]


def test_match_batches_cap_and_cover(tmp_path):
    files = []
    for i, size in enumerate([300, 250, 100, 60, 40, 10]):
        f = tmp_path / f"s{i}.gz"
        f.write_bytes(b"x" * size)
        files.append(str(f))
    batches = pgs.match_batches(files, cap=320)
    assert sorted(f for b in batches for f in b) == sorted(files)
    assert all(sum(Path(f).stat().st_size for f in b) <= 320 or len(b) == 1 for b in batches)
    assert len(batches) == 3
