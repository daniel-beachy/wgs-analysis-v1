"""GIAB HG002 benchmark and second caller (ADR-023): intervals, scoring, verdicts, hard-to-sequence regions."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from wgs import evidence as ev
from wgs.knowledge import giab
from wgs.stages import benchmark as bm
from wgs.stages import recount as rc
from wgs.stages import secondcall as sc


def test_merge_joins_overlapping_and_touching():
    iv = [("2", 10, 20), ("1", 5, 10), ("1", 0, 5), ("1", 30, 40), ("1", 35, 50), ("10", 0, 1)]
    assert giab.merge(iv) == [("1", 0, 10), ("1", 30, 50), ("2", 10, 20), ("10", 0, 1)]


def test_read_bed_strips_chr(tmp_path):
    p = tmp_path / "x.bed"
    p.write_text("#header\nchr1\t5\t9\textra\n2\t1\t3\n")
    assert giab.read_bed(p) == [("1", 5, 9), ("2", 1, 3)]


def test_difficult_reasons_overlap(tmp_path):
    d = tmp_path / "d.parquet"
    pd.DataFrame([{"chrom": "1", "start": 99, "end": 110, "reason": "repeat"},
                  {"chrom": "1", "start": 100, "end": 200, "reason": "segdup"}]).to_parquet(d)
    sites = pd.DataFrame([{"chrom": "1", "pos": 100, "ref": "A"},     # 0-based 99: inside the repeat only
                          {"chrom": "1", "pos": 105, "ref": "AT"},    # inside both
                          {"chrom": "1", "pos": 99, "ref": "A"},      # 0-based 98: just before both
                          {"chrom": "2", "pos": 105, "ref": "A"}])
    assert rc.difficult_reasons(sites, d) == [["repeat"], ["repeat", "segdup"], [], []]
    assert rc.difficult_reasons(sites, None) == [[], [], [], []]


def _site(kind, truth, status, wrong, observed, clear, difficult=()):
    return {"kind": kind, "truth": truth, "status": status, "wrong": json.dumps(wrong), "observed": observed,
            "clear": clear, "difficult": list(difficult)}


def test_recount_metrics_counts_right_and_wrong_claims():
    df = pd.DataFrame([
        _site("snv", 1, "confirmed", {"0": "disagrees", "2": "disagrees"}, 1, True),
        _site("snv", 0, "confirmed", {"1": "disagrees"}, 0, True),
        _site("indel", 2, "disagrees", {"1": "confirmed"}, 1, True, ["repeat"]),   # the re-count got it wrong
        _site("indel", 1, "weak", {"0": "unclear", "2": "unclear"}, 1, False, ["repeat"]),
    ])
    m = bm.recount_metrics(df)
    assert m["sites"] == 4 and m["right_confirmed"] == 0.5 and m["right_false_alarm"] == 0.25
    assert m["wrong"] == {"disagrees": 3, "confirmed": 1, "unclear": 2}
    assert m["wrong_caught"] == 0.5 and m["clear_accuracy"] == pytest.approx(2 / 3, abs=1e-4)
    g = bm.by_groups(df, bm.recount_metrics)
    easy = next(x for x in g["groups"] if x["kind"] == "all" and x["region"] == "easy")
    assert easy["sites"] == 2 and easy["right_confirmed"] == 1.0
    assert [r["reason"] for r in g["by_reason"]] == ["repeat"]


def test_caller_metrics_precision_recall():
    df = pd.DataFrame([{"label": "tp", "recount": "confirmed"}] * 8 + [{"label": "fp", "recount": "confirmed"}]
                      + [{"label": "fp", "recount": "disagrees"}] + [{"label": "tp_base", "recount": None}] * 8
                      + [{"label": "fn", "recount": None}] * 2)
    m = bm.caller_metrics(df)
    assert (m["precision"], m["recall"]) == (0.8, 0.8)
    assert m["precision_if_confirmed"] == pytest.approx(8 / 9, abs=1e-4) and m["fp_confirmed"] == 1


@pytest.mark.parametrize("reads,lod,observed,flags,want", [
    (30, 8.0, 1, [], "reads_support"), (30, 8.0, 1, ["look_alike"], "reads_flagged"),
    (30, 8.0, 0, [], "reads_disagree"), (5, 1.0, 1, [], "unclear"),
])
def test_second_caller_verdict(reads, lod, observed, flags, want):
    assert sc.verdict({"reads": reads, "lod": lod, "observed": observed, "flags": flags}) == want


def test_difficult_reason_is_explained_not_graded():
    base = {"ref": "A", "alt": "G", "ref_reads": 15, "alt_reads": 14, "expected": 1, "observed": 1, "flags": [],
            "status": "confirmed"}
    level, why = ev.recount_call("High", [], {**base, "difficult": ["mhc", "other", "repeat"]})
    assert level == "High"
    assert "HLA" in why[-1] and "repeated" in why[-1] and "other reasons" not in why[-1]


@pytest.mark.parametrize("row,want", [
    ({"gnomad_filtered": True, "gnomad_af": 0.5, "af_1kg": 0.6}, "gnomad_filtered"),
    ({"gnomad_filtered": False, "gnomad_af": 0.002, "af_1kg": 0.03}, "common"),
    ({"gnomad_filtered": False, "gnomad_af": 0.002, "af_1kg": float("nan")}, "rare"),
    ({"gnomad_filtered": False, "gnomad_af": None, "af_1kg": None}, "not_seen"),
])
def test_population_class(row, want):
    assert bm.population_class(row) == want


def test_read_gnotated_spots_gnomad_filter_flag(tmp_path):
    p = tmp_path / "g.vcf"
    p.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr19>\n"
                 '##INFO=<ID=gnomad_popmax_af,Number=1,Type=Float,Description="x">\n'
                 '##INFO=<ID=gnomad_popmax_af_filter,Number=0,Type=Flag,Description="x">\n'
                 '##INFO=<ID=gnomad_nhomalt,Number=1,Type=Integer,Description="x">\n'
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "chr19\t100\t.\tT\tC\t.\t.\tgnomad_popmax_af=0.49;gnomad_popmax_af_filter;gnomad_nhomalt=0\n"
                 "chr19\t200\t.\tG\tA\t.\t.\tgnomad_popmax_af=0.2;gnomad_nhomalt=30\n"
                 "chr19\t300\t.\tG\tA\t.\t.\t.\n"
                 "chr19\t400\t.\tG\tA\t.\t.\tgnomad_popmax_af=-1;gnomad_nhomalt=-1\n")   # slivar: not in gnomAD
    g = bm.read_gnotated(p)
    assert g["chrom"].tolist() == ["19"] * 4 and g["gnomad_filtered"].tolist() == [True, False, False, False]
    assert g["gnomad_af"].iloc[0] == pytest.approx(0.49) and g["gnomad_af"].iloc[2:].isna().all()
    assert bm.population_class(g.iloc[3].to_dict() | {"af_1kg": None}) == "not_seen"


def test_population_metrics_precision_per_class():
    recs = pd.DataFrame([{"label": "fp", "population": "gnomad_filtered"}] * 3
                        + [{"label": "tp", "population": "gnomad_filtered"}]
                        + [{"label": "tp", "population": "common"}] * 4
                        + [{"label": "fn", "population": "common"}])
    m = {r["population"]: r for r in bm.population_metrics(recs)}
    assert m["gnomad_filtered"]["precision"] == 0.25 and m["common"]["calls"] == 4 and "rare" not in m


def test_viewer_hard_region_words_match_python():
    import re
    from pathlib import Path

    ts = (Path(__file__).parents[1] / "viewer/src/lib/evidence.ts").read_text()
    block = ts[ts.index("export const DIFFICULT"):].split("};", 1)[0]
    words = {k: v.replace("’", "'") for k, v in re.findall(r"^\s*(\w+): '(.+)',$", block, re.M)}
    assert words == ev.DIFFICULT_WORD
