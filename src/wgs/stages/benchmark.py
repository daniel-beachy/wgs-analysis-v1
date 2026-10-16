"""Accuracy benchmark on Genome in a Bottle's HG002 (ADR-023).

Genome in a Bottle knows HG002's true genotype almost everywhere. HG002's DNBSEQ reads, thinned to 30× like a typical
consumer genome, are run through two of this pipeline's checks and compared with that truth, in the ACMG genes'
exons and at the medicine-defining positions:

- the read re-count (ADR-022): when a call is right, does the re-count confirm it? When a call is wrong (the wrong
  number of copies, or a variant that isn't there), does the re-count catch it?
- bcftools, a second variant caller: what share of its calls are right (precision), and what share of the true
  variants does it find (recall)? Scored with RTG vcfeval, which matches variants written in different but
  equivalent ways (Krusche P, et al. Nat Biotechnol. 2019;37:555-560. PMID 30858580).

Every result is split into ordinary positions and the ones GIAB marks as hard to sequence, by reason. This measures
the methods, not the owner's genome: nothing here reads the owner's data.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import replace
from pathlib import Path

import duckdb
import pandas as pd

from ..knowledge import Store
from ..pipeline import Context, Stage, run, write_json
from . import recount as rc
from . import reference

DECOY_SNVS = 2500      # ClinVar single-letter variants HG002 doesn't carry, used as wrong "you carry one copy" claims
DECOY_INDELS = 1500    # the same for insertions and deletions (sampled with a fixed hash, so reruns pick the same)
DECOY_GAP = 25         # decoys this close to a true variant are skipped: may be the same change written differently
MIN_QUAL = 30          # bcftools calls with a quality below this (1 in 1,000 chance of being wrong) are dropped
COMMON_AF = 0.01       # the usual line between a common and a rare variant: 1 in 100 copies in some population
RTG_MEM = "4g"         # RTG Tools otherwise reserves 90% of the computer's memory
WRONG = {0: [1], 1: [0, 2], 2: [1]}   # wrong claims tested for each true number of copies
REASONS = ["repeat", "segdup", "low_mappability", "gc", "mhc", "other"]
SETTINGS = {"decoy_snvs": DECOY_SNVS, "decoy_indels": DECOY_INDELS, "decoy_gap": DECOY_GAP, "min_qual": MIN_QUAL,
            "common_af": COMMON_AF, "recount": rc.SETTINGS}


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "benchmark"
    p.mkdir(parents=True, exist_ok=True)
    return p


def out_sites(ctx: Context) -> Path:
    return out_dir(ctx) / "hg002_sites.parquet"


def out_calls(ctx: Context) -> Path:
    return out_dir(ctx) / "hg002_calls.parquet"


def out_summary(ctx: Context) -> Path:
    return out_dir(ctx) / "benchmark.json"


def giab_files(st: Store) -> dict[str, Path | None]:
    return {"bam": st.file("giab", "hg002.bam"), "truth": st.file("giab", "truth.vcf.gz"),
            "regions": st.table("giab", "regions"), "confident": st.table("giab", "confident"),
            "difficult": st.table("giab", "difficult")}


# ---------------------------------------------------------------------------------------------------------------
# Test sites
# ---------------------------------------------------------------------------------------------------------------

def eval_intervals(regions: Path, confident: Path) -> list[tuple[str, int, int]]:
    """Where the truth is reliable inside the streamed regions (0-based, half-open)."""
    return duckdb.connect().execute(f"""SELECT r.chrom, greatest(r.start, c.start), least(r."end", c."end")
        FROM '{regions}' r JOIN '{confident}' c ON r.chrom = c.chrom AND c.start < r."end" AND c."end" > r.start
        ORDER BY 1, 2""").fetchall()


def truth_sites(vcf: Path, intervals: list[tuple[str, int, int]]) -> list[dict]:
    """Every true variant allele inside `intervals`, with HG002's true number of copies."""
    import pysam

    out = []
    with pysam.VariantFile(str(vcf)) as vf:
        for ch, s, e in intervals:
            for r in vf.fetch(ch, s, e):
                gt = r.samples[0]["GT"]
                if r.start < s or r.stop > e or gt is None or None in gt:
                    continue
                for i, alt in enumerate(r.alts or (), start=1):
                    if not set(alt) <= set("ACGT") or not set(r.ref) <= set("ACGT"):
                        continue
                    out.append({"chrom": ch, "pos": r.pos, "ref": r.ref, "alt": alt, "truth": sum(a == i for a in gt),
                                "other": sum(a not in (0, i) for a in gt), "source": "truth"})
    return out


def decoy_sites(clinvar: Path, intervals: list[tuple[str, int, int]], truth: list[dict]) -> list[dict]:
    """ClinVar variants HG002 does not carry (absent from the truth inside its reliable regions)."""
    con = duckdb.connect()
    con.register("iv", pd.DataFrame(intervals, columns=["chrom", "start", "end"]))
    con.register("t", pd.DataFrame(truth or [{"chrom": "", "pos": 0, "ref": "N"}])[["chrom", "pos", "ref"]])
    rows = con.execute(f"""WITH c AS (
          SELECT DISTINCT v.chrom, v.pos, v.ref, v.alt, length(v.ref) = 1 AND length(v.alt) = 1 AS snv
          FROM '{clinvar}' v JOIN iv ON v.chrom = iv.chrom AND v.pos - 1 >= iv.start
            AND v.pos - 1 + length(v.ref) <= iv."end"
          WHERE regexp_full_match(v.ref, '[ACGT]+') AND regexp_full_match(v.alt, '[ACGT]+') AND v.ref <> v.alt
            AND NOT EXISTS (SELECT 1 FROM t WHERE t.chrom = v.chrom
              AND t.pos + length(t.ref) + {DECOY_GAP} >= v.pos AND t.pos <= v.pos + length(v.ref) + {DECOY_GAP})),
        ranked AS (SELECT *, row_number() OVER (PARTITION BY snv
                   ORDER BY hash(chrom || ':' || pos || ref || '>' || alt)) AS k FROM c)
        SELECT chrom, pos, ref, alt FROM ranked WHERE (snv AND k <= {DECOY_SNVS}) OR (NOT snv AND k <= {DECOY_INDELS})
        ORDER BY chrom, pos""").fetchall()
    return [{"chrom": c, "pos": p, "ref": r, "alt": a, "truth": 0, "other": 0, "source": "decoy"}
            for c, p, r, a in rows]


def count_sites(bam: Path, fasta: Path, sites: list[dict]) -> pd.DataFrame:
    """Re-count each site, then judge the right claim and each wrong one against the same counts."""
    import pysam

    recs = []
    with pysam.AlignmentFile(str(bam), "rb") as af, pysam.FastaFile(str(fasta)) as fa:
        for s in sites:
            t = rc.Target(f"{s['chrom']}:{s['pos']}:{s['ref']}>{s['alt']}", s["chrom"], s["pos"], s["ref"], s["alt"],
                          expected=s["truth"], other_copies=s["other"], source=s["source"])
            c = rc.count_site(af, fa, f"chr{s['chrom']}", t)
            rec = rc.site_record(t, c, rc.classify(c, t))
            rec["truth"] = s["truth"]
            rec["wrong"] = json.dumps({str(w): rc.classify(c, replace(t, expected=w))["status"]
                                       for w in WRONG[min(s["truth"], 2)]})
            recs.append(rec)
    df = pd.DataFrame.from_records(recs)
    df["flags"] = df["flags"].map(list)
    return df


# ---------------------------------------------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------------------------------------------

def _share(n: int, d: int) -> float | None:
    return round(n / d, 4) if d else None


def _group_kind(k: str) -> str:
    return "snv" if k == "snv" else "indel"


def recount_metrics(df: pd.DataFrame) -> dict:
    """Counts and shares for one group of re-counted sites."""
    right = df["status"].value_counts().to_dict()
    wrong: dict[str, int] = {}
    for w in df["wrong"]:
        for st in json.loads(w).values():
            wrong[st] = wrong.get(st, 0) + 1
    clear = df[df["clear"]]
    n_right, n_wrong = len(df), sum(wrong.values())
    return {"sites": n_right, "right": right, "wrong": wrong,
            "right_confirmed": _share(right.get("confirmed", 0), n_right),
            "right_false_alarm": _share(right.get("disagrees", 0), n_right),
            "wrong_caught": _share(wrong.get("disagrees", 0), n_wrong),
            "wrong_confirmed": _share(wrong.get("confirmed", 0), n_wrong),
            "clear": len(clear), "clear_correct": int((clear["observed"] == clear["truth"]).sum()),
            "clear_accuracy": _share(int((clear["observed"] == clear["truth"]).sum()), len(clear))}


def by_groups(df: pd.DataFrame, metrics) -> dict:
    """metrics() for SNVs and indels × all / easy / difficult, plus each difficulty reason."""
    df = df.assign(group=df["kind"].map(_group_kind), easy=df["difficult"].map(len) == 0)
    out = {"groups": [], "by_reason": []}
    for g in ("all", "snv", "indel"):
        sub = df if g == "all" else df[df["group"] == g]
        for region, part in (("all", sub), ("easy", sub[sub["easy"]]), ("difficult", sub[~sub["easy"]])):
            if len(part):
                out["groups"].append({"kind": g, "region": region, **metrics(part)})
    for reason in REASONS:
        part = df[df["difficult"].map(lambda r, reason=reason: reason in r)]
        if len(part):
            out["by_reason"].append({"reason": reason, **metrics(part)})
    return out


def caller_metrics(df: pd.DataFrame) -> dict:
    """Precision and recall from vcfeval's labels (tp / fp on the calls, tp_base / fn on the truth)."""
    n = df["label"].value_counts().to_dict()
    tp, fp, tpb, fn = (n.get(k, 0) for k in ("tp", "fp", "tp_base", "fn"))
    p, r = _share(tp, tp + fp), _share(tpb, tpb + fn)
    calls = df[df["label"].isin(["tp", "fp"])]
    conf = calls[calls["recount"] == "confirmed"]
    return {"tp": tp, "fp": fp, "tp_base": tpb, "fn": fn, "precision": p, "recall": r,
            "f1": round(2 * p * r / (p + r), 4) if p and r else None,
            "fp_confirmed": int((conf["label"] == "fp").sum()), "tp_confirmed": int((conf["label"] == "tp").sum()),
            "precision_if_confirmed": _share(int((conf["label"] == "tp").sum()), len(conf))}


# ---------------------------------------------------------------------------------------------------------------
# Second caller: bcftools, scored with RTG vcfeval
# ---------------------------------------------------------------------------------------------------------------

def tool_version(cmd: list[str]) -> str | None:
    try:
        return run(cmd, capture=True, check=False).stdout.splitlines()[0].strip()
    except (OSError, IndexError):
        return None


def bed_file(intervals: list[tuple[str, int, int]], path: Path, prefix: str = "chr") -> Path:
    path.write_text("".join(f"{prefix}{c}\t{s}\t{e}\n" for c, s, e in intervals))
    return path


def bcftools_calls(aln: Path, fasta: Path, bed: Path, out: Path, log: Path, threads: int) -> Path:
    """bcftools' multiallelic caller on the reads in `bed`, left-aligned, calls of quality ≥ MIN_QUAL."""
    run(f"bcftools mpileup -f '{fasta}' -R '{bed}' -a FORMAT/AD,FORMAT/DP -q {rc.MIN_MAPQ} -Q {rc.MIN_BASEQ} "
        f"--max-depth 1000 -Ou '{aln}' | bcftools call -m -v -Ou | bcftools norm -f '{fasta}' -Ou | "
        f"bcftools view -i 'QUAL>={MIN_QUAL}' --threads {threads} -Oz -o '{out}'", shell=True, log=log)
    run(["bcftools", "index", "-f", "-t", str(out)], log=log)
    return out


def reference_sdf(ctx: Context, fasta: Path) -> Path:
    """RTG's own copy of the reference (about 1 GB), kept in fast scratch and rebuilt if missing or stale."""
    d = ctx.scratch.parent / "rtg" / f"{fasta.stem}.sdf"
    done = d / "wgs.done"
    stamp = f"{fasta}|{fasta.stat().st_size}|{fasta.stat().st_mtime_ns}"
    if done.exists() and done.read_text() == stamp:
        return d
    shutil.rmtree(d, ignore_errors=True)
    d.parent.mkdir(parents=True, exist_ok=True)
    run(["rtg", "format", "-o", str(d), str(fasta)], env={"RTG_MEM": RTG_MEM}, log=ctx.logs / "benchmark.log")
    done.write_text(stamp)
    return d


def vcfeval(base: Path, calls: Path, sdf: Path, bed: Path, out: Path, log: Path, threads: int,
            squash: bool = False) -> Path:
    shutil.rmtree(out, ignore_errors=True)
    cmd = ["rtg", "vcfeval", "-b", str(base), "-c", str(calls), "-t", str(sdf), "-e", str(bed), "-o", str(out),
           "--ref-overlap", "--vcf-score-field", "QUAL", "-T", str(threads)]
    run(cmd + (["--squash-ploidy"] if squash else []), env={"RTG_MEM": RTG_MEM}, log=log)
    return out


def vcfeval_records(out: Path) -> pd.DataFrame:
    """One row per record in vcfeval's output files: label, position, alleles, and the call's genotype."""
    import pysam

    rows = []
    for label, name in (("tp", "tp.vcf.gz"), ("fp", "fp.vcf.gz"), ("tp_base", "tp-baseline.vcf.gz"),
                        ("fn", "fn.vcf.gz")):
        p = out / name
        if not p.exists():
            continue
        with pysam.VariantFile(str(p)) as vf:
            for r in vf:
                gt = r.samples[0].get("GT") if len(r.samples) else None
                alts = [a for a in (r.alts or ()) if set(a) <= set("ACGT")]
                if not alts:
                    continue
                rows.append({"label": label, "chrom": r.chrom.removeprefix("chr"), "pos": r.pos, "ref": r.ref,
                             "alt": alts[0], "alts": list(r.alts or ()), "qual": r.qual,
                             "copies": sum(a == 1 for a in gt) if gt and None not in gt else None,
                             "other": sum(a not in (0, 1) for a in gt) if gt and None not in gt else 0,
                             "kind": "snv" if all(len(a) == 1 for a in (r.ref, *r.alts)) else "indel"})
    return pd.DataFrame(rows, columns=["label", "chrom", "pos", "ref", "alt", "alts", "qual", "copies", "other",
                                       "kind"])


def recount_calls(bam: Path, fasta: Path, calls: pd.DataFrame, contig) -> list[dict]:
    """Re-count each caller's variant with the call's own genotype as the claim (ADR-022)."""
    import pysam

    out = []
    with pysam.AlignmentFile(str(bam), "rc" if bam.suffix == ".cram" else "rb",
                             reference_filename=str(fasta)) as af, pysam.FastaFile(str(fasta)) as fa:
        for r in calls.to_dict("records"):
            t = rc.Target("", r["chrom"], int(r["pos"]), r["ref"], r["alt"], expected=r["copies"],
                          ploidy=r.get("ploidy", 2), other_copies=int(r["other"] or 0))
            c = rc.count_site(af, fa, contig(r["chrom"]), t)
            v = rc.classify(c, t)
            out.append({"recount": v["status"], "observed": v["observed"], "lod": v["lod"], "reads": v["reads"],
                        "ref_reads": c.ref, "alt_reads": c.alt, "low_mapq_reads": c.low_mapq, "flags": v["flags"]})
    return out


# ---------------------------------------------------------------------------------------------------------------
# Seen in other people (gnomAD, 1000 Genomes)
# ---------------------------------------------------------------------------------------------------------------

POPULATION = ["gnomad_filtered", "common", "rare", "not_seen"]


def population_class(r: dict) -> str:
    """How often a call is seen in other people. gnomAD marks a site whose calls failed its own quality checks
    (gnomAD v2.1.1's random-forest and other filters; Karczewski KJ, et al. Nature. 2020;581:434-443. PMID 32461654);
    such sites usually hold reads from a look-alike copy, in gnomAD's data as in yours."""
    if r.get("gnomad_filtered"):
        return "gnomad_filtered"
    af = max((x for x in (r.get("gnomad_af"), r.get("af_1kg")) if x is not None and x == x), default=None)
    if af is None:
        return "not_seen"
    return "common" if af >= COMMON_AF else "rare"


def read_gnotated(path: Path) -> pd.DataFrame:
    """gnomAD fields slivar added to each record: highest population frequency, homozygous people, filtered flag."""
    import pysam

    def known(x):   # slivar writes -1 when gnomAD has no value: that is "unknown", not "ultra-rare"
        return x if x is not None and x >= 0 else None

    rows = []
    with pysam.VariantFile(str(path)) as vf:
        for r in vf:
            i = r.info
            rows.append({"chrom": r.chrom.removeprefix("chr"), "pos": r.pos, "ref": r.ref, "alt": r.alts[0],
                         "gnomad_af": known(i.get("gnomad_popmax_af")),
                         "gnomad_nhomalt": known(i.get("gnomad_nhomalt")),
                         "gnomad_filtered": "gnomad_popmax_af_filter" in i})
    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "gnomad_af", "gnomad_nhomalt",
                                       "gnomad_filtered"])


def population_metrics(recs: pd.DataFrame) -> list[dict]:
    """Precision of the caller's calls in each population class (how often a call there was real)."""
    calls = recs[recs["label"].isin(["tp", "fp"])]
    out = []
    for p in POPULATION:
        sub = calls[calls["population"] == p]
        if len(sub):
            tp = int((sub["label"] == "tp").sum())
            out.append({"population": p, "calls": len(sub), "tp": tp, "fp": len(sub) - tp,
                        "precision": _share(tp, len(sub))})
    return out


def population(st: Store, calls: pd.DataFrame, work: Path, log: Path) -> pd.DataFrame:
    """Add gnomAD (v2.1.1 exomes + genomes, via slivar) and 1000 Genomes frequencies and the population class."""
    keys = calls[["chrom", "pos", "ref", "alt"]].drop_duplicates()
    found = pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "gnomad_af", "gnomad_nhomalt", "gnomad_filtered"])
    gn = st.table("gnomad", "fields")
    zip_ = gn.parent / "gnomad.hg37.zip" if gn else None
    if len(keys) and zip_ and zip_.exists() and shutil.which("slivar"):
        vcf = work / "population.vcf"
        order = sorted(keys.itertuples(index=False), key=lambda r: (r.chrom.zfill(2), r.pos))
        vcf.write_text("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID=chr{c}>\n" for c in sorted(
            {r.chrom for r in order}, key=lambda c: c.zfill(2))) + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "".join(f"chr{r.chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\t.\t.\n" for r in order))
        run(["slivar", "expr", "--vcf", str(vcf), "--gnotate", str(zip_), "-o", str(work / "population.out.vcf")],
            log=log)
        found = read_gnotated(work / "population.out.vcf")
    out = calls.merge(found, on=["chrom", "pos", "ref", "alt"], how="left")
    kg = st.table("1000genomes", "af")
    if kg:
        con = duckdb.connect()
        con.register("k", keys)
        af = con.execute(f"SELECT k.*, a.af::DOUBLE AS af_1kg FROM k JOIN '{kg}' a USING (chrom, pos, ref, alt)").df()
        out = out.merge(af, on=["chrom", "pos", "ref", "alt"], how="left")
    else:
        out["af_1kg"] = None
    out["gnomad_filtered"] = out["gnomad_filtered"].fillna(False).astype(bool)
    out["population"] = [population_class(r) for r in out.to_dict("records")]
    return out


# ---------------------------------------------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------------------------------------------

def _benchmark(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    g = giab_files(st)
    fasta = reference.reference_fasta(ctx)
    log = ctx.logs / "benchmark.log"
    work = ctx.scratch / "benchmark"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    iv = eval_intervals(g["regions"], g["confident"])
    truth = truth_sites(g["truth"], iv)
    decoys = decoy_sites(st.table("clinvar", "variants"), iv, truth)
    sites = count_sites(g["bam"], fasta, truth + decoys)
    sites["difficult"] = rc.difficult_reasons(sites, g["difficult"])
    sites.to_parquet(out_sites(ctx), index=False)

    summary: dict = {"recount": by_groups(sites, recount_metrics)}
    reads = (st.current("giab") or {}).get("notes", {}).get("reads", {})
    summary["data"] = {"sample": "HG002", "platform": "DNBSEQ (BGI), 2×150", "depth": reads.get("depth"),
                       "full_depth": reads.get("full_depth"), "regions": reads.get("regions"),
                       "streamed_bases": reads.get("bases"), "eval_bases": sum(e - s for _, s, e in iv),
                       "truth_alleles": len(truth), "decoys": len(decoys), "giab": st.versions().get("giab")}

    if shutil.which("rtg"):
        bed = bed_file(iv, work / "eval.bed")
        calls = bcftools_calls(g["bam"], fasta, bed_file(iv, work / "regions.bed"), work / "bcftools.vcf.gz", log,
                               ctx.threads)
        truth_chr = work / "truth.vcf.gz"
        (work / "chroms.txt").write_text("".join(f"{i}\tchr{i}\n" for i in range(1, 23)))
        run(f"bcftools view -R '{bed_file(iv, work / 'eval_nochr.bed', prefix='')}' '{g['truth']}' -Ou | "
            f"bcftools annotate --rename-chrs '{work / 'chroms.txt'}' -Oz -o '{truth_chr}'", shell=True, log=log)
        run(["bcftools", "index", "-f", "-t", str(truth_chr)], log=log)
        ev_out = vcfeval(truth_chr, calls, reference_sdf(ctx, fasta), bed, work / "vcfeval", log, ctx.threads)
        recs = vcfeval_records(ev_out)
        is_call = recs["label"].isin(["tp", "fp"])
        counted = recount_calls(g["bam"], fasta, recs[is_call], lambda ch: f"chr{ch}")
        recs.loc[is_call, "recount"] = [c["recount"] for c in counted]
        recs["difficult"] = rc.difficult_reasons(recs, g["difficult"])
        recs = population(st, recs, work, log)
        recs.to_parquet(out_calls(ctx), index=False)
        summary["caller"] = {"tool": tool_version(["bcftools", "--version"]),
                             "scorer": tool_version(["rtg", "version"]), "min_qual": MIN_QUAL,
                             **by_groups(recs, caller_metrics), "by_population": population_metrics(recs)}
    else:
        summary["caller"] = {"skipped": "RTG Tools not installed (pixi environment 'default' provides it)"}
    summary["settings"] = SETTINGS
    write_json(out_summary(ctx), summary)
    top = next(x for x in summary["recount"]["groups"] if x["kind"] == "all" and x["region"] == "all")
    return {"sites": top["sites"], "right_confirmed": top["right_confirmed"], "wrong_caught": top["wrong_caught"]}


def _inputs(ctx: Context) -> list[Path]:
    st = Store(ctx.cfg)
    return [p for p in (*giab_files(st).values(), st.table("clinvar", "variants"), st.table("gnomad", "fields"),
                        st.table("1000genomes", "af"), reference.reference_fasta(ctx)) if p]


def _available(ctx: Context) -> str | None:
    g = giab_files(Store(ctx.cfg))
    if missing := [k for k, v in g.items() if not v]:
        return f"Genome in a Bottle files missing ({', '.join(missing)}): run `wgs knowledge refresh --source giab`"
    fasta = reference.reference_fasta(ctx)
    if not fasta:
        return "no verified GRCh37 reference (reference stage)"
    if "chr1" not in reference.fai_contigs(fasta):
        return "the reference names chromosomes without 'chr'; the HG002 reads need 'chr1'…"
    return None


STAGE = Stage(
    name="benchmark", version="2", title="Accuracy check on a reference genome (GIAB HG002)", fn=_benchmark,
    inputs=_inputs, outputs=lambda ctx: [out_sites(ctx), out_summary(ctx)], available=_available,
    params=lambda ctx: {"settings": SETTINGS},
)
