"""Re-count the reads behind every finding (ADR-022).

The genotypes in the VCF come from tellmeGen's variant caller. This stage goes back to the raw reads in the CRAM at
every position a finding rests on, plus every position that defines a medicine-gene result, and counts read by
read which version each one shows. A deliberately simple genotype model (each read is a coin flip between the two
copies, with a small error rate; Li 2011) turns the counts into 0, 1 or 2 copies independently of the caller's own
model, and the finding's call confidence is lowered when the two disagree. It is never raised.

What this does and doesn't check: it catches a caller misreading the reads it had (wrong number of copies, a call
the reads don't support, a variant missed at a position we rely on). It cannot catch errors in the reads themselves;
the HG002 benchmark measures that side.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

from .. import evidence as ev
from ..knowledge import Store
from ..pipeline import Context, Stage, write_json
from . import annotate, pgx, qc, reference, special

MIN_MAPQ = 20          # reads placed with less certainty than 1 in 100 are counted separately, not used
MIN_BASEQ = 20         # a base called with less certainty than 1 in 100 is not used for single-letter variants
ERROR = 0.02           # chance a read shows the other version by error (sequencing plus alignment)
MIN_READS = 8          # fewer usable reads than this cannot confirm a genotype on their own
MIN_LOD = 2.0          # best genotype must be 100× likelier than the next to count as a clear answer
LOW_MAPQ_SHARE = 0.30  # this share of reads placed uncertainly marks a look-alike (paralogous) region
STRAND_MIN_ALT = 6     # strand-bias test needs this many variant reads
STRAND_P = 1e-3        # Fisher exact p below which the variant reads sit suspiciously on one strand
PAD = 3                # reference bases kept on each side of an indel's repeat span
MAX_SPAN = 100         # indel windows wider than this can't be spanned by 150-base reads reliably
FLANK = 300            # reference fetched around an indel to find its repeat span
FIT_BASEQ = 10         # bases below this quality are ignored when fitting a whole read
FIT_SPREAD = 12        # how far (bases) from its aligned start a read may be laid when fitting it
FIT_MARGIN = 2         # a read must fit one version by at least this many fewer mismatches to be counted
FIT_MAX_MISMATCH = 8   # a read with more mismatches than this against both versions shows something else
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800   # unmapped, secondary, QC-fail, duplicate, supplementary
PAR = {"X": [(60_001, 2_699_520), (154_931_044, 155_260_560)]}   # GRCh37 pseudo-autosomal regions

SETTINGS = {"min_mapq": MIN_MAPQ, "min_baseq": MIN_BASEQ, "error": ERROR, "min_reads": MIN_READS, "min_lod": MIN_LOD,
            "low_mapq_share": LOW_MAPQ_SHARE, "strand_p": STRAND_P, "max_span": MAX_SPAN, "fit_margin": FIT_MARGIN,
            "fit_max_mismatch": FIT_MAX_MISMATCH}
STATUSES = ["confirmed", "weak", "unclear", "disagrees", "not_called", "not_checkable"]


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "recount"
    p.mkdir(parents=True, exist_ok=True)
    return p


def out_sites(ctx: Context) -> Path:
    return out_dir(ctx) / "recount.parquet"


def out_claims(ctx: Context) -> Path:
    return out_dir(ctx) / "claims.parquet"


def out_pgx_genes(ctx: Context) -> Path:
    return out_dir(ctx) / "pgx_genes.parquet"


def out_summary(ctx: Context) -> Path:
    return out_dir(ctx) / "recount.json"


def claim_sources(ctx: Context) -> list[Path]:
    """The per-stage claim files, before re-counting."""
    return [p for p in (annotate.claims_path(ctx), pgx.out_paths(ctx)["claims"], special.out_claims(ctx))
            if p.exists()]


def fresh(ctx: Context, out: Path, sources: list[Path]) -> bool:
    """True when `out` was written after every file it was derived from (else publish uses the originals)."""
    return out.exists() and all(out.stat().st_mtime >= s.stat().st_mtime for s in sources if s.exists())


# ---------------------------------------------------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class Target:
    key: str
    chrom: str                 # without "chr"
    pos: int                   # 1-based, VCF style
    ref: str
    alt: str
    expected: int | None       # copies of `alt` the call says you have; None when the call has no genotype
    ploidy: int = 2
    other_copies: int = 0      # copies of a third version (a second alternative allele) in the call
    source: str = "claim"
    gene: str | None = None
    rsid: str | None = None
    why_uncalled: str | None = None   # when `expected` is None: why the variant file gave no usable genotype


@dataclass
class Counts:
    ref: int = 0
    alt: int = 0
    other: int = 0
    ref_fwd: int = 0
    alt_fwd: int = 0
    low_mapq: int = 0
    low_baseq: int = 0
    not_spanning: int = 0
    note: str | None = None    # why the site could not be counted, if it couldn't
    window: tuple[int, int] | None = field(default=None)


def kind(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "snv"
    if len(ref) == len(alt):
        return "mnv"
    return "indel"


def indel_span(seq: str, offset: int, ref: str, alt: str) -> tuple[int, int] | None:
    """Range of reference positions (0-based, in `seq` coordinates) over which an indel can slide unchanged.

    `offset` is the 0-based index of the VCF anchor base in `seq`. An indel inside a repeat (say one extra "A" in a
    run of A's) can be written at any point of the run, and aligners don't always pick the same point as the
    caller, so reads are compared over the whole run."""
    if not ref or not alt or ref[0] != alt[0]:
        return None
    if len(ref) > len(alt):                      # deletion of seq[a:b]
        a = offset + len(alt)
        b = a + len(ref) - len(alt)
        while b < len(seq) and seq[a] == seq[b]:
            a, b = a + 1, b + 1
        hi = b
        a = offset + len(alt)
        b = a + len(ref) - len(alt)
        while a - 1 >= 0 and seq[a - 1] == seq[b - 1]:
            a, b = a - 1, b - 1
        return a, hi
    ins = alt[len(ref):]                          # insertion before seq[a]
    a = offset + len(ref)
    s = ins
    while a < len(seq) and seq[a] == s[0]:
        s, a = s[1:] + s[0], a + 1
    hi = a
    a, s = offset + len(ref), ins
    while a - 1 >= 0 and seq[a - 1] == s[-1]:
        s, a = s[-1] + s[:-1], a - 1
    return a, hi


def haplotypes(seq: str, seq_start: int, t: Target) -> tuple[int, int, str, str] | None:
    """(window start, window end, reference version, variant version) for comparing whole read segments."""
    off = t.pos - 1 - seq_start
    if seq[off:off + len(t.ref)].upper() != t.ref.upper():
        return None
    span = indel_span(seq, off, t.ref, t.alt) if kind(t.ref, t.alt) == "indel" else (off, off + len(t.ref))
    if span is None:
        span = (off, off + len(t.ref))
    lo, hi = min(span[0], off) - PAD, max(span[1], off + len(t.ref)) + PAD
    if lo < 0 or hi > len(seq) or hi - lo > MAX_SPAN:
        return None
    ref_hap = seq[lo:hi].upper()
    alt_hap = (seq[lo:off] + t.alt + seq[off + len(t.ref):hi]).upper()
    return seq_start + lo, seq_start + hi, ref_hap, alt_hap


def best_fit(seq: str, quals, hap: str, start: int, spread: int) -> int | None:
    """Fewest mismatches with the whole read laid (without gaps) along `hap`, its first base anywhere within
    `spread` of `start`. Bases called with low confidence don't count."""
    q = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    lo, hi = max(0, start - spread), min(len(hap) - len(q), start + spread)
    if hi < lo:
        return None
    h = np.frombuffer(hap.upper().encode(), dtype=np.uint8)[lo:hi + len(q)]
    good = np.asarray(quals) >= FIT_BASEQ
    return int(((sliding_window_view(h, len(q)) != q) & good).sum(axis=1).min())


def fit_read(read, ref_flank: str, alt_flank: str, flank_start: int, shift: int) -> str | None:
    """'ref' or 'alt' when the whole read (soft-clipped ends included) fits one version clearly better than the
    other, 'other' when it fits neither, None when it doesn't reach far enough across the variant to tell."""
    seq, quals = read.query_sequence, read.query_qualities
    if seq is None or quals is None:
        return None
    lead = read.cigartuples[0][1] if read.cigartuples and read.cigartuples[0][0] == 4 else 0
    est = read.reference_start - lead - flank_start
    spread = abs(shift) + FIT_SPREAD
    fr, fa = best_fit(seq, quals, ref_flank, est, spread), best_fit(seq, quals, alt_flank, est, spread)
    if fr is None or fa is None:
        return None
    if min(fr, fa) > FIT_MAX_MISMATCH:
        return "other"
    if fa <= fr - FIT_MARGIN:
        return "alt"
    if fr <= fa - FIT_MARGIN:
        return "ref"
    return None


def count_site(af, fasta, contig: str, t: Target) -> Counts:
    """Count reads showing each version at one site. `af` is an open pysam AlignmentFile, `fasta` a FastaFile."""
    c = Counts()
    k = kind(t.ref, t.alt)
    if k == "snv":
        p0 = t.pos - 1
        if fasta.fetch(contig, p0, p0 + 1).upper() != t.ref.upper():
            c.note = "the reference letter here doesn't match the call"
            return c
        c.window = (p0, p0 + 1)
        for read in af.fetch(contig, p0, p0 + 1):
            if read.flag & SKIP_FLAGS:
                continue
            if read.mapping_quality < MIN_MAPQ:
                c.low_mapq += 1
                continue
            hit = next(((q, r) for q, r in read.get_aligned_pairs(matches_only=True) if r == p0), None)
            if hit is None:                       # a deletion or a splice across the site
                c.other += 1
                continue
            base, qual = read.query_sequence[hit[0]], read.query_qualities[hit[0]]
            if qual < MIN_BASEQ:
                c.low_baseq += 1
            elif base == t.ref.upper():
                c.ref += 1
                c.ref_fwd += not read.is_reverse
            elif base == t.alt.upper():
                c.alt += 1
                c.alt_fwd += not read.is_reverse
            else:
                c.other += 1
        return c
    start = max(0, t.pos - 1 - FLANK)
    seq = fasta.fetch(contig, start, t.pos - 1 + len(t.ref) + FLANK)
    h = haplotypes(seq, start, t)
    if h is None:
        if seq[t.pos - 1 - start:t.pos - 1 - start + len(t.ref)].upper() != t.ref.upper():
            c.note = "the reference letters here don't match the call"
        elif abs(len(t.ref) - len(t.alt)) + 1 > MAX_SPAN:
            c.note = (f"the change is {abs(len(t.ref) - len(t.alt)):,} letters long, too long for reads to show it "
                      "whole; read depth or long reads are needed")
        else:
            c.note = "the variant sits in a repeat too long for single reads to cross"
        return c
    ws, we, _, _ = h
    c.window = (ws, we)
    off = t.pos - 1 - start
    alt_flank = seq[:off] + t.alt + seq[off + len(t.ref):]
    shift = len(t.alt) - len(t.ref)
    for read in af.fetch(contig, max(0, ws - FIT_SPREAD), we + FIT_SPREAD):   # soft-clipped ends can reach in
        if read.flag & SKIP_FLAGS:
            continue
        if read.mapping_quality < MIN_MAPQ:
            c.low_mapq += 1
            continue
        a = fit_read(read, seq, alt_flank, start, shift)
        if a is None:
            c.not_spanning += 1
        elif a == "ref":
            c.ref += 1
            c.ref_fwd += not read.is_reverse
        elif a == "alt":
            c.alt += 1
            c.alt_fwd += not read.is_reverse
        else:
            c.other += 1
    return c


# ---------------------------------------------------------------------------------------------------------------
# Genotype from counts
# ---------------------------------------------------------------------------------------------------------------

def likelihoods(k: int, n: int, ploidy: int) -> list[float]:
    """log10 likelihood of k variant reads out of n for 0..ploidy copies (binomial; Li 2011)."""
    ps = [ERROR, 1 - ERROR] if ploidy == 1 else [ERROR, 0.5, 1 - ERROR]
    return [k * math.log10(p) + (n - k) * math.log10(1 - p) for p in ps]


def call(k: int, n: int, ploidy: int) -> tuple[int, float]:
    """(most likely copies, how many times likelier than the runner-up, as log10)."""
    ll = likelihoods(k, n, ploidy)
    order = sorted(range(len(ll)), key=lambda i: -ll[i])
    return order[0], ll[order[0]] - ll[order[1]]


def fisher_p(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact test on [[a, b], [c, d]]."""
    n, r1, c1 = a + b + c + d, a + b, a + c

    def prob(x: int) -> float:
        return math.comb(r1, x) * math.comb(n - r1, c1 - x) / math.comb(n, c1)

    p0 = prob(a)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    return min(1.0, sum(p for x in range(lo, hi + 1) if (p := prob(x)) <= p0 * (1 + 1e-9)))


def classify(c: Counts, t: Target) -> dict:
    """Turn read counts into a verdict on the call (ADR-022)."""
    n = c.ref + c.alt + (c.other if t.other_copies else 0)
    flags: list[str] = []
    seen = c.ref + c.alt + c.other + c.low_mapq
    if seen and c.low_mapq >= 3 and c.low_mapq / seen >= LOW_MAPQ_SHARE:
        flags.append("look_alike")
    if c.alt >= STRAND_MIN_ALT and c.ref >= STRAND_MIN_ALT:
        p = fisher_p(c.ref_fwd, c.ref - c.ref_fwd, c.alt_fwd, c.alt - c.alt_fwd)
        if p < STRAND_P:
            flags.append("one_strand")
    out = {"reads": n, "observed": None, "lod": None, "flags": flags}
    if c.note:
        return {**out, "status": "not_checkable", "note": c.note}
    if n == 0:
        return {**out, "status": "not_checkable", "note": "no reads could be placed here with confidence"}
    best, lod = call(c.alt, n, t.ploidy)
    out.update(observed=best, lod=round(lod, 1))
    if t.expected is None:
        return {**out, "status": "not_called", "note": t.why_uncalled or "the variant file has no genotype here"}
    clear = lod >= MIN_LOD and n >= MIN_READS
    if best == t.expected:
        return {**out, "status": "confirmed" if clear else "weak"}
    return {**out, "status": "disagrees" if clear else "unclear"}


# ---------------------------------------------------------------------------------------------------------------
# Which sites
# ---------------------------------------------------------------------------------------------------------------

def ploidy_of(chrom: str, pos: int, sex: str | None) -> int:
    if chrom == "Y":
        return 1
    if chrom == "X" and sex and sex.startswith("XY"):
        return 2 if any(a <= pos <= b for a, b in PAR["X"]) else 1
    return 2


def _gt_alleles(gt: str | None, alleles: list[str]) -> list[str] | None:
    if not gt:
        return None
    idx = [x for x in gt.replace("|", "/").split("/")]
    if any(i == "." for i in idx):
        return None
    return [alleles[int(i)] for i in idx]


def _other_copies(rows: pd.DataFrame, ref: str, alt: str) -> int:
    """Copies of a different alternative allele at the same position (a split multi-allelic site such as A/C at a
    site where the reference is G): those reads count against neither version, so they are kept in the total."""
    same = rows[(rows.ref == ref) & (rows.alt == alt)]
    if not len(same) or not bool(same.iloc[0].multiallelic):
        return 0
    sib = rows[(rows.ref != ref) | (rows.alt != alt)]
    return sum(str(g).replace("|", "/").split("/").count("1") for g in sib["gt"] if g and "." not in str(g))


def claim_targets(claims: pd.DataFrame, variants: pd.DataFrame, sex: str | None) -> list[Target]:
    """One target per distinct variant a claim rests on; expected copies are the claim's own."""
    out: dict[str, Target] = {}
    rows = claims.dropna(subset=["chrom", "pos", "ref", "alt"])
    for r in rows.itertuples():
        chrom, pos, ref, alt = str(r.chrom), int(r.pos), str(r.ref), str(r.alt)
        key = f"{chrom}:{pos}:{ref}:{alt}"
        if key in out:
            continue
        if chrom == "MT":
            out[key] = Target(key, chrom, pos, ref, alt, None, source="claim", gene=r.gene, rsid=r.rsid)
            continue
        copies = getattr(r, "copies", None)
        expected = None if copies is None or pd.isna(copies) else int(copies)
        other = _other_copies(variants[(variants.chrom == chrom) & (variants.pos == pos)], ref, alt)
        ploidy = ploidy_of(chrom, pos, sex)
        if expected is not None and ploidy == 1:
            expected = min(expected, 1)
        out[key] = Target(key, chrom, pos, ref, alt, expected, ploidy, other, "claim", r.gene, r.rsid)
    return list(out.values())


UNCALLED = {"low_quality": "the variant file's call was too uncertain to use",
            "no_call": "the variant file has no genotype",
            "not_callable": "the variant file had too few good reads to make a call"}


def pgx_targets(positions: pd.DataFrame, fasta, contig_name, sex: str | None) -> list[Target]:
    """PharmCAT's defining positions, re-expressed against GRCh37 (where the reference letter can differ)."""
    out = []
    for r in positions.itertuples():
        chrom = str(r.chrom37).removeprefix("chr")
        if not r.pos37 or chrom == "M":
            continue
        pos = int(r.pos37)
        alleles = [str(r.ref), *[str(a) for a in (r.alts if r.alts is not None else [])]]
        key = f"pgx:{chrom}:{pos}:{r.rsid or ''}"
        called = _gt_alleles(r.gt, alleles)
        ref37 = fasta.fetch(contig_name(chrom), pos - 1, pos - 1 + len(alleles[0])).upper()
        if ref37 in alleles:
            others = [a for a in alleles if a != ref37]
            ref = ref37
        elif len(alleles[0]) == 1 and len(ref37) == 1:   # GRCh37 carries a letter PharmCAT doesn't list
            others, ref = alleles, ref37
        else:
            continue
        alt = others[0] if others else None
        if alt is None:
            continue
        expected = None if called is None else sum(a == alt for a in called)
        other = 0 if called is None else sum(a not in (ref, alt) for a in called)
        ploidy = ploidy_of(chrom, pos, sex)
        if expected is not None and ploidy == 1:
            expected = min(expected, 1)
        out.append(Target(key, chrom, pos, ref, alt, expected, ploidy, other, "pgx", r.gene, r.rsid,
                          UNCALLED.get(str(r.status))))
    return out


# ---------------------------------------------------------------------------------------------------------------
# Applying the verdicts
# ---------------------------------------------------------------------------------------------------------------

def site_record(t: Target, c: Counts, v: dict) -> dict:
    return {"key": t.key, "source": t.source, "gene": t.gene, "rsid": t.rsid, "chrom": t.chrom, "pos": t.pos,
            "ref": t.ref, "alt": t.alt, "kind": kind(t.ref, t.alt), "ploidy": t.ploidy, "expected": t.expected,
            "other_copies": t.other_copies, "ref_reads": c.ref, "alt_reads": c.alt, "other_reads": c.other,
            "alt_forward": c.alt_fwd, "ref_forward": c.ref_fwd, "low_mapq_reads": c.low_mapq,
            "low_baseq_reads": c.low_baseq, "not_spanning_reads": c.not_spanning, "reads": v["reads"],
            "observed": v["observed"], "lod": v["lod"], "status": v["status"], "flags": v["flags"],
            "clear": bool(v["reads"] >= MIN_READS and (v["lod"] or 0) >= MIN_LOD), "note": v.get("note")}


def adjust_claims(claims: pd.DataFrame, sites: pd.DataFrame, gene_verdicts: dict[str, dict]) -> pd.DataFrame:
    """Add the re-count to each claim and lower its call confidence where the reads say so."""
    by_key = {r["key"]: r for r in sites.to_dict("records")}
    out = claims.copy()
    recs, levels, reasons, overall = [], [], [], []
    for r in out.to_dict("records"):
        rc = None
        if r.get("chrom") is not None and not pd.isna(r.get("pos")) and r.get("ref") and r.get("alt"):
            rc = by_key.get(f"{r['chrom']}:{int(r['pos'])}:{r['ref']}:{r['alt']}")
        elif r.get("section") == "pgx" and r.get("group") in ("pgx_gene", "pgx_drug"):
            rc = _worst_gene([gene_verdicts.get(g.strip()) for g in str(r.get("gene") or "").split(",")])
        level, why = ev.recount_call(r["call_level"], _as_list(r.get("call_reasons")), rc)
        recs.append(json.dumps(_public(rc)) if rc else None)
        levels.append(level)
        reasons.append(why)
        overall.append(ev.overall(r["evidence_level"], level) if r["evidence_level"] in ev.EVIDENCE else r["overall"])
    out["call_level"], out["call_reasons"], out["overall"] = levels, reasons, overall
    out["recount"] = recs
    out["recount_status"] = [json.loads(x)["status"] if x else None for x in recs]
    return out


def difficult_reasons(sites: pd.DataFrame, difficult: Path | None) -> list[list[str]]:
    """GIAB's reasons each position is hard to sequence (ADR-023), one list per row of `sites` (chrom, pos, ref)."""
    if difficult is None or sites.empty:
        return [[] for _ in range(len(sites))]
    s = pd.DataFrame({"chrom": sites["chrom"].astype(str).to_numpy(), "pos": sites["pos"].astype(int).to_numpy(),
                      "span": sites["ref"].fillna("N").str.len().clip(lower=1).to_numpy(), "i": range(len(sites))})
    con = duckdb.connect()
    con.register("s", s)
    got = dict(con.execute(f"""SELECT s.i, list(DISTINCT d.reason ORDER BY d.reason) FROM s JOIN '{difficult}' d
        ON d.chrom = s.chrom AND d.start < s.pos - 1 + s.span AND d."end" > s.pos - 1 GROUP BY 1""").fetchall())
    return [list(got.get(i, [])) for i in range(len(sites))]


def _as_list(x) -> list:
    return [] if x is None else list(x)


def _public(rc: dict) -> dict:
    keep = ("status", "ref", "alt", "ref_reads", "alt_reads", "other_reads", "low_mapq_reads", "expected", "observed",
            "ploidy", "lod", "flags", "note", "positions", "problems", "kind", "uncalled", "difficult")
    return {k: (list(v) if k in ("flags", "problems", "difficult") and v is not None else v) for k, v in rc.items()
            if k in keep}


STATUS_RANK = {"disagrees": 0, "unclear": 1, "weak": 2, "confirmed": 3}


def _worst_gene(verdicts: list[dict | None]) -> dict | None:
    vs = [v for v in verdicts if v]
    return min(vs, key=lambda v: STATUS_RANK.get(v["status"], 9)) if vs else None


def gene_verdicts(sites: pd.DataFrame, genes: pd.DataFrame) -> dict[str, dict]:
    """One verdict per medicine gene from all of its defining positions. Genes typed by a specialised caller that
    reads the raw reads itself (Cyrius, T1K) are left as they are."""
    out: dict[str, dict] = {}
    src = dict(zip(genes.gene, genes.call_source, strict=True)) if len(genes) else {}
    for gene, g in sites[sites.source == "pgx"].groupby("gene"):
        if not gene or not str(src.get(gene, "")).startswith("PharmCAT"):
            continue
        judged = g[g.status.isin(STATUS_RANK)]
        if not len(judged):
            continue
        worst = min(judged.status, key=lambda s: STATUS_RANK[s])
        problems = judged[judged.status.isin(("disagrees", "unclear", "weak"))]
        flags = sorted({f for fl in judged["flags"] for f in (fl if fl is not None else [])})
        uncalled = [{"where": p.rsid or f"{p.chrom}:{p.pos}", "why": p.note, "alt": p.alt,
                     "alt_reads": int(p.alt_reads), "reads": int(p.reads),
                     "observed": None if pd.isna(p.observed) else int(p.observed),
                     "clear": bool(p.clear)}
                    for p in g[g.status == "not_called"].itertuples()]
        out[gene] = {"status": worst, "positions": len(judged), "flags": flags, "kind": "gene", "uncalled": uncalled,
                     "problems": [f"{p.rsid or f'{p.chrom}:{p.pos}'} ({p.status})" for p in problems.itertuples()]}
    return out


def adjust_genes(genes: pd.DataFrame, verdicts: dict[str, dict]) -> pd.DataFrame:
    out = genes.copy()
    levels, reasons, overall, recs = [], [], [], []
    for r in out.to_dict("records"):
        rc = verdicts.get(r["gene"])
        level, why = ev.recount_call(r["call_level"], _as_list(r.get("call_reasons")), rc)
        levels.append(level)
        reasons.append(why)
        overall.append(ev.overall(r["evidence_level"], level) if r["evidence_level"] in ev.EVIDENCE else r["overall"])
        recs.append(json.dumps(_public(rc)) if rc else None)
    out["call_level"], out["call_reasons"], out["overall"], out["recount"] = levels, reasons, overall, recs
    return out


# ---------------------------------------------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------------------------------------------

def _sex(ctx: Context) -> str | None:
    p = qc.out_path(ctx)
    if not p.exists():
        return None
    return (json.loads(p.read_text()).get("inferred_sex") or {}).get("call")


def _recount(ctx: Context) -> dict:
    import pysam

    from . import variants

    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    fasta_path = reference.reference_fasta(ctx)
    con = duckdb.connect()
    claims = con.sql(" UNION ALL BY NAME ".join(f"SELECT * FROM '{p}'" for p in claim_sources(ctx))).df()
    vp = variants.variants_parquet(ctx)
    pos_set = claims.dropna(subset=["chrom", "pos"])[["chrom", "pos"]].drop_duplicates()
    con.register("want", pos_set)
    vrows = con.sql(f"SELECT v.chrom, v.pos, v.ref, v.alt, v.gt, v.multiallelic FROM '{vp}' v "
                    "JOIN want w ON v.chrom = w.chrom AND v.pos = w.pos").df()
    sex = _sex(ctx)
    contigs = reference.cram_contigs(ctx)
    prefix = "chr" if "chr1" in contigs else ""

    def contig(ch: str) -> str:
        return f"{prefix}{'M' if ch == 'MT' and prefix else ch}"

    paths = pgx.out_paths(ctx)
    positions = pd.read_parquet(paths["positions"]) if paths["positions"].exists() else pd.DataFrame()
    genes = pd.read_parquet(paths["genes"]) if paths["genes"].exists() else pd.DataFrame()
    records = []
    with pysam.AlignmentFile(str(aln.path), "rc" if aln.kind == "cram" else "rb",
                             reference_filename=str(fasta_path)) as af, pysam.FastaFile(str(fasta_path)) as fa:
        targets = claim_targets(claims, vrows, sex)
        if len(positions):
            targets += pgx_targets(positions, fa, contig, sex)
        for t in targets:
            if t.chrom == "MT":
                c = Counts(note="mitochondrial DNA comes in hundreds of copies per cell; it isn't re-counted this way")
            else:
                c = count_site(af, fa, contig(t.chrom), t)
            records.append(site_record(t, c, classify(c, t)))
    sites = pd.DataFrame.from_records(records)
    sites["flags"] = sites["flags"].map(lambda f: list(f))
    sites["difficult"] = difficult_reasons(sites, Store(ctx.cfg).table("giab", "difficult"))
    sites.to_parquet(out_sites(ctx), index=False)
    verdicts = gene_verdicts(sites, genes)
    adjust_claims(claims, sites, verdicts).to_parquet(out_claims(ctx), index=False)
    if len(genes):
        adjust_genes(genes, verdicts).to_parquet(out_pgx_genes(ctx), index=False)
    summary = {
        "sites": len(sites),
        "by_source": {s: {k: int(v) for k, v in g.status.value_counts().items()} for s, g in sites.groupby("source")},
        "flags": {f: int(sum(f in fl for fl in sites["flags"])) for f in ("look_alike", "one_strand")},
        "genes": verdicts,
        "settings": SETTINGS,
    }
    write_json(out_summary(ctx), summary)
    return {"sites": len(sites), **{s: int((sites.status == s).sum()) for s in STATUSES}}


def _inputs(ctx: Context) -> list[Path]:
    from . import variants

    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    paths = pgx.out_paths(ctx)
    return [p for p in (aln.path if aln else None, *claim_sources(ctx), paths["positions"], paths["genes"],
                        variants.variants_parquet(ctx), qc.out_path(ctx), Store(ctx.cfg).table("giab", "difficult"))
            if p and p.exists()]


def _available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "no CRAM/BAM found (re-counting needs the reads)"
    if not reference.reference_fasta(ctx):
        return "CRAM reference not verified yet (reference stage)"
    if not claim_sources(ctx):
        return "no findings yet (annotate stage)"
    return None


STAGE = Stage(
    name="recount", version="2", title="Re-count the reads behind every finding", fn=_recount,
    inputs=_inputs, outputs=lambda ctx: [out_sites(ctx), out_claims(ctx), out_summary(ctx)],
    available=_available, params=lambda ctx: {"evidence_model": ev.MODEL_VERSION, "settings": SETTINGS},
)
