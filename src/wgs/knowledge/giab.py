"""Genome in a Bottle: the HG002 truth set, the genome's hard-to-sequence regions, and HG002 reads (ADR-023).

Genome in a Bottle (GIAB, run by the US National Institute of Standards and Technology) sequenced one person, HG002,
with many technologies and worked out the true letters at about 92% of the genome. Running a method on HG002's reads
and comparing with that truth shows how often the method is right, the standard way to measure accuracy.

What is downloaded:
- the HG002 truth set (v4.2.1, GRCh37): every true variant, plus the regions where the truth is reliable;
- GIAB's genome stratifications (v3.6): stretches of GRCh37 that short reads handle badly, each with its reason
  (repeated letters, near-identical copies elsewhere, extreme GC content, other known problem regions);
- HG002 reads sequenced on DNBSEQ (the same platform family as the owner's genome), streamed only for the ACMG genes'
  exons and the medicine-defining positions, and thinned to 30× so they match a typical consumer genome.

Wagner J, et al. Benchmarking challenging small variants with linked and long reads. Cell Genom. 2022;2:100128.
PMID 36452119. Dwarshuis N, et al. The GIAB genomic stratifications resource for human reference genomes.
Nat Commun. 2024;15:9029. PMID 39424793.
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import duckdb
import pandas as pd

from . import Built, Source, Store, Upstream
from .http import fetch, head

FTP = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/"
TRUTH = FTP + "release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark"
STRAT = FTP + "release/genome-stratifications/v3.6/GRCh37@all/"
READS = FTP + ("data/AshkenazimTrio/HG002_NA24385_son/NIST_BGIseq_2x150bp_100x/GRCh37/"
               "HG002_GRCh37_BGIseq-2x150-100x_NIST_20211126.bam")
# reason -> stratification file. GIAB's own "all difficult regions" union is these plus difficult X/Y regions (the
# truth set covers chromosomes 1-22 only), so the union of these files is used and every position keeps its reasons.
REASONS = {
    "repeat": "LowComplexity/GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
    "segdup": "SegmentalDuplications/GRCh37_segdups.bed.gz",
    "low_mappability": "Mappability/GRCh37_lowmappabilityall.bed.gz",
    "gc": "GCcontent/GRCh37_gclt25orgt65_slop50.bed.gz",
    "mhc": "OtherDifficult/GRCh37_MHC.bed.gz",
    "other": "OtherDifficult/GRCh37_allOtherDifficultregions.bed.gz",
}
PAD = 500            # bases kept around each exon and medicine position, so nearby non-coding variants count too
DEPTH = 30           # HG002 reads are thinned to this mean depth, to match a typical 30× consumer genome
SEED = 42            # thinning keeps whole read pairs, chosen by a hash of the read name with this seed
AUTOSOMES = {str(i) for i in range(1, 23)}


def merge(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Merge overlapping or touching intervals (0-based, half-open), sorted by chromosome then start."""
    out: list[list] = []
    for ch, s, e in sorted(intervals, key=lambda x: (x[0].zfill(2), x[1])):
        if out and out[-1][0] == ch and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([ch, s, e])
    return [tuple(x) for x in out]


def regions(workspace: Path) -> list[tuple[str, int, int]]:
    """ACMG SF exons (MANE Select) and PharmCAT defining positions, ±PAD, chromosomes 1-22 (the truth set's span)."""
    from ..stages.special import panel_exons

    st = Store(SimpleNamespace(workspace=workspace))
    iv = [(e["chrom"], int(e["start"]) - 1 - PAD, int(e["end"]) + PAD) for e in panel_exons(st)]
    pos = st.table("pharmcat", "positions")
    if pos:
        for ch, p, ref in duckdb.connect().execute(
                f"SELECT chrom37, pos37, ref FROM '{pos}' WHERE pos37 IS NOT NULL").fetchall():
            iv.append((ch.removeprefix("chr"), p - 1 - PAD, p - 1 + len(ref or "N") + PAD))
    return merge([(c, max(s, 0), e) for c, s, e in iv if c in AUTOSOMES])


def _fingerprint(iv: list[tuple[str, int, int]]) -> str:
    return hashlib.sha256(repr(iv).encode()).hexdigest()[:12]


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    out = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            ch, s, e = line.split("\t")[:3]
            out.append((ch.removeprefix("chr"), int(s), int(e)))
    return out


def _samtools() -> str:
    exe = shutil.which("samtools")
    if not exe:
        raise RuntimeError("samtools not found (run through pixi)")
    return exe


def stream_reads(iv: list[tuple[str, int, int]], scratch: Path, out: Path) -> dict:
    """Stream HG002 reads for `iv` from GIAB's server (only those regions are transferred), then thin to DEPTH."""
    sam = _samtools()
    bed = out / "regions.bed"
    bed.write_text("".join(f"chr{c}\t{s}\t{e}\n" for c, s, e in iv))
    full = scratch / "hg002_full.bam"
    # htslib downloads the remote index into the working directory, so run there
    subprocess.run([sam, "view", "-b", "-M", "-L", str(bed), "-o", str(full), READS], cwd=scratch, check=True)
    subprocess.run([sam, "index", str(full)], check=True)
    cov = subprocess.run([sam, "bedcov", str(bed), str(full)], capture_output=True, text=True, check=True).stdout
    total = length = 0
    for line in cov.splitlines():
        _, s, e, d = line.split("\t")[:4]
        total, length = total + int(d), length + int(e) - int(s)
    mean = total / length if length else 0
    if mean < DEPTH:
        raise RuntimeError(f"HG002 slice is only {mean:.0f}× deep; expected about 100×")
    frac = round(DEPTH / mean, 4)
    thin = out / "hg002.bam"
    subprocess.run([sam, "view", "-b", "--subsample", str(frac), "--subsample-seed", str(SEED), "-o", str(thin),
                    str(full)], check=True)
    subprocess.run([sam, "index", str(thin)], check=True)
    return {"regions": len(iv), "bases": length, "full_depth": round(mean, 1), "kept_fraction": frac,
            "depth": DEPTH}


class Giab(Source):
    id = "giab"
    title = "Genome in a Bottle HG002 (accuracy benchmark)"
    homepage = "https://www.nist.gov/programs-projects/genome-bottle"
    licence = "Public domain (NIST data)"
    cadence = "yearly"
    description = ("A reference person whose DNA letters are known almost exactly. Running the methods on their reads "
                   "shows how often they are right, and GIAB's map of hard-to-sequence regions says where calls are "
                   "least reliable.")
    schema = 1

    def probe(self, pin: str | None = None) -> list[Upstream]:
        ups = [head(TRUTH + ".vcf.gz"), head(TRUTH + ".vcf.gz.tbi"), head(TRUTH + "_noinconsistent.bed")]
        ups += [head(STRAT + f) for f in REASONS.values()]
        reads = head(READS)
        reads.etag = f"{_fingerprint(regions(self.workspace))}|{DEPTH}x|{PAD}"
        return ups + [reads]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        vcf, tbi, bed = upstream[:3]
        fetch(vcf.url, out / "truth.vcf.gz", vcf.size)
        fetch(tbi.url, out / "truth.vcf.gz.tbi", tbi.size)
        fetch(bed.url, out / "confident.bed", bed.size)
        rows = []
        for (reason, rel), up in zip(REASONS.items(), upstream[3:3 + len(REASONS)], strict=True):
            local = fetch(up.url, scratch / rel.rsplit("/", 1)[-1], up.size)
            rows += [(c, s, e, reason) for c, s, e in merge(read_bed(local)) if c in AUTOSOMES | {"X", "Y"}]
        con = duckdb.connect()
        d = pd.DataFrame(rows, columns=["chrom", "start", "end", "reason"])  # noqa: F841 (read by DuckDB)
        difficult = out / "difficult.parquet"
        con.execute(f"COPY (SELECT * FROM d ORDER BY chrom, start) TO '{difficult}' (FORMAT parquet)")
        confident = out / "confident.parquet"
        con.execute(f"""COPY (SELECT replace(chrom, 'chr', '') AS chrom, start, "end"
            FROM read_csv('{out / "confident.bed"}', delim='\t', header=false,
                          columns={{'chrom': 'VARCHAR', 'start': 'BIGINT', 'end': 'BIGINT'}}))
            TO '{confident}' (FORMAT parquet)""")
        iv = regions(self.workspace)
        if upstream[-1].etag and not upstream[-1].etag.startswith(_fingerprint(iv)):
            raise RuntimeError("the gene regions changed between checking and building; refresh again")
        if not iv:
            raise RuntimeError("no ACMG exons or PharmCAT positions found (refresh ensembl, mane and pharmcat first)")
        slice_ = stream_reads(iv, scratch, out)
        r = pd.DataFrame(iv, columns=["chrom", "start", "end"])  # noqa: F841 (read by DuckDB)
        reg = out / "regions.parquet"
        con.execute(f"COPY (SELECT * FROM r) TO '{reg}' (FORMAT parquet)")
        return Built(version="4.2.1", tables={"difficult": difficult, "confident": confident, "regions": reg},
                     upstream=upstream,
                     notes={"reads": slice_, "truth": "NISTv4.2.1 GRCh37, chromosomes 1-22",
                            "stratifications": "v3.6",
                            "citation": "Wagner J, et al. Benchmarking challenging small variants with linked and "
                                        "long reads. Cell Genom. 2022;2:100128. Dwarshuis N, et al. The GIAB genomic "
                                        "stratifications resource for human reference genomes. Nat Commun. "
                                        "2024;15:9029."})
