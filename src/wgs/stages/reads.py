"""Read-level QC that is cheap enough to run routinely.

* alignment_stats: `samtools stats` on chr20 (~2% of the genome) — representative mapping rate, duplicate rate,
  insert size, error rate, without decoding all 53 GB.
* fastq_stats: first 200k read pairs of each FASTQ — read length, base quality by cycle, GC, N rate, instrument.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from ..pipeline import Context, Stage, run, write_json
from .reference import reference_fasta

SAMPLE_READS = 200_000


def aln_out(ctx: Context) -> Path:
    return ctx.work / "qc" / "alignment_stats.json"


def fq_out(ctx: Context) -> Path:
    return ctx.work / "qc" / "fastq_stats.json"


def _alignment(ctx: Context) -> dict:
    cram = ctx.inv.first("cram") or ctx.inv.first("bam")
    import pysam

    with pysam.AlignmentFile(str(cram.path), "rc" if cram.kind == "cram" else "rb", check_sq=False) as af:
        region = next((r for r in ("chr20", "20") if r in af.references), af.references[0])
    args = ["samtools", "stats", "-@", "4"]
    if cram.kind == "cram":
        args += ["-r", str(reference_fasta(ctx)), "--reference", str(reference_fasta(ctx))]
    res = run([*args, str(cram.path), region], capture=True, log=ctx.logs / "reads.log")
    sn, is_hist, rl_hist = {}, [], []
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if parts[0] == "SN":
            sn[parts[1].rstrip(":")] = float(parts[2]) if re.match(r"^-?[\d.]+(e[-+]?\d+)?$", parts[2]) else parts[2]
        elif parts[0] == "IS":
            is_hist.append((int(parts[1]), int(parts[2])))
        elif parts[0] == "RL":
            rl_hist.append((int(parts[1]), int(parts[2])))
    total = sn.get("raw total sequences", 0) or 1
    out = {
        "region": region,
        "reads": int(total),
        "mapped_pct": round(100 * sn.get("reads mapped", 0) / total, 3),
        "properly_paired_pct": round(100 * sn.get("reads properly paired", 0) / total, 3),
        "duplicate_pct": round(100 * sn.get("reads duplicated", 0) / total, 3),
        "mapq0_pct": round(100 * sn.get("reads MQ0", 0) / total, 3),
        "error_rate": sn.get("error rate"),
        "insert_size_mean": sn.get("insert size average"),
        "insert_size_sd": sn.get("insert size standard deviation"),
        "average_quality": sn.get("average quality"),
        "average_length": sn.get("average length"),
        "insert_size_hist": [h for h in is_hist if h[0] <= 1000],
    }
    write_json(aln_out(ctx), out)
    return {k: v for k, v in out.items() if k != "insert_size_hist"}


def _fastq(ctx: Context) -> dict:
    out: dict = {"sampled_reads_per_file": SAMPLE_READS, "files": {}}
    for f in ctx.inv.of("fastq"):
        qsum: list[int] = []
        qn: list[int] = []
        lengths: dict[int, int] = {}
        gc = n_bases = total = 0
        q30 = 0
        ids = []
        with gzip.open(f.path, "rb") as fh:
            for i in range(SAMPLE_READS):
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().rstrip()
                fh.readline()
                qual = fh.readline().rstrip()
                if i < 3:
                    ids.append(header.decode().split()[0])
                lengths[len(seq)] = lengths.get(len(seq), 0) + 1
                gc += seq.count(b"G") + seq.count(b"C")
                n_bases += seq.count(b"N")
                total += len(seq)
                if len(qsum) < len(qual):
                    qsum.extend([0] * (len(qual) - len(qsum)))
                    qn.extend([0] * (len(qual) - len(qn)))
                for j, q in enumerate(qual):
                    qsum[j] += q - 33
                    qn[j] += 1
                    q30 += q - 33 >= 30
        m = re.match(r"@([A-Z]?\d+)L(\d)C(\d+)R(\d+)", ids[0]) if ids else None
        out["files"][f"R{f.meta.get('mate')}"] = {
            "file": f.path.name,
            "read_lengths": lengths,
            "gc_pct": round(100 * gc / max(total, 1), 2),
            "n_pct": round(100 * n_bases / max(total, 1), 4),
            "q30_pct": round(100 * q30 / max(total, 1), 2),
            "mean_quality_by_cycle": [round(s / n, 2) for s, n in zip(qsum, qn, strict=True)],
            "example_read_ids": ids,
            "flowcell": m.group(1) if m else None,
            "platform_guess": "MGI DNBSEQ" if m else None,
        }
    write_json(fq_out(ctx), out)
    return {k: {kk: vv for kk, vv in v.items() if kk not in ("mean_quality_by_cycle",)}
            for k, v in out["files"].items()}


def _aln_available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "no CRAM/BAM found"
    if ctx.inv.first("cram") and not reference_fasta(ctx):
        return "CRAM reference not verified yet (reference stage)"
    return None


ALIGNMENT = Stage(
    name="alignment_stats", version="3", title="Alignment stats (samtools stats, chr20)", fn=_alignment,
    inputs=lambda ctx: [f.path for f in (ctx.inv.of("cram") or ctx.inv.of("bam"))[:1]],
    outputs=lambda ctx: [aln_out(ctx)], available=_aln_available,
)
FASTQ = Stage(
    name="fastq_stats", version="1", title="FASTQ sample stats", fn=_fastq,
    inputs=lambda ctx: [f.path for f in ctx.inv.of("fastq")],
    outputs=lambda ctx: [fq_out(ctx)],
    available=lambda ctx: None if ctx.inv.of("fastq") else "no FASTQ found",
)
