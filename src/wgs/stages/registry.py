"""Ordered list of pipeline stages. Each stage skips itself (with a reason) when its inputs are absent."""

from __future__ import annotations

from . import coverage, genotypes, qc, reads, reference, variants

STAGES = [
    variants.STAGE,
    genotypes.STAGE,
    reference.STAGE,
    reads.FASTQ,
    reads.ALIGNMENT,
    coverage.STAGE,
    reference.FULL,
    qc.STAGE,
]
