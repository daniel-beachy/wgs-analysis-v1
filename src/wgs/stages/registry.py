"""Ordered list of pipeline stages. Each stage skips itself (with a reason) when its inputs are absent."""

from __future__ import annotations

from . import annotate, coverage, genotypes, pgs, pgx, publish, qc, reads, reference, traits, variants

STAGES = [
    variants.STAGE,
    genotypes.STAGE,
    reference.STAGE,
    reads.FASTQ,
    reads.ALIGNMENT,
    coverage.STAGE,
    reference.FULL,
    qc.STAGE,
    annotate.STAGE,
    pgx.CYRIUS,
    pgx.HLA,
    pgx.PGX,
    traits.STAGE,
    pgs.STAGE,
    publish.STAGE,
]
