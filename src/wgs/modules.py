"""Analysis modules and what they need.

Each dashboard section is produced by a module that declares the inputs it
*requires* (any one group must be satisfied) and the inputs that *enhance* it.
Missing enhancers degrade a section to "partial" with an explanation instead
of failing the run.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import Config
from .discovery import Inventory

VARIANTS = ("vcf", "gvcf")
GENOTYPES = ("vcf", "gvcf", "raw_genotypes")


@dataclass
class Enhancer:
    kinds: tuple[str, ...]  # all must be present
    benefit: str


@dataclass
class Module:
    id: str
    title: str
    requires_any: tuple[str, ...]
    enhancers: list[Enhancer] = field(default_factory=list)
    note: str = ""


MODULES: list[Module] = [
    Module("qc", "Data quality & QC", GENOTYPES, [
        Enhancer(("cram", "reference"), "genome-wide coverage, callable regions, sex/contamination checks"),
        Enhancer(("raw_genotypes",), "concordance with provider's genotype file"),
        Enhancer(("fastq",), "read-level metrics (base quality, read length, instrument)"),
    ]),
    Module("pgx", "Pharmacogenomics", VARIANTS, [
        Enhancer(("cram", "reference"), "CYP2D6 copy number / hybrid alleles from reads"),
        Enhancer(("external:pharmcat",), "ingest a PharmCAT report produced elsewhere"),
        Enhancer(("external:cyrius",), "ingest a Cyrius CYP2D6 call produced elsewhere"),
    ]),
    Module("health", "Health & disease risk", VARIANTS, [
        Enhancer(("cram", "reference"), "distinguish 'reference' from 'not covered' at key positions"),
    ], note="raw genotype files alone cover too few positions for monogenic screening"),
    Module("carrier", "Carrier status", VARIANTS, [
        Enhancer(("cram", "reference"), "distinguish 'reference' from 'not covered' at key positions"),
    ]),
    Module("traits", "Physical & behavioral traits", GENOTYPES, [
        Enhancer(("vcf",), "full polygenic scores rather than single-SNP traits"),
    ]),
    Module("prs", "Polygenic scores & percentiles", VARIANTS, [
        Enhancer(("gvcf",), "explicit reference calls instead of assuming reference when absent"),
        Enhancer(("cram", "reference"), "check coverage at score positions absent from the VCF"),
        Enhancer(("external:pgsc_calc",), "ingest pgsc_calc results produced elsewhere"),
    ]),
    Module("ancestry", "Ancestry & haplogroups", GENOTYPES, [
        Enhancer(("external:haplogrep",), "ingest Haplogrep mtDNA result"),
        Enhancer(("external:yhaplo",), "ingest yhaplo Y-chromosome result"),
    ]),
    Module("provider", "Comparison with provider reports", ("report",)),
    Module("benchmark", "Accuracy vs Genome in a Bottle", ("truth_vcf",), [
        Enhancer(("external:happy",), "ingest hap.py results"),
        Enhancer(("external:vcfeval",), "ingest RTG vcfeval results"),
    ], note="runs when the sample is a GIAB reference sample (e.g. HG002); truth sets are downloaded"),
]


def kinds_present(inv: Inventory, cfg: Config) -> set[str]:
    kinds = {f.kind for f in inv.files}
    cram = inv.first("cram")
    ref_dir = cfg.cache_dir / "reference"
    if cram and (ref_dir / "reference.ready").exists():
        kinds.add("reference")
    elif cram and "reference" in kinds and not any(f.meta.get("matches_cram") for f in inv.of("reference")):
        # A FASTA was found, but it only counts if it matches the CRAM's contigs.
        kinds.discard("reference")
    if inv.sample and inv.sample.upper().startswith(("HG00", "NA1", "NA2")) and \
            (cfg.cache_dir / "giab").exists():
        kinds.add("truth_vcf")
    return kinds


def evaluate(inv: Inventory, cfg: Config) -> list[dict]:
    present = kinds_present(inv, cfg)
    rows = []
    for m in MODULES:
        satisfied = [k for k in m.requires_any if k in present]
        if not satisfied:
            status, detail = "unavailable", "needs one of: " + ", ".join(m.requires_any)
        else:
            missing = [e for e in m.enhancers if not all(k in present for k in e.kinds)]
            got = [e for e in m.enhancers if all(k in present for k in e.kinds)]
            # External-tool outputs are pure bonuses; only missing primary data makes a section partial.
            status = "partial" if any(not k.startswith("external:") for e in missing for k in e.kinds) \
                else "ready"
            parts = []
            if got:
                parts.append("+ " + "; ".join(e.benefit for e in got))
            if missing:
                parts.append("- missing " + "; ".join(f"{'+'.join(e.kinds)} ({e.benefit})" for e in missing))
            detail = " | ".join(parts) or f"using {satisfied[0]}"
        rows.append({"id": m.id, "title": m.title, "status": status, "detail": detail, "note": m.note})
    return rows
