#!/usr/bin/env python3
"""Refuse to commit anything that could be personal genomic data.

Checks staged files (or every tracked file with --all) for:
  * genomic file types (VCF/BAM/CRAM/FASTQ/FASTA/Parquet/DuckDB/...)
  * files larger than 2 MB (except vendored WebAssembly modules under viewer/public/duckdb-ext/)
  * genotype-looking rows (rsID or chrom/pos followed by a genotype)
  * any string listed in [guard].forbidden_strings of wgs.local.toml (e.g. your kit ID)
Synthetic test fixtures under tests/fixtures may opt out with a first line containing
"wgs-guard: synthetic".
"""
from __future__ import annotations

import re
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BLOCKED = re.compile(r"\.(vcf|bcf|tbi|csi|cram|crai|bam|bai|fq|fastq|fa|fasta|fai|parquet|duckdb|gvcf)(\.b?gz)?$",
                     re.I)
GENOTYPE_ROW = re.compile(r"^(rs\d+|chr[\dXYM]+)[\t,;](\w+)?[\t,;]?\d{3,}[\t,;]([ACGT]{1,2}|[01][/|][01])\s*$",
                          re.M)
MAX_BYTES = 2 * 1024 * 1024
VENDORED_WASM = "viewer/public/duckdb-ext/"


def forbidden_strings() -> list[str]:
    cfg = ROOT / "wgs.local.toml"
    if not cfg.exists():
        return []
    data = tomllib.loads(cfg.read_text())
    return [s for s in data.get("guard", {}).get("forbidden_strings", []) if s]


def files(all_files: bool) -> list[str]:
    cmd = ["git", "ls-files"] if all_files else ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"]
    return [f for f in subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True).stdout.splitlines() if f]


def blob(path: str, all_files: bool) -> bytes:
    if all_files:
        return (ROOT / path).read_bytes()
    return subprocess.run(["git", "show", f":{path}"], cwd=ROOT, capture_output=True).stdout


def main() -> int:
    all_files = "--all" in sys.argv
    secrets = forbidden_strings()
    problems = []
    for f in files(all_files):
        if BLOCKED.search(f):
            problems.append(f"{f}: genomic file type")
            continue
        data = blob(f, all_files)
        if f.startswith(VENDORED_WASM) and f.endswith(".wasm") and data[:4] == b"\0asm":
            continue
        if len(data) > MAX_BYTES:
            problems.append(f"{f}: larger than 2 MB")
            continue
        text = data.decode("utf-8", errors="ignore")
        synthetic = "wgs-guard: synthetic" in text.splitlines()[0] if text else False
        if any(s in text for s in secrets):
            problems.append(f"{f}: contains a forbidden string from wgs.local.toml")
        if not synthetic and len(GENOTYPE_ROW.findall(text)) >= 5:
            problems.append(f"{f}: contains genotype-like rows")
    if problems:
        print("data_guard: commit blocked — possible personal genomic data:\n  " + "\n  ".join(problems))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
