# Architecture decision record

Each entry gives the decision, the alternatives considered, and why. Newest at the bottom.

## ADR-001 Code, data and results live in separate sibling folders

- **Decision:** `Genomics/WGS` (inputs, treated as read-only), `Genomics/wgs-analysis-v1` (code, git), `Genomics/wgs-data` (workspace).
- **Why:** inputs are never modified, so the provider checksums stay valid. Code can be public without risking data leaks. The workspace can be deleted and rebuilt.
- **Inputs are read-only:** where a tool needs a file next to an input (for example a `.tbi` index), we write it to the workspace. htslib's `file.vcf.gz##idx##index.tbi` syntax points the tool at it.

## ADR-002 Inputs are discovered by content, not by path

- **Decision:** walk the configured search roots and classify files by magic bytes and headers (VCF `##contig` lengths identify the build; CRAM headers provide the sample, M5 tags and programs). Configuration is layered: CLI > environment > `wgs.local.toml` (with per-host sections) > defaults relative to the repo.
- **Alternatives:** a fixed folder layout (brittle across machines); absolute paths in a config file (breaks when the drive letter or mount point changes).
- **Why:** the same drive is mounted as `/Volumes/Daniel_SSD` on macOS and as `D:\` on Windows. The GIAB demo sample has a completely different layout.

## ADR-003 Reproducible environment with pixi (conda-forge + bioconda)

- **Decision:** a single `pixi.toml` and lock file provide Python, DuckDB, pysam, bcftools/samtools/htslib and OpenJDK (for PharmCAT, RTG vcfeval, snpEff and Haplogrep, which are Java and therefore cross-platform).
- **Alternatives:** Homebrew plus pip (not reproducible, and pins drift); Docker (heavy on macOS, and the user has no Linux machine); uv only (cannot provide bcftools or samtools).
- **Environment location:** the env is *detached* (`.pixi/config.toml`) and lives on the host disk. The exFAT drive has no hardlinks and uses 256 KB clusters, so a conda env with tens of thousands of small files would waste gigabytes and run slowly. Envs are architecture-specific anyway.

## ADR-004 No Linux machine required

- **Decision:** everything runs natively on Apple Silicon. Benchmarking uses RTG `vcfeval` (Java) rather than hap.py (Linux and Python 2 only). CYP2D6 uses Cyrius (pure Python). Polygenic scores are computed in DuckDB rather than with pgsc_calc (Nextflow). If outputs from those Linux tools ever exist, they are ingested.

## ADR-005 exFAT-aware storage

- macOS writes `._*` AppleDouble files on exFAT, so discovery, git and data loading ignore them.
- 256 KB clusters make many small files expensive, so results are stored as a few large columnar files rather than thousands of JSON shards.
- exFAT reports every file as executable, so the git `core.fileMode=false` setting and ruff EXE rules are disabled.

## ADR-006 Privacy model

- Local only for now: no hosting, no portfolio link. The GitHub repo is public and contains code only.
- `scripts/data_guard.py` runs as a pre-commit hook. It blocks genomic file types, files over 2 MB, genotype-like rows and strings listed in `[guard].forbidden_strings` (for example the kit ID).
- In the UI, serious findings (ACMG SF genes, APOE ε4) are shown behind a click-to-reveal panel with context.

## ADR-007 Columnar variant store (Parquet + DuckDB)

- **Decision:** the VCF is normalised once (`bcftools norm -m -any`) and written to zstd Parquet, sorted by chromosome and position, in 100k-row groups. All analyses and the viewer query it with DuckDB (native in Python, WASM in the browser).
- **Why:** a single portable file, about 85 MB for 7.3M records, that any OS or architecture can read. Sorted row groups let DuckDB prune by region, so a full-genome aggregate takes under a second. Parquet is an open standard that will still be readable in 10 years.
- **Detail:** when a multiallelic site is split, the allele not carried (e.g. allele 1 at a `0/2` site) produces a meaningless `0/0` PASS row. These rows are dropped. RefCall and NoCall rows are kept, so "tested and reference" can be told apart from "not called".
- **Alternatives:** keeping the bgzipped VCF with tabix (needs htslib, and is slow for aggregates); SQLite (row-oriented, about 10× larger, slow for scans); a DuckDB database file (one writer, and its format still changes between versions).

## ADR-008 Rebuilding the CRAM reference

- The provider CRAM has no `M5`/`UR` tags, and tellmeGen did not supply its reference (`GRCh37.primary_assembly.par_y_masked`). Its contig names and lengths match Ensembl GRCh37 release 75 `primary_assembly` with a `chr` prefix, `chrM` (rCRS, 16569) and `GL*.1` scaffolds.
- **Decision:** download Ensembl r75 and check it against Ensembl's `CHECKSUMS`. Rename the contigs and N-mask the chrY PARs. Write the contigs in CRAM header order.
- **Verification (both must pass before anything uses the reference):**
  1. Decode reads from 3 random 200 kb windows on every main contig and from every scaffold. htslib verifies each CRAM slice's reference MD5 while reconstructing sequence.
  2. `bcftools norm --check-ref` must find 0 VCF REF mismatches.
- **Full proof (optional, one-time):** the `reference_full` stage decodes every read (`samtools view -u | samtools view -c`) and writes `full_decode.json`. On this dataset it decoded 868,564,724 reads with 0 MD5 mismatches.
- If verification fails, the CRAM-dependent stages (coverage, alignment stats) are skipped. Everything else still works from the VCF.

## ADR-009 Re-runnable stages with stamps

- `wgs run` runs a fixed list of stages. Each one writes a stamp: stage version, input size and mtime, and parameters. Unchanged stages are skipped, `--force` reruns them, and failures write tracebacks to `logs/<stage>.log`.
- Large temporary files go to internal-disk scratch (`$WGS_SCRATCH` or `$TMPDIR`), not the exFAT drive.
- Coverage uses `mosdepth --fast-mode` with quantized bins (0, 1–4, 5–9, 10–149, 150+). The VCF is not a gVCF, so the callable BED is what separates "homozygous reference" from "not covered" at any position.

## ADR-010 Learning first

- **Principle:** the tool is for someone who is curious and learning about their own genome. That shapes every section and every UI decision.
- **In practice:** every metric, term and finding gets a plain-language explanation that you can reach in one hover or click: tooltips, "what does this mean?" expanders, small diagrams, and links into `docs/GENETICS_PRIMER.md`. Jargon is always defined where it first appears. Each claim shows its evidence and confidence in words, not just as a score. Reports explain *why* something matters, or why it probably doesn't.
