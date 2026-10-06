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
