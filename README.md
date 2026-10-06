# wgs-analysis-v1

A personal whole-genome exploration tool. It turns a 30x WGS dataset into a fast, portable, **evidence-graded** dashboard, and it can be **refreshed** as genomic knowledge changes.

> **Privacy:** this repository contains code only. Genomic data, intermediate files and results live outside the repo (by default `../wgs-data/`). A pre-commit guard blocks genomic file types, large files, genotype-like rows and your sample ID.

> **Not medical advice:** results are for education and exploration. Confirm anything clinically relevant with a clinical-grade test and a genetic counsellor.

## How it is laid out on the drive

```
Genomics/
├── WGS/                 your provider delivery (read-only; any layout works)
├── wgs-analysis-v1/     this repo (code)
└── wgs-data/            workspace: caches, intermediate files, dated result releases, dashboard
```

Nothing is tied to fixed paths. Inputs are found by **content** (VCF/CRAM headers, magic bytes) anywhere under the search folders. Paths in `wgs.local.toml` are resolved relative to the file, so the drive can be mounted anywhere. You can override paths per machine (`[host."name"]`), with environment variables (`WGS_DATA`, `WGS_WORKSPACE`, `WGS_CONFIG`) or with CLI flags (`--data`, `--workspace`, `--sample`).

## Setup (analysis machine: macOS or Linux)

```sh
brew install pixi            # or: curl -fsSL https://pixi.sh/install.sh | sh
cd Genomics/wgs-analysis-v1
pixi install                 # creates the environment (stored on the host disk, not the exFAT drive)
pixi run install-hooks       # enable the data-leak guard
pixi run wgs init            # write wgs.local.toml with portable, relative paths
pixi run wgs doctor          # show what was found and which dashboard sections can run
```

## Commands

| Command | What it does |
|---|---|
| `wgs doctor` | Lists the inputs it found, which sections are ready, partial or unavailable (and why), tool versions and workspace health |
| `wgs init` | Writes `wgs.local.toml` |
| `wgs verify-inputs` | Checks inputs against `MANIFEST.sha256` (reads everything, so it is slow) |

More commands are added as each development step lands (see `docs/PLAN.md`).

## Graceful degradation

Each dashboard section declares the inputs it **requires** and the inputs that **enhance** it. With only a raw genotype file you still get traits and ancestry. With a VCF you get everything except read-based checks. With CRAM plus a matching reference you also get coverage, callability and CYP2D6 structural calls. Outputs from optional external tools (PharmCAT, Cyrius, hap.py/vcfeval, pgsc_calc, Haplogrep, yhaplo, mosdepth, ...) are ingested when present but never required.

See [docs/DECISIONS.md](docs/DECISIONS.md) for why things are built this way.
