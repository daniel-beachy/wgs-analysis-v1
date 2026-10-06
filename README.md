# wgs-analysis-v1

A personal whole-genome exploration tool. It turns a 30x WGS dataset into a fast, portable, **evidence-graded** dashboard, and it can be **refreshed** as genomic knowledge changes.

> **Privacy:** this repository contains code only. Genomic data, intermediate files and results live outside the repo (by default `../wgs-data/`). A pre-commit guard blocks genomic file types, large files, genotype-like rows and your sample ID.

> **Not medical advice:** results are for education and exploration. Confirm anything clinically relevant with a clinical-grade test and a genetic counsellor.

**New to genetics?** Start with the [friendly introduction to your genome](docs/GENETICS_PRIMER.md). It explains every term, file and report in this project, with diagrams.

## How it is laid out on the drive

```
Genomics/
├── WGS/                     your provider delivery (read-only; any layout works)
├── wgs-analysis-v1/         this repo (code)
├── wgs-data/                workspace: caches, intermediate files, dated result releases
├── Genome Dashboard.app     double-click on a Mac (Apple Silicon or Intel)
└── Genome Dashboard.exe     double-click on Windows (x64)
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
| `wgs run` | Runs or resumes every pipeline stage and publishes a new dashboard release |
| `wgs knowledge status` / `refresh` | Shows or updates the public knowledge versions (ClinVar, gnomAD, ClinGen, …) |
| `wgs query "<SQL>"` | Read-only SQL against the newest release (`claims`, `annotations`, `variants`, …); `-` reads SQL from stdin |

More commands are added as each development step lands (see `docs/PLAN.md`).

### Exploring with an AI assistant

Open an AI coding assistant (GitHub Copilot, Claude Code, …) in this repo and
ask about your results, e.g. *"Walk me through my carrier findings and how sure
we are about each."* [`AGENTS.md`](AGENTS.md) tells it where the data lives, how
to query it with `wgs query`, what every column means, and how to talk about
results responsibly (grades and reasons, not a diagnosis, ask before sensitive
topics, never upload raw files).

## Running the pipeline

```bash
pixi run wgs run                  # all stages; anything already up to date is skipped
pixi run wgs run --only qc        # a single stage
pixi run wgs run --force coverage # recompute one stage
```

Stages: `variants` (VCF → Parquet), `genotypes` (provider raw genotypes), `reference` (rebuild and verify the CRAM reference), `fastq_stats`, `alignment_stats`, `coverage` (mosdepth), `qc` and `publish`. Working files go to `wgs-data/work/<sample>/`. `publish` then writes a dated, read-only **release** to `wgs-data/releases/` — that is all the dashboard reads.

## Keeping the knowledge current

Interpretation uses public databases downloaded to `wgs-data/knowledge/`. Your variants are never sent online ([ADR-012](docs/DECISIONS.md)).

```bash
pixi run wgs knowledge status                 # which versions you have and when they were checked
pixi run wgs knowledge refresh                # download only sources that changed upstream (about 4 GB the first time)
pixi run wgs knowledge refresh --source clinvar --pin clinvar=20251006   # reproduce an archived ClinVar
pixi run wgs run                              # re-annotate and publish; a new release only if something changed
```

The dashboard's **What changed** tab then shows which sources moved to a new version, which findings were added, removed or regraded, and every variant whose ClinVar classification changed. A sensible routine is to refresh monthly. Every finding carries an evidence grade and a call-confidence grade, each with reasons ([ADR-013](docs/DECISIONS.md), and the "How this project grades its confidence" section of the [primer](docs/GENETICS_PRIMER.md)).

### Knowledge sources

| Source | Used for | Licence |
|---|---|---|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (variants + per-lab submissions) | Variant classifications, review stars, who says so | Public domain |
| [ClinGen](https://clinicalgenome.org/) | Gene–disease validity, dosage sensitivity, actionability | CC0 |
| [ACMG SF v3.3](https://www.gimjournal.org/) | 84 medically actionable genes and their reporting rules | Gene list from the published guideline |
| [Mondo](https://mondo.monarchinitiative.org/) | Disease names, definitions, cross-references | CC BY 4.0 |
| [HPO](https://hpo.jax.org/) | Inheritance per disease, gene–disease links | HPO licence (free, attribution) |
| [Orphanet](https://www.orphadata.com/) | Inheritance, onset, prevalence of rare diseases | CC BY 4.0 |
| [gnomAD](https://gnomad.broadinstitute.org/) v2.1.1 / v4.1 constraint | Population frequency; how well genes tolerate broken copies | CC0 |
| [1000 Genomes](https://www.internationalgenome.org/) | Population frequency fallback | Open (Fort Lauderdale) |
| [AlphaMissense](https://github.com/google-deepmind/alphamissense) | Missense damage prediction | CC BY-NC-SA 4.0 (personal, non-commercial use) |
| [REVEL](https://sites.google.com/site/revelgenomics/) | Missense damage prediction | Free for non-commercial use |

AlphaMissense and REVEL are licensed for non-commercial use, which this personal project is. The downloaded files stay in `wgs-data/knowledge/` and are never committed to the repo.

## Opening the dashboard (any machine, nothing to install)

Plug in the drive and double-click **`Genome Dashboard`** in the `Genomics` folder (`.app` on a Mac, `.exe` on Windows). Your browser opens on a local page (`http://127.0.0.1:8787`); the program quits by itself about a minute after you close the tab. Everything runs on your computer — no internet needed, nothing is uploaded.

- **First launch on a Mac:** if macOS says it can't verify the developer, right-click the app → **Open** → **Open** (only needed once per machine). Files copied by Finder from the internet get this flag; files built locally don't.
- **First launch on Windows:** if SmartScreen appears, click **More info → Run anyway**.
- The launcher finds `wgs-data/releases` by looking next to and above itself, so the drive letter or mount point doesn't matter. To point it elsewhere: `"Genome Dashboard.exe" --data D:\path\to\wgs-data` (or set `WGS_RELEASES`). A log is written to the temp folder as `genome-dashboard.log`.

**Rebuilding the dashboard** (after changing `viewer/` or `launcher/`, on the Mac): `pixi run build-dashboard`. It builds the viewer and both launchers in internal scratch space (exFAT can't hold npm's symlinks) and copies the results next to `wgs-data/`. For live development: `cd $TMPDIR/wgs-build/viewer && WGS_RELEASES=<path to wgs-data/releases> npx vite` after one build.

## About the reference genome

A CRAM file stores only how your reads differ from a reference genome. Each block also records an MD5 fingerprint of the reference used to compress it. tellmeGen did not ship that reference, so `wgs run --only reference` rebuilds it from Ensembl GRCh37 release 75. It checks the download against Ensembl's published checksum and renames contigs to match the CRAM, and N-masks the chrY PARs. It then proves the result:

- every CRAM block it decodes must match its stored fingerprint (on this dataset, 17.5M reads across all 84 contigs, 0 errors);
- every VCF REF allele must match (0 mismatches).

As a negative control, changing a single base made decoding fail. The optional one-time stage `wgs run --only reference_full` then decodes **every** read; on this dataset that was 868,564,724 reads with 0 mismatches (about 26 minutes). So the rebuilt reference matches the one tellmeGen used everywhere the CRAM depends on it. Details are in [ADR-008](docs/DECISIONS.md).

## Graceful degradation

Each dashboard section declares the inputs it **requires** and the inputs that **enhance** it. With only a raw genotype file you still get traits and ancestry. With a VCF you get everything except read-based checks. With CRAM plus a matching reference you also get coverage, callability and CYP2D6 structural calls. Outputs from optional external tools (PharmCAT, Cyrius, hap.py/vcfeval, pgsc_calc, Haplogrep, yhaplo, mosdepth, ...) are ingested when present but never required.

See [docs/DECISIONS.md](docs/DECISIONS.md) for why things are built this way.
