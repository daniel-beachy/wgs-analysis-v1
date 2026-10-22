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
| `wgs refresh [--pgs]` | Updates the knowledge, re-runs what it affects and publishes a release with a "What changed" report; polygenic-score sources only with `--pgs` (hours) |
| `wgs audit [--only curated]` | Re-checks the newest release against live ClinVar, GWAS Catalog, PGS Catalog, PubMed, source links and every cited baseline quote (needs internet; summary shown on the Data quality tab) |
| `wgs query "<SQL>"` | Read-only SQL against the newest release (`claims`, `annotations`, `variants`, …); `-` reads SQL from stdin |
| `wgs mcp` | The same read-only access for an AI assistant that speaks MCP (Claude Code, Claude Desktop, …); the assistant starts it, see below |

More commands are added as each development step lands (see `docs/PLAN.md`).

### Exploring with an AI assistant

Open an AI coding assistant (GitHub Copilot, Claude Code, …) in this repo and
ask about your results, e.g. *"Walk me through my carrier findings and how sure
we are about each."* [`AGENTS.md`](AGENTS.md) tells it where the data lives, how
to query it with `wgs query`, what every column means, and how to talk about
results responsibly (grades and reasons, not a diagnosis, ask before sensitive
topics, never upload raw files).

Assistants that speak the Model Context Protocol can read the results without a
shell through `wgs mcp`, a read-only server over the newest release with five
tools: `guide` (AGENTS.md), `overview`, `highlights`, `result` (one row by table
and key) and `query` (one SELECT, at most 1,000 rows). It uses the same sandbox
as `wgs query`: it can read only the releases folder, never the FASTQ, CRAM or
VCF, and can't write anything ([ADR-027](docs/DECISIONS.md)). Claude Code
offers it automatically in this repo (`.mcp.json`; approve it once). Elsewhere:

```bash
claude mcp add wgs -- pixi run --manifest-path /path/to/wgs-analysis-v1/pixi.toml -q wgs mcp
```

For Claude Desktop, add the same command under `mcpServers` in its
`claude_desktop_config.json`. Whatever the assistant reads goes to wherever its
model runs, so with a cloud model, ask narrow questions rather than for whole
tables.

In the dashboard, every highlight and every finding's details have an **✨ Ask
AI** button. It copies that one result, its grades with their reasons, its
sources and the rules above as a prompt you can paste into any assistant.
Nothing is sent until you paste it.

## Running the pipeline

```bash
pixi run wgs run                  # all stages; anything already up to date is skipped
pixi run wgs run --only qc        # a single stage
pixi run wgs run --force coverage # recompute one stage
```

Stages: `variants` (VCF → Parquet), `genotypes` (provider raw genotypes), `reference` (rebuild and verify the CRAM reference), `fastq_stats`, `alignment_stats`, `coverage` (mosdepth), `qc`, the interpretation stages (`annotate`, `claims`, `cyp2d6`, `hla`, `pgx`, `traits`, `pgs`) and `publish`. Working files go to `wgs-data/work/<sample>/`. `publish` then writes a dated, read-only **release** to `wgs-data/releases/` — that is all the dashboard reads.

## Keeping the knowledge current

Interpretation uses public databases downloaded to `wgs-data/knowledge/`. Your variants are never sent online ([ADR-012](docs/DECISIONS.md)).

```bash
pixi run wgs refresh                          # the routine: update the knowledge, re-run, publish (about 20 minutes)
pixi run wgs knowledge status                 # which versions you have and when they were checked
pixi run wgs knowledge refresh                # only download sources that changed upstream (about 4 GB the first time)
pixi run wgs knowledge refresh --source clinvar --pin clinvar=20251006   # reproduce an archived ClinVar
pixi run wgs run                              # re-annotate and publish; a new release only if something changed
pixi run wgs audit                            # confirm every cited record and quote still says what we show
```

`wgs refresh` leaves the polygenic-score sources (PGS Catalog and its 1000 Genomes reference scoring) alone unless you add `--pgs`, because re-scoring the reference panel takes hours. The dashboard's **What changed** tab then shows which sources moved to a new version; which findings were added, removed or regraded, each with the source updates behind it ("Because: ClinVar 2026-10-04 → 2026-10-19"); other results that moved (medicine genes, polygenic percentiles, ancestry shares, haplogroups, hard-to-read genes, repeats); and every variant whose ClinVar classification changed ([ADR-024](docs/DECISIONS.md)). A sensible routine is to refresh monthly. Every finding carries an evidence grade and a call-confidence grade, each with reasons ([ADR-013](docs/DECISIONS.md), and the "How this project grades its confidence" section of the [primer](docs/GENETICS_PRIMER.md)).

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
| [PharmCAT](https://pharmcat.clinpgx.org/) 3.4 (CPIC, DPWG and FDA guidance via [ClinPGx](https://www.clinpgx.org/)) | Star-allele calling and prescribing guidelines | MPL-2.0 tool; ClinPGx data CC BY-SA 4.0 |
| [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) | HLA allele sequences for T1K | CC BY-ND 4.0 |
| [MedlinePlus Genetics](https://medlineplus.gov/genetics/) | Plain-language gene descriptions (about 1,500 genes) | Public domain (NLM) |
| [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene) | Official gene names and summaries for genes MedlinePlus doesn't cover | Public domain (NCBI) |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | Single-variant trait associations (featured list) | EMBL-EBI open data |
| [PGS Catalog](https://www.pgscatalog.org/) | Polygenic scores, their publications and independent evaluations | EMBL-EBI open data; per-score licence recorded |
| [1000 Genomes reference panel](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/) (PGS Catalog build) | Ancestry-matched percentiles for polygenic scores (≈ 7 GB download) | Open access |
| [SFARI Gene](https://gene.sfari.org/) | Genes linked to autism and neurodevelopment | Free for research use |
| [STRchive](https://strchive.org) | Repeat-expansion disease loci: typical and disease ranges, evidence, references | MIT |
| [gnomAD-SV v2.1](https://gnomad.broadinstitute.org/) | How common each deletion/duplication is (GRCh37, ≈10,800 genomes) | CC0 |
| [IGSR populations](https://www.internationalgenome.org/data-portal/population) | Names and sampling locations of the 26 1000 Genomes populations (Ancestry tab) | Free to use without restriction |
| [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) HG002 v4.2.1 + stratifications v3.6 | Accuracy benchmark: a reference person's true variants, their DNBSEQ reads for the panel regions (≈ 55 MB after thinning to 30×), and the map of hard-to-sequence regions | Public domain (NIST data) |

Every release also runs **self-checks**: consistency rules over its own tables (every finding has a source and grade, every matched guideline has text, buckets follow the curated flags, and so on). Results are on the Data quality tab and in `checks.json`; a failure is shown loudly rather than blocking the release. The rules live in `src/wgs/checks.py`. See [ADR-016](docs/DECISIONS.md).

**Risk in absolute numbers.** Where it can be traced to sources, a polygenic score also shows "≈ X in 100 people with your score vs Y in 100 typical". The typical risk and any per-SD effect come from `src/wgs/knowledge/data/baseline_risk.toml`, a hand-curated table where every row carries its source URL, PubMed id and the exact quoted sentence; `wgs audit` confirms each quote is still on its page. Otherwise the card says "not established" and why ([ADR-018](docs/DECISIONS.md)).

AlphaMissense and REVEL are licensed for non-commercial use, which this personal project is. The downloaded files stay in `wgs-data/knowledge/` and are never committed to the repo.

## Opening the dashboard (any machine, nothing to install)

Plug in the drive and double-click **`Genome Dashboard`** in the `Genomics` folder (`.app` on a Mac; the Windows `.exe` is built too but hasn't been tested on a Windows machine). Your browser opens on a local page (`http://127.0.0.1:8787`); the program quits by itself about a minute after you close the tab. Everything runs on your computer — no internet needed, nothing is uploaded.

- **First launch on a Mac:** if macOS says it can't verify the developer, right-click the app → **Open** → **Open** (only needed once per machine). Files copied by Finder from the internet get this flag; files built locally don't.
- **First launch on Windows:** if SmartScreen appears, click **More info → Run anyway**.
- The launcher finds `wgs-data/releases` by looking next to and above itself, so the drive letter or mount point doesn't matter. To point it elsewhere: `"Genome Dashboard.exe" --data D:\path\to\wgs-data` (or set `WGS_RELEASES`). A log is written to the temp folder as `genome-dashboard.log`.

**Rebuilding the dashboard** (after changing `viewer/` or `launcher/`, on the Mac): `pixi run build-dashboard`. It builds the viewer and both launchers in internal scratch space (exFAT can't hold npm's symlinks) and copies the results next to `wgs-data/`. For live development: `cd $TMPDIR/wgs-build/viewer && WGS_RELEASES=<path to wgs-data/releases> npx vite` after one build.

## Demo release (shareable, no personal data)

`config/demo.toml` builds a second workspace, `wgs-data-demo/`, from Genome in a Bottle's HG002 truth VCF (v4.2.1, GRCh37) instead of your genome. HG002 is a reference genome shared openly through the Personal Genome Project, so the demo can be shown to anyone. Every page carries a banner saying it shows HG002, not the dashboard owner. It borrows the knowledge store (read-only) and the rebuilt reference from `wgs-data/`, so it needs a normal run and the `giab` knowledge source first.

```bash
pixi run demo                                                      # prepare the input, then run every stage
pixi run wgs --config config/demo.toml demo --only annotate --force annotate   # one stage
pixi run wgs --config config/demo.toml query                       # query the demo release
```

`wgs demo` keeps only the truth set's confident regions, renames contigs to `chr…` and drops every per-sample field except the genotype and its quality: the truth set's depths are summed over many sequencing runs and would read as one run's reads. With no reads, the callable map is the confident regions, and call confidence follows a truth-set rule instead of read depth (ADR-023 addendum). The truth set covers chromosomes 1–22 only, so there is no sex inference and no haplogroups. Stages that need reads or the provider's genotype file skip with a reason: `genotypes`, `fastq_stats`, `alignment_stats`, `reference_full`, `cyp2d6`, `hla`, `special`, `recount` and `secondcall`. Polygenic scores (and the ancestry estimate built on their PCA) are left out by default because the 1000 Genomes panel has to be scored for this genome's positions, which takes hours; add `--pgs` to include them.

Open it with the same launcher, pointed at the demo folder. It picks a different port when your own dashboard is already open, so the two never mix:

```bash
open -a "../Genome Dashboard.app" --args -data "$PWD/../wgs-data-demo"
```

On Windows: `"Genome Dashboard.exe" -data wgs-data-demo` from the `Genomics` folder.

## About the reference genome

A CRAM file stores only how your reads differ from a reference genome. Each block also records an MD5 fingerprint of the reference used to compress it. tellmeGen did not ship that reference, so `wgs run --only reference` rebuilds it from Ensembl GRCh37 release 75. It checks the download against Ensembl's published checksum and renames contigs to match the CRAM, and N-masks the chrY PARs. It then proves the result:

- every CRAM block it decodes must match its stored fingerprint (on this dataset, 17.5M reads across all 84 contigs, 0 errors);
- every VCF REF allele must match (0 mismatches).

As a negative control, changing a single base made decoding fail. The optional one-time stage `wgs run --only reference_full` then decodes **every** read; on this dataset that was 868,564,724 reads with 0 mismatches (about 26 minutes). So the rebuilt reference matches the one tellmeGen used everywhere the CRAM depends on it. Details are in [ADR-008](docs/DECISIONS.md).

## Graceful degradation

Each dashboard section declares the inputs it **requires** and the inputs that **enhance** it. With only a raw genotype file you still get traits and ancestry. With a VCF you get everything except read-based checks. With CRAM plus a matching reference you also get coverage, callability and CYP2D6 structural calls. Outputs from optional external tools (PharmCAT, Cyrius, hap.py/vcfeval, pgsc_calc, Haplogrep, yhaplo, mosdepth, ...) are ingested when present but never required.

### Medicines (pharmacogenomics)

`wgs run` includes three stages: `cyp2d6` (Cyrius, about 20 minutes, from the CRAM), `hla` (T1K, a few minutes) and `pgx` (PharmCAT, about 2 minutes). Cyrius needs its own Python environment, installed once with `pixi install -e cyrius`. T1K is installed for Apple Silicon only. On other machines, run [Cyrius](https://github.com/Illumina/Cyrius) or [T1K](https://github.com/mourisl/T1K) yourself and put the outputs (a Cyrius `.tsv` with "cyrius" or "cyp2d6" in its name; T1K's `*_genotype.tsv`) in the WGS folder; they are picked up automatically. Without either, PharmCAT still runs and reports *CYP2D6*/HLA as not typed. Design: [ADR-015](docs/DECISIONS.md).

### Traits, polygenic scores and Brain & mind

`traits` (seconds) reads featured single-variant traits and rare variants in SFARI genes. `pgs` computes polygenic score percentiles the way the PGS Catalog's pgsc_calc does, with `pgscatalog-utils` and `plink2` from a separate environment (`pixi install -e pgs`, macOS and Linux). It needs the `pgs_catalog` and `pgs_reference` knowledge sources. The first run takes 2–3 hours on a 16 GB laptop (reference QC, score normalising and reference-panel scores are cached afterwards, so re-runs take minutes) and needs about 10 GB of temporary space on the internal disk. Without the `pgs` environment the Traits page still shows the single-variant traits. Design: [ADR-017](docs/DECISIONS.md).

### Hard-to-read genes, repeat expansions and HLA

`special` (about 5 minutes) re-reads what a VCF can't show from the CRAM: *SMN1*/*SMN2* copy number (SMNCopyNumberCaller), *GBA1* (Gauchian), repeat-expansion disease loci (ExpansionHunter against STRchive's ranges), HLA disease associations (T1K types against a quoted, curated table), blood group (ABO, RhD) and deletions/duplications in the 84 ACMG genes (Delly plus GC-corrected exon depth, compared with gnomAD-SV). It also writes the "how your genome compares" counts. The tools live in their own environment: `pixi install -e special`. Illumina's callers are pinned to their last Apache-2.0 releases (see ADR-021). Needs the `strchive` knowledge source. Design: [ADR-021](docs/DECISIONS.md).

### Ancestry and haplogroups

`recount` (under a minute) goes back to the raw reads in the CRAM at every position a finding or a medicine result rests on and counts which version each read shows, separately from the variant caller. Where the reads fit a different genotype the finding's call drops to Low; weak, unclear or look-alike-region support caps it at Medium. Design: [ADR-022](docs/DECISIONS.md).

`benchmark` (about 3 minutes) measures the methods on Genome in a Bottle's HG002, whose true variants are known: the re-count is given HG002's true genotypes and deliberately wrong ones, and bcftools makes its own calls from HG002's reads, scored with RTG Tools `vcfeval` (BSD-2-Clause). Every result is split into ordinary and hard-to-sequence stretches. `secondcall` (about 3 minutes) runs the same bcftools command on your reads in the ACMG genes and medicine positions and lists what it finds that your VCF lacks, each re-counted and looked up in ClinVar, gnomAD and 1000 Genomes. Neither changes a finding; findings in hard-to-sequence stretches get a call reason that says so. Needs the `giab` knowledge source (`wgs knowledge refresh --source giab`: about 0.6 GB downloaded, 300 MB kept; RTG's 1 GB copy of the reference goes to temporary storage). Design: [ADR-023](docs/DECISIONS.md).

`ancestry` (under a minute) reuses the polygenic-score step's ancestry PCA to write your DNA as the closest mix of the five 1000 Genomes continental groups and to find your nearest reference population. Both are calibrated by scoring every reference person with themselves left out, so the page shows each group's noise floor and which populations the method can't tell apart. It also reads your maternal lineage (mtDNA, haplogrep3 + PhyloTree) from the VCF and, for XY samples, your paternal lineage (Y chromosome, Yleaf + YFull tree) from the Y reads in the CRAM, or from the VCF when there is no alignment. The haplogroup tools live in their own environment: `pixi install -e haplo` (Yleaf is installed from its pinned GitHub release because the bioconda build leaves out its marker data). Design: [ADR-020](docs/DECISIONS.md).

See [docs/DECISIONS.md](docs/DECISIONS.md) for why things are built this way.
