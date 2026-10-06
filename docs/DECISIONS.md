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

## ADR-011 Releases, an in-browser query engine and a tiny launcher

- **Releases:** the `publish` stage turns working files into an immutable release: `wgs-data/releases/<YYYY-MM-DD_HHMMSS>/manifest.json` plus content-keyed files in `releases/objects/` (`<name>-<hash>.parquet|json`, each with its sha256). Keys hash the table name, transform version and input fingerprint, so unchanged tables are shared between releases, not copied. A new release is written only when something changed, and each manifest links to its `previous` one. This is the basis for "What changed" (Steps 3 and 10): two releases can be diffed table by table. `index.json` lists releases newest first.
- **Viewer:** a static Svelte 5 + TypeScript app. Parquet tables are queried in the browser with **DuckDB-WASM** over HTTP range requests, so even the 7.3M-row variant table is never loaded whole. Page loads take about 1–2 s, filtered queries about 0.1–0.8 s. JSON documents (QC, coverage…) feed the charts. Charts are hand-written SVG components (no chart library). The genetics primer is bundled and its glossary drives every tooltip (ADR-010). Mermaid loads only on the Learn tab. Two details matter for this. First, this DuckDB-WASM build defaults to downloading whole files, so `db.ts` opens the database with `forceFullHTTPReads: false` and `reliableHeadRequests: true`. Without these, every page load fetched about 350 MB; with them, Health opens in under 1 s. Second, the Parquet reader extension is vendored in `viewer/public/duckdb-ext/` and loaded through `custom_extension_repository`, so the dashboard works offline and does not depend on the DuckDB CDN staying the same for a decade. When upgrading `@duckdb/duckdb-wasm`, download the matching `v<duckdb>/wasm_eh/parquet.duckdb_extension.wasm`.
  - *Why:* the same build runs on any OS and CPU with no install. Results stay as open formats (Parquet/JSON) that remain queryable outside the app, and the SQL tab exposes them directly.
  - *Alternatives:* Streamlit/Dash/Shiny (need a Python or R runtime on every laptop; not architecture-neutral); Electron/Tauri (heavy, or platform-specific toolchains); pre-rendered static HTML (no search over 7M rows); Datasette (Python runtime).
- **Launcher:** a dependency-free Go program, cross-compiled from the Mac for Windows x64 (`-H windowsgui`, icon via go-winres) and macOS (universal arm64 + x86_64 `.app`, ad-hoc signed). It embeds the built viewer, finds `wgs-data/releases` by walking up from its own location (or `--data` / `WGS_RELEASES`), serves both on `127.0.0.1` with Range support, opens the default browser, and refuses dotfiles and AppleDouble `._*` files. Port 8787 is preferred; if another instance already answers, it just opens a browser tab. The page sends a heartbeat every 10 s and the launcher exits about 75 s after the last one.
  - *Why a local server at all:* browsers block `fetch`/WebAssembly workers from `file://` pages, so double-clicking an `index.html` cannot work reliably. A single static binary is the simplest thing that works on both laptops with nothing installed.
  - DuckDB-WASM's "eh" bundle is used, which needs no cross-origin isolation headers.
- **Builds** run on the analysis Mac (`pixi run build-dashboard`) in internal scratch space, because exFAT has no symlinks for `node_modules/.bin`. Node and Go come from the pixi environment.

## ADR-012 Local knowledge store, refreshed by version

- **Privacy rule:** knowledge is downloaded and joined to your variants **on this computer**. Variants are never sent to web APIs (no VEP REST, no MyVariant.info). Links out (gnomAD, ClinVar, dbSNP) are only opened if you click them.
- **Sources (all free; versions are recorded per release):**

  | Source | Used for | Size | Update cadence |
  |---|---|---|---|
  | ClinVar (GRCh37 VCF) | Classifications, review stars, conditions | ~180 MB | Weekly |
  | Ensembl GFF3 GRCh37.87 | Gene models for `bcftools csq`; gene search | ~50 MB | Frozen (GRCh37) |
  | ClinGen gene validity + dosage | How strong each gene–disease link is; inheritance | ~2 MB | Daily export, rebuilt only on content change |
  | 1000 Genomes phase 3 sites (v5b, NCBI mirror — EBI is ~30× slower) | Overall and continental allele frequencies | ~2 GB | Frozen |
  | gnomAD v2.1.1 (slivar gnotate) | Popmax AF across about 140k people | ~2.2 GB | Frozen for GRCh37 |

- **Store:** `wgs-data/knowledge/<source>/<version>/` holds Parquet tables. `knowledge/index.json` records the upstream URL, Last-Modified/ETag/size fingerprint, build time and table checksums of each version, plus which version is *current*. A refresh `HEAD`s each source and rebuilds only those whose fingerprint changed. It builds into `.building/` and renames on success, so a failed download never corrupts the current version. The last 4 versions are kept, and the current one is never pruned. `--pin clinvar=YYYYMMDD` reproduces an archived ClinVar release, which makes past results reproducible.
- **Consequences:** `bcftools csq` (haplotype-aware, about 1 minute genome-wide) with Ensembl GFF3. *Rejected:* VEP with a 24 GB cache (too heavy for a portable drive) and snpEff (needs Java plus a 1 GB database, and its HGVS doesn't add much here). echtvar would be ideal for gnomAD, but it has no GRCh37 archive and no osx-arm64 build. slivar's gnotate zip provides the same lookup.
- **Release diff:** `publish` records the knowledge versions and evidence-model version in each manifest, and computes `changes.json` against the previous release. It lists sources whose version changed, claims added, removed or changed (field by field), and every ClinVar reclassification among your alleles. The diff is computed *after* the release digest, so it can never trigger a release on its own.
- **Why:** this fulfils "not a snapshot". Each finding traces back to a dated source version, and each change has a visible cause.

## ADR-013 Evidence model: two grades, weaker wins

- **Problem:** a variant "meaning something" depends on two independent questions. How strong is the published evidence for its effect? And is your genotype really there? A strong ClinVar classification on a low-quality call is not a strong finding, and neither is a perfect call of a 0★ assertion.
- **Decision:** every claim (`claims` table) carries an **evidence grade** (Strong / Moderate / Limited), a **call confidence** (High / Medium / Low / Not callable) and an **overall** grade = the weaker of the two. Each grade has a list of plain-language **reasons** and a list of **sources** with versions. Rules live in `src/wgs/evidence.py` (`MODEL_VERSION`). Changing a rule bumps the version, and the release diff shows its effect like a data update.
- **Evidence grade (v1, ClinVar-based claims):**
  - Review stars: 3–4★ → Strong; 2★ → Moderate; 0–1★ → Limited.
  - *Drug response* caps at Moderate until checked against CPIC/DPWG (Step 5). *Conflicting* → Limited.
  - Pathogenic or likely pathogenic, with population allele frequency (gnomAD popmax, otherwise 1000 Genomes) ≥ 5% → Limited (ACMG BA1). ≥ 1% → down one step (BS1-like).
  - Gene–disease validity in ClinGen *Disputed/Refuted* → Limited. *Limited* → down one step.
- **Call confidence:**
  - Non-PASS → Low.
  - GQ ≥ 30 and DP ≥ 15 → High. GQ ≥ 20 and DP ≥ 10 → Medium. Otherwise Low.
  - Unexpected read balance (het outside 20–80%, hom < 80%) → down one step.
  - Disagreement with the provider's genotype file (SNVs) → Low.
- **Which variants become claims:** ClinVar pathogenic, likely pathogenic, risk factor, association, protective and drug response; conflicting only when a submitter said pathogenic. VUS and benign are never claims; they stay visible in the variant explorer.
- **Routing to sections:**
  - Drug response → Medicines.
  - Association → Traits.
  - Pathogenic, likely pathogenic or conflicting in one copy of a gene whose ClinGen inheritance is recessive-only → Carrier.
  - Everything else → Health.
  - Steps 4–6 add curated claims (ACMG SF genes, CPIC diplotypes, GWAS/PGS traits) with their own evidence rules into the same table.
- **Sensitive findings:** pathogenic or likely pathogenic, and anything in genes such as *BRCA1/2*, *APOE*, *HTT* and *PRNP*, are masked until you click to reveal (per browser session). Each comes with a note that most such findings in healthy people are carrier status, low penetrance or misclassification, and need clinical confirmation.
- **Alternatives considered:** a single numeric score (hides *why*); full ACMG/AMP automated classification such as InterVar (it over-calls without the lab evidence, and is still a black box to a learner); copying ClinVar labels verbatim (ignores call quality and population frequency, both of which caught real problems in testing).

## ADR-014 Health and carrier: inheritance, per-lab conditions, actionability and predictions (evidence model v2)

- **Problem:** with v1, a ClinVar "pathogenic" label went straight into a claim. Three things made that misleading for one person's genome. (1) ClinVar's merged record names every condition any lab mentioned, and some labs list every disease of the gene, so *GJB2* W77R showed dominant skin syndromes next to the recessive hearing loss it is actually known for. (2) Whether one copy matters depends on how the *specific* condition is inherited, not on the gene. (3) Rare variants no lab has classified were invisible, even when every predictor agrees they are damaging.
- **Decision:**
  - **Curated knowledge added** (all downloaded and joined locally): Mondo (disease names and cross-references), HPO and Orphanet (inheritance per disease, gene–disease links), ClinGen dosage and actionability, the ACMG SF v3.3 gene list with each gene's reporting rule (`src/wgs/knowledge/data/acmg_sf_v3.3.csv`), and ClinVar's per-lab `submission_summary` (about 6.8M submissions).
  - **Which condition is meant:** for each variant, each lab submission votes for the conditions it names. Submissions naming one or two conditions count; gene-wide lists do not. Conditions with at least 20% of the top condition's support are kept; the rest are shown as "also listed". The "Who says so" table shows every lab, its verdict and date.
  - **Inheritance:** per condition from HPO/Orphanet via Mondo cross-references, falling back to the gene's modes. Combined with zygosity and sex chromosome this gives a **role**: *carrier*, *affected* (both copies, or the only X copy), *possible* (dominant), or *unknown*. Mixed genes (e.g. *SDHA*: recessive complex II deficiency and dominant paraganglioma) get a per-condition role and a statement that explains both. Only "disease-like" claims (P/LP, conflicting, predicted) get a role; risk factors and associations don't.
  - **Routing:** disease-like and carrier role → Carrier; otherwise Health groups *monogenic*, *compound* (two P/LP variants in one recessive gene — flagged as "possible", since we can't tell if they're on different copies), *predicted*, *uncertain* (conflicting) and *risk* (risk factor, protective, low penetrance). Blood-group "phenotypes" go to Traits.
  - **ACMG SF:** a P/LP claim in a listed gene is *reportable* only if it meets the gene's rule (biallelic for recessive genes, *HFE* only homozygous C282Y, *TTN* only truncating). The `acmg_genes` release table also records the callable fraction of each gene, so "no findings" comes with how much was actually readable.
  - **Computer predictions** (AlphaMissense, REVEL, gnomAD v4.1 constraint) apply only to rare (< 0.1% in every population), PASS variants in known disease genes that ClinVar hasn't classified (or calls benign). Flag if LoF in a gene with LOEUF < 0.6 or a recessive gene, or if REVEL ≥ 0.644 and AlphaMissense > 0.564 agree, or if REVEL alone reaches the "strong" level (≥ 0.932), or (when REVEL has no score) AlphaMissense ≥ 0.9. Always **Limited**: predictors are about 90% accurate on known variants, which still yields many false alarms across a genome. Thresholds follow Pejaver et al. 2022 (ClinGen SVI calibration).
- **Rebuilds:** changing a knowledge source's parsed `schema` forces a rebuild on the next refresh; changing claims logic bumps the annotate stage version so `wgs run` re-annotates automatically.
- **Alternatives considered:** OMIM (best inheritance data, but licensed — HPO/Orphanet cover the same diseases openly); automated ACMG classification such as InterVar or CPSR (over-calls without lab evidence and hides the reasoning); CADD instead of REVEL (genome-wide but less calibrated for missense; its licence also forbids redistribution); SpliceAI (useful, but the precomputed scores need an Illumina login — a candidate for later).
