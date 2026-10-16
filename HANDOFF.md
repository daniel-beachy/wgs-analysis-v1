# Handoff (temporary)

Written 2026-10-15 by the previous agent (GitHub Copilot) for the next one
(Claude Code). It captures everything that lived only in that agent's session
memory. **Once you have absorbed it into your own notes / `CLAUDE.md` /
`docs/`, delete this file.**

This repo is **public**. Never commit personal results or the sample ID (the
pre-commit data guard blocks the ID; `pixi run guard` checks everything).
Personal context lives outside the repo in `../wgs-notes/` (see its README).

## 1. Read in this order

1. This file.
2. `AGENTS.md`: data layout, tables, query tool, how to talk about results.
3. `docs/PLAN.md`: steps and checkpoints (what's done, what's next).
4. `docs/DECISIONS.md`: ADR-001…018. ADR-016 (no hand-written medical text) and ADR-018 (audit, absolute risk) govern almost every change.
5. `README.md`: setup, CLI, knowledge sources and licences.
6. `../wgs-notes/session-history/checkpoints/` (private): full history, newest 024–029 most relevant.

## 2. The original brief (Daniel, 2026-10-05), condensed

- Turn a personal 30× tellmeGen WGS (FASTQ, CRAM, VCF, CSVs in `Genomics/WGS`) into something he can explore and learn from **for a decade**.
- Point at files anywhere (NVMe, SSD, any path per machine). **Degrade gracefully** when files are missing. Ingest outputs of Linux-only tools if present, never require them.
- Fast repeated queries; a **durable, portable on-disk artifact**.
- **Interpretation is the hard part.** Explicit about what evidence underpins each claim and how confident it is. Every finding links to its exact source.
- **Not a snapshot.** Refresh evidence as knowledge changes and **show what changed** in his results.
- Dashboard sections: pharmacogenomics, health/disease risk, carrier status, physical and behavioural traits, ancestry, data quality/QC.
- Heavy compute once on the fast Mac, written to the external SSD. The dashboard must be viewable from a **Windows x86 laptop and an Apple Silicon Mac** on that same drive; nothing architecture-locked; setup genuinely simple.
- **Fun infographic** ("you vs the average human"): body map with blurbs pointing at heart/brain/muscles, showing the variant and his percentile.
- Must also run on **Genome in a Bottle (HG002)** to check accuracy against truth.
- Research best practice, document why, deliver runnable code, **visually verify the app yourself** (screenshots) and iterate.
- GitHub account `daniel-beachy`. Repo name `wgs-analysis-v[n]` (he typed "wsg", meaning wgs; v1 is the current one).
- Hosting: **decided local-only for now**. If ever hosted: public GitHub Pages; if dynamic, Cloudflare free tier only (Pages/Workers/D1, no billing, use the locally authenticated wrangler, never interactive auth). Portfolio card on https://daniel-beachy.github.io (dark thumbnail, same size as the others, link back from the app) is **deferred** until hosting is decided.

Later standing requests:
- Learning-first: explanations, tooltips and glossary everywhere (ADR-010, `docs/GENETICS_PRIMER.md`).
- Sensitive findings are click-to-reveal.
- Brain & mind section (done in Step 6).
- **Verification standard (Step 6b):** every finding links to its exact study; impact shown in absolute terms ("26 in 100 vs 13 typical") where sourced, otherwise "not established".

## 3. How Daniel likes to work

- Sequential steps, each ending in a **checkpoint** he reviews before moving on. Ask before big scope changes. Present choices as a short list with a recommendation.
- He's a layperson in genetics but technical. Use plain language, short answers and concrete numbers. He often asks "explain the methodology / how sure is this".
- **ADR-016 is strict.** No hand-written medical text and no one-off spot fixes. Every statement is generated from a cited source by a rule. Every bug fix gets a self-check (`src/wgs/checks.py`, runs at publish) or a test that fails on the old code. Document rule changes as an ADR or addendum.
- Verify UI changes with real screenshots (Playwright) before reporting done.
- Commits: personal identity (already configured). The previous agent added a `Co-authored-by` trailer for itself; use whatever trailer is appropriate for you. Delete macOS `._*` files before committing (exFAT drive): `find . -name '._*' -not -path './.pixi/*' -delete`.

## 4. Environment facts and gotchas

- Machine: M1 Pro MacBook, **16 GB RAM**, macOS. No Linux machine (ADR-004). The SSD is **exFAT** (`._*` files appear, no symlinks or Unix permissions; build scratch goes to `$TMPDIR/wgs-build`).
- Env: `pixi` (conda-forge + bioconda). The envs live in `.pixi/envs/` (default, pgs, cyrius, fraposa) on the SSD. Shortcut: `export PATH=$PWD/.pixi/envs/default/bin:$PATH`. plink2 and the `pgscatalog-*` tools are in `.pixi/envs/pgs/bin`.
- Hooks: `git config core.hooksPath scripts/hooks` (already set; `pixi run install-hooks` on a fresh clone).
- Local config `wgs.local.toml` (git-ignored) holds paths and the sample ID. Template: `config/wgs.example.toml`.
- Data: `../wgs-data/` (knowledge store, work dirs, immutable releases). Raw `../WGS/` is read-only. Never write there.
- Commands: `wgs doctor`, `wgs run [--only X] [--force X]`, `wgs knowledge status|refresh`, `wgs query "<SQL>"` (`-f json`; use `COLUMNS=250` for wide tables), `wgs audit [--only links|curated]`, `pixi run build-dashboard`, `pixi run test` (103 tests), `pixi run lint`, `pixi run guard`.
- Dashboard: double-click `../Genome Dashboard.app` / `.exe` (Go launcher serving the static Svelte + DuckDB-WASM viewer). Headless: `"../Genome Dashboard.app/Contents/MacOS/genome-dashboard" -port 8791 -no-browser -idle 0`. Stop it with `kill <PID>`.
- Screenshot verification: the previous agent used a throwaway Node Playwright script (lost from `$TMPDIR`). Recreate one, or use a Playwright MCP, against the headless launcher.
- Long jobs: annotate takes ~15 min. The PGS stage's 1000 Genomes reference scoring can take hours (plink2 over 223 scores × 2,504 samples). Run in the background, and confirm with Daniel before starting multi-hour runs.
- DuckDB query gotchas:
  - `pgx_drugs.action` is BOOLEAN.
  - `trait_snps` has no `trait` column (use `label` / `trait_reported`).
  - `brain_variants` has no `overall` column.
  - QC metrics are in the `doc_qc` JSON (`sections.vcf/concordance/coverage`).
  - No `ilike any(array)`; use `regexp_matches`.
- PGS internals (`src/wgs/stages/pgs.py`): only the IDs in `knowledge/pgs_catalog/<date>/selected.parquet` are scored, and the target genotype file covers only their positions. `_score` deletes cached files matching `b0_…`, so don't share `score-ref-cache` between side runs. A side run that scores extra IDs by reusing `_reference_qc → _target → _pca → _match → _score → _adjust` with an isolated work dir worked (≈10 min for 3 scores).

## 5. Status (2026-10-15)

- **Done:** Steps 0–6 and 6b. Latest release `2026-10-07_214605`: 28/28 self-checks, 103 tests, live audit clean (one link blocked by Science's 403). HEAD `cc2cd08`, pushed, tree clean.
- **Rough share of the remaining work** (project ≈ 65–70% done): 7 Ancestry ~8%, 7b Score spread ~3%, 8 Body map ~7%, 9 GIAB ~8%, 10 Refresh end-to-end ~5%, 11 Polish/Windows/AI helpers ~7%.

## 6. Open decisions (ask Daniel)

1. **(Decided 2026-10-16, ADR-019.) Step 7b headline value.** Once all scores per trait are shown as a range, what single number drives absolute risk, body-map ranking and "what changed" diffs? The proposal was the median percentile for the headline and ranking, and the most-tested score for absolute risk (labelled), not yet confirmed. A related, not yet approved rule fix in `src/wgs/knowledge/pgs.py::select`: prefer scores with more independent evaluations over one large self-evaluation, and flag a reported R² above the trait's SNP heritability (needs a cited heritability source). Example that triggered this: PGS004427 (1 evaluation, R² 0.22 implausible) was picked over PGS003510 (9 evaluations).
2. **Step 7 ancestry design.** Global ancestry from the existing 1000 Genomes PCA (already computed for PGS) and/or HGDP + 1000G for finer resolution; admixture method (supervised ADMIXTURE-like, or rye / projection); haplogroups via haplogrep3 (mtDNA) and yhaplo (Y). Check that both run on macOS arm64 via pixi. Decide how fine-grained to claim and how to grade confidence (same evidence model).
3. **Step 8 body map without hand-written text (ADR-016).** Which findings qualify (featured PGS percentiles, PGx, monogenic), how organ placement is sourced (e.g. a trait→organ mapping from EFO/UBERON or HPO anatomy, not hand-assigned), and how ranges from 7b are shown.
4. **Step 9 GIAB.** There's no variant caller on this Mac (DeepVariant is Linux/x86), so "accuracy" likely means: (a) run the full interpretation pipeline on the HG002 truth VCF as a demo mode, and (b) compare tellmeGen-style calls only if a caller output exists. Use `rtg vcfeval` (works on macOS) rather than hap.py. Settle what is being benchmarked.
5. **Known limitation:** Leiden-type variants where GRCh37's reference carries the disease allele are absent from ClinVar's GRCh37 VCF, and the VCF is not a gVCF. A fix would need ClinVar `variant_summary` plus callable/CRAM evidence. Currently documented, not solved.
6. **Hosting / portfolio card:** still local-only. Revisit at Step 11.

## 7. Next action

Ask Daniel about decision 1, then 2. Then build Step 7 (Ancestry). Step 7b can come before or after Step 7, but must land before Step 8.
