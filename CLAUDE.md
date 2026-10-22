# Claude Code instructions

@AGENTS.md

`AGENTS.md` is the shared guide for any AI assistant, so keep it tool-neutral and put Claude-specific notes here.
Steps and their status are in `docs/PLAN.md`; every method and threshold is justified in `docs/DECISIONS.md`
(ADR-001 onwards). Read both before changing anything.

## The brief, condensed

Turn one person's 30× whole genome (tellmeGen: FASTQ, CRAM, VCF, CSVs in `../WGS`) into something he can explore
and learn from for a decade:
- Point at files anywhere; degrade gracefully when files or tools are missing; never require Linux-only tools.
- A durable, portable on-disk artifact (immutable releases) and a dashboard that runs from the external drive on an
  Apple Silicon Mac (Windows was dropped on 2026-10-22; the launcher still cross-compiles an untested `.exe`).
- Interpretation is the hard part: every finding says what evidence underpins it, how sure we are, and links to its
  exact source; impact in absolute terms ("26 in 100 vs 13 typical") where sourced, otherwise "not established".
- Not a snapshot: knowledge refreshes, and "What changed" shows the difference.
- Learning-first: plain language, tooltips and the primer/glossary everywhere. The owner is technical but not a
  geneticist.
- Every place the genome reliably differs from a typical one, what it means and how sure we are, Strong to Limited.
- Also run on Genome in a Bottle (HG002) to check accuracy.
- Hosting: the dashboard itself is local-only; the HG002 demo is public on GitHub Pages (ADR-028).

## How to work here

- **The repo is public.** Never commit personal results or the sample ID; the pre-commit hook and `pixi run guard`
  block the ID. Personal context lives outside the repo in `../wgs-notes/`. That includes docs, tests and code
  comments: never use the owner's own genes, variants or genotypes as examples or as an ADR's "first result"
  (aggregate counts are fine). The guard can't spot these, so check by eye.
- Work in sequential steps; each ends with a checkpoint the owner reviews. Commit and push verified work without
  asking. Ask before big scope changes and multi-hour runs; present choices as a short list with a recommendation.
- **ADR-016 is strict:** no hand-written medical text; every statement comes from a cited source by a rule. Curated
  tables (`src/wgs/knowledge/data/*.toml`) quote their source word for word and are re-checked by
  `wgs audit --only curated` (numbers and first author). Every bug fix gets a self-check (`src/wgs/checks.py`) or a
  test that fails on the old code. Document rule changes as an ADR or addendum.
- Verify citations (PMID, first author, volume) against PubMed before writing them; never from memory.
- Verify UI changes in the built-in browser pane (screenshots, dark and light, phone width) before reporting done.
- Before committing: `find . -name '._*' -not -path './.pixi/*' -delete` (exFAT), `pixi run lint`, `pixi run test`,
  `pixi run guard`.

## Environment

- M1 Pro MacBook, 16 GB RAM, macOS; the SSD is exFAT (`._*` files appear; no symlinks; build scratch goes to
  `$TMPDIR/wgs-build`). Data: `../wgs-data/` (knowledge store, work dirs, immutable releases); raw `../WGS/` is
  read-only.
- pixi environments in `.pixi/envs/`: `default`, `pgs` (plink2, pgscatalog-utils), `fraposa`, `cyrius`, `haplo`
  (haplogrep3, Yleaf from its GitHub tag) and `special` (SMNCopyNumberCaller, ExpansionHunter, Delly, Gauchian from
  its GitHub tag). Illumina's callers are pinned to their pre-December-2022 releases (Apache-2.0/GPL-3.0; ADR-021).
- Local config `wgs.local.toml` (git-ignored) holds paths and the sample ID; template `config/wgs.example.toml`.
- Commands: `wgs doctor`, `wgs run [--only X] [--force X]`, `wgs knowledge status|refresh [--source X]`,
  `wgs query - -f csv <<'SQL' … SQL` (stdin avoids quoting problems), `wgs audit [--only …]`,
  `pixi run build-dashboard`.
- Dashboard: `../Genome Dashboard.app` (Go launcher serving the Svelte + DuckDB-WASM viewer, which is embedded, so
  restart after each build). Headless:
  `"../Genome Dashboard.app/Contents/MacOS/genome-dashboard" -port 8791 -no-browser -idle 0`, then open
  `http://127.0.0.1:8791/` in the built-in browser pane (use 127.0.0.1, not localhost). Stop with `pkill -f
  "genome-dashboard -port 8791"`.
- Stage run times: annotate ≈ 7–15 min; pgs regrade seconds (the 1000 Genomes reference scoring is hours: confirm
  first); ancestry ≈ 40 s; special ≈ 5 min; publish ≈ 10 s.
- Pace work to the Pro plan: check `get_usage` at step boundaries and stop at committed points.

## Gotchas

- DuckDB: `pgx_drugs.action` is BOOLEAN; `trait_snps` has no `trait` column (use `label` / `trait_reported`);
  `brain_variants` has no `overall`; QC metrics live in `doc_qc` (`sections.vcf/concordance/coverage`); no
  `ilike any(array)` (use `regexp_matches`); `->>` precedence breaks WHERE clauses (use `json_extract_string`);
  `QUALIFY` needs a window function; pandas reads empty CSV fields as NaN (clean before writing parquet).
- The VCF names contigs `chrM`/`chrY`/`chr1`…; its index lives in `work/<sample>/variants/input.vcf.gz.tbi`
  (`variants.vcf_with_index`). `variants.parquet` uses names without `chr`. The CRAM is `chr`-prefixed.
- GRCh37 reference alleles that are the rarer or disease version: F5 Leiden (rs6025) and the ABO O allele (rs8176719).
  "No variant" can mean "has the reference's unusual allele".
- `samtools idxstats` on a CRAM reads the whole 53 GB file: cut regions through the index first (as the Y and
  ACMG-panel steps do).
- PGS internals (`src/wgs/stages/pgs.py`): only IDs in `knowledge/pgs_catalog/<date>/selected.parquet` are scored and
  the target genotype file covers only their positions; `_score` deletes cached files matching `b0_…`, so don't share
  `score-ref-cache` between side runs.

## Open decisions and known gaps

- **Demo (HG002):** polygenic scores only run with `--pgs` (hours of 1000 Genomes reference scoring: confirm
  first).
- Leiden-type variants where GRCh37's reference carries the disease allele are absent from ClinVar's GRCh37 VCF, and
  the VCF is not a gVCF; a fix needs ClinVar `variant_summary` plus callable/CRAM evidence.
- Polygenic scores: the selection-rule change (prefer more independent evaluations) and an R²-above-heritability flag
  need a new multi-hour scoring run; scoring every catalog score per trait likewise (ADR-019).
- Finer-scale ancestry (HGDP + 1000 Genomes) is a later, separate step (ADR-020).
- The public demo (ADR-028): `pixi run demo` rebuilds the HG002 release, then `pixi run demo-site` exports a trimmed
  copy (truth-set releases only), uploads it as a `demo-data-<id>` release asset and deploys
  https://daniel-beachy.github.io/wgs-analysis-v1/. The portfolio card lives in `~/AI_PROJECTS/portfolio-site`
  (follow its AGENTS.md).
