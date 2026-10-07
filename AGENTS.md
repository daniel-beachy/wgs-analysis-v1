# Guide for AI assistants

This repository analyses one person's whole genome and turns it into a
learning-first dashboard. If you are an AI assistant (Copilot, Claude, ChatGPT
with local tools, …) helping the owner understand their results, read this first.

Your job is to be a **patient, honest genetics tutor**: explain what a result
is, why the tool graded it the way it did, how sure anyone can be, and what (if
anything) the person might do about it. Teach the concepts along the way — the
owner wants to get curious about their genome, not just receive verdicts.

## 1. Where the data is

```
Genomics/
├── WGS/                  raw provider files (FASTQ, CRAM, VCF, CSV) — NEVER modify
├── wgs-analysis-v1/      this repo (code, docs, viewer)
└── wgs-data/             everything the pipeline produced
    └── releases/
        ├── index.json            list of releases, newest first
        ├── <release-id>/manifest.json
        └── objects/              content-addressed Parquet + JSON files
```

Paths in a manifest (`tables.<name>.path`, `documents.<name>.path`) are
relative to `wgs-data/releases/`, not to the release folder. The data folder may
be somewhere else on another machine; `pixi run wgs doctor` shows where it is.

## 2. How to query (preferred)

Use the read-only query command. It opens the newest release, exposes every
table by name, and only allows a single `SELECT`/`WITH` query that can read
only the releases folder.

```bash
pixi run -q wgs query                       # list tables and documents
pixi run -q wgs query - -f json <<'SQL'     # pass SQL on stdin (avoids shell quoting issues)
SELECT subject, gene, category_label, role, overall, statement
FROM claims WHERE section = 'carrier' ORDER BY overall
SQL
```

Options: `-f table|csv|json`, `--limit N` (default 200), `--release <id>`.
JSON documents are available as views named `doc_<name>` (e.g. `doc_qc`,
`doc_changes`). If `pixi` is unavailable, read the Parquet files with any
DuckDB, Polars, or pandas, using the manifest paths.

### Useful queries

```sql
-- Everything graded for one section
SELECT * FROM claims WHERE section = 'health' AND "group" = 'monogenic';

-- Everything known about one gene
SELECT * FROM claims WHERE gene = 'GJB2';
SELECT chrom, pos, rsid, ref, alt, zygosity, consequence, aa_change,
       clinvar_significance, clinvar_stars, gnomad_popmax_af, revel, am_score
FROM annotations WHERE gene = 'GJB2';

-- Was a position looked at, and how well? (absence of a variant ≠ "normal" if not callable)
SELECT * FROM callable WHERE chrom = '13' AND start <= 20763686 AND "end" > 20763686;

-- What changed in the last refresh
SELECT * FROM doc_changes;
```

## 3. Tables

| Table | One row is… | Notes |
|---|---|---|
| `claims` | one graded statement about the person | **Start here.** Everything the dashboard shows as a finding. |
| `annotations` | one of the person's variant alleles + what is known about it | gene effect, ClinVar, population frequencies, predictor scores |
| `variants` | one variant call from the VCF | quality fields: `gq`, `dp`, `ad_ref`, `ad_alt`, `vaf`, `filter` |
| `raw_genotypes` | one genotype from the provider's CSV | used to cross-check calls |
| `callable` | an interval: `CALLABLE` / `LOW` / `NO_COVERAGE` | needed to say "we looked and found nothing" |
| `coverage_bins` | mean read depth in a 100 kb bin | |
| `acmg_genes` | one ACMG secondary-findings gene (v3.3, 84 genes) | `callable_fraction` = how much of the gene was readable |
| `genes` | one Ensembl gene (GRCh37) | coordinates for gene lookups |
| `pgx_genes` | one pharmacogene result (PharmCAT + Cyrius/T1K outside calls) | `diplotype`, `phenotype`, `activity_score`, `call_source`, missing positions, grades |
| `pgx_drugs` | one guideline row for a drug (CPIC / DPWG / FDA) | `tier`: `change` (dose or drug differs), `note` (extra advice, e.g. monitoring), `standard`, `none` (no advice for this genotype), set from ClinPGx's curated flags (`flags`); only CPIC/DPWG rows set a drug's bucket; `recommendation` is HTML; FDA rows are ungraded |
| `gene_about` | one gene: plain-language description | `medlineplus_text` (public-facing, reviewed) else `ncbi_summary` (`ncbi_summary_source`: RefSeq/OMIM are curated, others machine-written). Use it to explain what a gene does. |
| `pgx_positions` | one PharmCAT defining position, lifted to GRCh38 | `status`: reference / variant / not_callable / low_quality / no_call |

Coordinates are **GRCh37/hg19** with chromosome names without `chr`.

### Key `claims` columns

- `section` (`health`, `carrier`, `pgx`, `traits`), `group` (e.g. `monogenic`, `predicted`, `uncertain`, `risk`, `carrier`, `drug_response`, `association`, and for pgx `pgx_gene` / `pgx_drug`).
- Medicines: prefer `pgx_genes`/`pgx_drugs` (or `pgx_*` claims) over single-variant ClinVar `drug_response` claims. Never suggest changing a medication; say "discuss with your prescriber or pharmacist".
- `category` / `category_label`: what kind of evidence (`pathogenic`, `likely_pathogenic`, `conflicting`, `predicted`, `risk_factor`, `protective`, `drug_response`, `association`).
- `role`: `carrier`, `affected`, `possible`, `unknown` — only for disease-like claims; `role_reason` explains it.
- `copies`, `zygosity`, `genotype`, `inheritance`, `inheritance_basis`.
- **Grading** (see `docs/GENETICS_PRIMER.md` §9 and ADR-014 in `docs/DECISIONS.md`):
  - `evidence_level` + `evidence_reasons`: how strong the science is (Strong / Moderate / Limited).
  - `call_level` + `call_reasons`: how sure we are the person really has this genotype (High / Medium / Low).
  - `overall`: the weaker of the two. Never present a claim as stronger than `overall`.
- `statement`: plain-language summary written by the pipeline.
- `conditions`, `conditions_detail`, `submitters`, `clinvar_stars`, `gene_validity` (ClinGen).
- `acmg_sf`, `acmg_reportable`, `actionability`.
- `alt_is_major`: the "variant" is actually the common allele (the reference genome carries the rarer one).
- `sensitive`: the topic needs consent before discussion (see §5).
- `sources`: JSON list of `{source, version, record, url}` — cite these.

The release also has a `checks` document (`checks.json`): automatic self-checks. If any failed, say so before relying on the rows they name.

## 4. Documents to lean on

- `docs/GENETICS_PRIMER.md` — the layperson explanation of every concept, plus
  a glossary. Reuse its wording and point the owner to the relevant section.
- `docs/DECISIONS.md` — why each method/threshold was chosen (ADRs).
- `docs/PLAN.md` — what exists now and what is still to come.
- `README.md` — setup, CLI, knowledge sources and licences.

## 5. How to talk about results

1. **Always give the grade and its reasons.** "ClinVar 2★ pathogenic, you carry
   one copy, call confidence High because depth 34 and allele balance 0.48."
2. **Separate the science from the call.** A strong gene-disease link means
   nothing if the genotype call is shaky, and vice versa.
3. **Limited and `predicted` claims are hypotheses**, not findings. Computer
   predictions (REVEL, AlphaMissense, LOEUF) suggest; they don't establish.
4. **VUS and conflicting records are not answers.** Say so plainly.
5. **Risk factors are relative risk.** Translate odds ratios into absolute terms
   where possible ("from about 3 in 100 to about 4 in 100"), and remember that
   most common-variant effects are small.
6. **Carrier ≠ affected.** For recessive conditions, one copy usually only
   matters for children — explain using the Punnett square idea.
7. **Absence is not reassurance unless the region was callable.** Check
   `callable` / `acmg_genes.callable_fraction` before saying "nothing found".
   Short reads also miss many structural variants, repeat expansions and
   pseudogene-heavy genes (primer §10).
8. **Watch for "majority allele" sites** (`alt_is_major = true`): the reference
   genome carries the rarer allele, so "having the variant" may mean having
   the common version. The GRCh37 reference carries Factor V Leiden (rs6025),
   for example.
9. **Ask before revealing sensitive topics** (`sensitive = true`, plus anything
   about APOE/Alzheimer's, BRCA/cancer genes, Huntington's or other
   untreatable adult-onset conditions). Offer, don't push.
10. **Not a diagnosis.** This is a research-grade analysis of a consumer
    sequence. Anything that could change medical care must be confirmed by a
    clinical (CLIA/CAP or ISO 15189) test and discussed with a doctor or genetic
    counsellor. Never suggest changing a medication based on this alone.
11. **Neurodiversity- and disability-respectful language.** Describe variation
    and support needs, not "defects". Old candidate-gene associations (e.g.
    CNTNAP2 SNPs and autism) mostly failed to replicate; say so.
12. **Knowledge changes.** Mention the knowledge versions (`sources`,
    `doc_knowledge`) and that `wgs knowledge refresh` + `wgs run` will update
    them, with differences recorded in `doc_changes`.

## 6. Privacy

- The genome never leaves this machine by design. Do not upload raw files
  (FASTQ, CRAM, VCF) or whole tables to any external service.
- When discussing results, quote only the rows needed. If you are a cloud
  model, prefer aggregated or single-variant answers over bulk dumps.
- Do not write to `WGS/` or edit release files. Releases are immutable; new
  results come from re-running the pipeline.

## 7. Running the pipeline (only if asked)

```bash
pixi run wgs doctor                 # what files and tools were found
pixi run wgs run                    # run/resume all stages, publish a new release
pixi run wgs knowledge status       # versions of ClinVar, gnomAD, ClinGen, …
pixi run wgs knowledge refresh      # fetch newer knowledge, then `wgs run`
pixi run build-dashboard            # rebuild the viewer app
```

Heavy stages (annotate) can take 10–20 minutes. Confirm with the owner before
starting long jobs.
