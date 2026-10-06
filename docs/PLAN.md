# Development plan

Each step ends with a **checkpoint** where Daniel tries the result before we move on.

| # | Step | Checkpoint |
|---|---|---|
| 0 | Scaffold: repo, pixi env, config and discovery, `wgs doctor`, data guard | `pixi run wgs doctor` lists the WGS files and section readiness |
| 1 | Ingest and QC: VCF → columnar store, raw genotype concordance, VCF QC metrics, reconstruct and verify the CRAM reference, coverage | QC numbers reviewed |
| 2 | Viewer shell and launchers (Windows x64, macOS universal), QC tab, variant explorer, SQL console, Learn tab | Dashboard opens by double-click on both laptops |
| 3 | Knowledge layer, evidence and confidence model, dated releases and diffs | Review the evidence model on real variants |
| 4 | Health risk and carrier status | Section review |
| 5 | Pharmacogenomics (PharmCAT, CYP2D6 via Cyrius) | Compare with tellmeGen |
| 6 | Traits and polygenic score percentiles, plus a "Brain & mind" deep-dive (autism/ADHD research explainer, SFARI gene check, polygenic percentiles with caveats) | Section review |
| 7 | Ancestry and haplogroups | Section review |
| 8 | "You vs the average human" body-map infographic | Visual review |
| 9 | GIAB HG002 benchmark and demo mode | Accuracy numbers |
| 10 | End-to-end refresh: "What changed" | Simulated knowledge update |
| 11 | Documentation, Windows test, polish; AI helpers: a local read-only MCP server (`wgs mcp`) and an "Ask AI" button that copies a claim plus its evidence as a prompt | Final sign-off |

## Future options (not scheduled)

- **Parent of origin ("who gave you what").** Only worthwhile once a parent is tested; a consumer array file (23andMe or AncestryDNA) from one parent is enough. Import it as an optional input and run trio phasing: wherever the parent has only one version, it shows which of your two copies came from them. That labels most variants on chromosomes 1–22 as maternal or paternal, settles cis/trans questions such as the *GALT* pair, and lets the Ancestry section paint each parent's side separately. Without a parent, read-backed phasing (WhatsHap on the CRAM) can only split variants into two unlabelled copies over short stretches.
