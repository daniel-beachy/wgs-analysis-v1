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
| 11 | Documentation, Windows test, polish | Final sign-off |
