// Dashboard sections, in navigation order. `step` = plan step that delivers it (docs/PLAN.md).
export interface SectionDef {
  id: string; title: string; icon: string; module?: string; step?: number; ready?: boolean;
  blurb: string; evidence?: string; learn?: string;
}

export const SECTION_DEFS: SectionDef[] = [
  { id: 'overview', title: 'Overview', icon: '🏠', ready: true, blurb: 'Start here: what is in this release and where to go next.' },
  { id: 'pgx', title: 'Medicines', icon: '💊', module: 'pgx', ready: true, blurb: 'How your genes affect the way you process common medicines (pharmacogenomics), with the published prescribing guideline behind each result.', evidence: 'CPIC and DPWG guidelines via ClinPGx; PharmCAT star-allele calling; CYP2D6 copy number from reads.', learn: 'Pharmacogenomics' },
  { id: 'health', title: 'Health', icon: '❤️', module: 'health', ready: true, blurb: 'Variants with documented links to disease: the 84 medically actionable genes, single-gene findings, computer predictions and common risk factors — each with its evidence grade.', evidence: 'ClinVar (per-lab submissions), ACMG SF v3.3, ClinGen validity & actionability, Mondo/HPO/Orphanet, gnomAD, AlphaMissense, REVEL.', learn: 'From variant to meaning' },
  { id: 'carrier', title: 'Carrier', icon: '🧬', module: 'carrier', ready: true, blurb: 'Recessive conditions where you carry one altered copy — usually no effect on you, but relevant for family planning.', evidence: 'ClinVar pathogenic / likely pathogenic with review stars and per-lab conditions; HPO/Orphanet inheritance; predictions for recessive genes.', learn: 'Inheritance patterns' },
  { id: 'traits', title: 'Traits', icon: '🧑', module: 'traits', step: 6, blurb: 'Physical and behavioural traits — eye colour, lactose tolerance, caffeine, height and more — with how strong the genetic evidence really is.', evidence: 'GWAS Catalog, PGS Catalog scores, well-replicated single-variant traits.', learn: 'Monogenic vs polygenic' },
  { id: 'ancestry', title: 'Ancestry', icon: '🌍', module: 'ancestry', step: 7, blurb: 'Where your ancestors likely lived, plus your maternal (mitochondrial) and paternal (Y) lineages.', evidence: '1000 Genomes + HGDP reference panel; PhyloTree (mtDNA) and ISOGG (Y) haplogroup trees.', learn: 'Ancestry' },
  { id: 'body', title: 'Body map', icon: '✨', step: 8, blurb: 'A fun infographic of the most notable ways your genome differs from the average human, pinned to the part of the body each gene affects.', learn: 'Percentile' },
  { id: 'changes', title: 'What changed', icon: '🔄', ready: true, blurb: 'How your results moved since the previous release when ClinVar, gene or population data were updated — and why.', learn: 'Keeping results current' },
  { id: 'qc', title: 'Data quality', icon: '✅', module: 'qc', ready: true, blurb: 'Can this data be trusted? Coverage, read quality, variant-call sanity checks, sex check and reference proof.', learn: 'Reading the QC report' },
  { id: 'variants', title: 'Variants', icon: '🔎', ready: true, blurb: 'Browse and search every variant in your genome, with a plain-language reading of each.' },
  { id: 'sql', title: 'SQL', icon: '⌨️', ready: true, blurb: 'Ask your own questions of the data with SQL, right in the browser.' },
  { id: 'learn', title: 'Learn', icon: '📖', ready: true, blurb: 'A friendly introduction to genetics and every term used in this dashboard.' },
];
