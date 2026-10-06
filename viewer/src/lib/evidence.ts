// Plain-language helpers for annotations and graded claims (evidence model: docs/DECISIONS.md ADR-013).
import type { Row } from './db';

export const OVERALL_PILL: Record<string, string> = { Strong: 'pass', Moderate: 'partial', Limited: '', 'Not callable': 'fail' };
export const CALL_PILL: Record<string, string> = { High: 'pass', Medium: 'partial', Low: 'warn', 'Not callable': 'fail' };
export const IMPACT_PILL: Record<string, string> = { HIGH: 'fail', MODERATE: 'partial', LOW: '', MODIFIER: 'soon' };

export const IMPACT_WORD: Record<string, string> = {
  HIGH: 'likely to break the protein',
  MODERATE: 'changes the protein',
  LOW: 'unlikely to change the protein',
  MODIFIER: 'outside the protein-coding sequence',
};

const CSQ: Record<string, string> = {
  missense: 'swaps one amino acid in the protein',
  synonymous: 'changes the DNA but not the protein (a “silent” change)',
  stop_gained: 'adds a premature stop signal, usually truncating the protein',
  stop_lost: 'removes the stop signal, so the protein runs on',
  start_lost: 'removes the start signal of the protein',
  frameshift: 'shifts the reading frame, scrambling everything after it',
  inframe_insertion: 'adds amino acids without shifting the reading frame',
  inframe_deletion: 'removes amino acids without shifting the reading frame',
  splice_donor: 'hits the edge of an exon, likely disrupting how the gene is spliced',
  splice_acceptor: 'hits the edge of an exon, likely disrupting how the gene is spliced',
  splice_region: 'sits near an exon edge and might affect splicing',
  '5_prime_utr': 'sits in the untranslated start of the gene’s message (can affect how much protein is made)',
  '3_prime_utr': 'sits in the untranslated tail of the gene’s message (can affect its stability)',
  intron: 'sits in an intron, the stretch cut out before the protein is made',
  non_coding: 'sits in a gene that makes RNA rather than protein',
  stop_retained: 'touches the stop signal but keeps it',
  start_retained: 'touches the start signal but keeps it',
  coding_sequence: 'sits in the coding sequence',
  intergenic: 'sits between genes',
};

/** "missense&NMD_transcript" → plain words. */
export function consequenceWord(c: string | null | undefined): string {
  if (!c) return 'sits between genes';
  const parts = c.split('&').filter((p) => p !== 'NMD_transcript');
  return parts.map((p) => CSQ[p] ?? p.replaceAll('_', ' ')).join('; and ');
}

export const consequenceLabel = (c: string | null | undefined) =>
  !c ? 'intergenic' : c.split('&').filter((p) => p !== 'NMD_transcript').map((p) => p.replaceAll('_', ' ').replace('prime utr', '′ UTR')).join(' + ');

/** bcftools csq "222A>222V" → "A222V". */
export function aaShort(aa: string | null | undefined): string {
  if (!aa) return '';
  const m = /^(\d+)([A-Z*]+)>(\d+)([A-Z*]+)$/.exec(aa);
  return m ? `${m[2]}${m[1]}${m[4]}` : aa;
}

export const CLASS_LABEL: Record<string, string> = {
  pathogenic: 'Pathogenic', likely_pathogenic: 'Likely pathogenic', uncertain: 'Uncertain significance',
  likely_benign: 'Likely benign', benign: 'Benign', conflicting: 'Conflicting', drug_response: 'Drug response',
  risk_factor: 'Risk factor', association: 'Association', protective: 'Protective', other: 'Other',
  predicted: 'Predicted damaging',
};
export const CLASS_PILL: Record<string, string> = {
  pathogenic: 'fail', likely_pathogenic: 'fail', conflicting: 'warn', uncertain: 'warn', risk_factor: 'partial',
  drug_response: 'partial', association: '', protective: 'pass', benign: 'pass', likely_benign: 'pass', other: 'soon',
  predicted: 'soon',
};

/** What your genotype means given the inheritance (evidence model v2): label, pill class. */
export const ROLE: Record<string, [string, string]> = {
  carrier: ['Carrier', 'pass'], affected: ['Both copies affected', 'fail'], possible: ['May matter', 'warn'],
  unknown: ['Unclear', 'soon'],
};

/** Inheritance mode codes → glossary terms. */
export const INH: Record<string, string> = {
  AD: 'Autosomal dominant', AR: 'Autosomal recessive', SD: 'Semi-dominant', XLR: 'X-linked recessive',
  XLD: 'X-linked dominant', XL: 'X-linked', YL: 'Y-linked', MT: 'Mitochondrial', MF: 'Multifactorial',
  DG: 'Digenic', SP: 'Sporadic',
};

export const STAR_WORD: Record<number, string> = {
  0: 'no assertion criteria provided', 1: 'one submitter with criteria, or conflicting submitters',
  2: 'multiple submitters agree, no conflicts', 3: 'reviewed by an expert panel', 4: 'practice guideline',
};
export const stars = (n: number | null | undefined) => (n == null ? '' : '★'.repeat(n) + '☆'.repeat(4 - n));

export function pct(af: number | null | undefined): string {
  if (af == null) return '—';
  if (af === 0) return '0%';
  if (af < 0.0001) return '<0.01%';
  if (af < 0.01) return `${(af * 100).toFixed(2)}%`;
  return `${(af * 100).toFixed(af < 0.1 ? 1 : 0)}%`;
}

/** How common an allele is, in words. */
export function rarity(af: number | null | undefined): string {
  if (af == null) return 'not seen in the reference panel';
  if (af >= 0.5) return 'the majority version';
  if (af >= 0.05) return 'common';
  if (af >= 0.01) return 'uncommon';
  if (af >= 0.001) return 'rare';
  return 'very rare';
}

/** For a carrier of one/two copies, the share of people with at least one copy (Hardy–Weinberg). */
export const carriersOf = (af: number) => 1 - (1 - af) ** 2;

export const CONTINENTS: [string, string][] = [
  ['afr', 'African'], ['amr', 'Admixed American'], ['eas', 'East Asian'], ['eur', 'European'], ['sas', 'South Asian'],
];

export function parseSources(r: Row): { source: string; version?: string; record?: string; url?: string }[] {
  try { return JSON.parse(r.sources ?? '[]'); } catch { return []; }
}

export const SECTION_OF_CLAIMS: Record<string, string> = { health: 'health', carrier: 'carrier', pgx: 'pgx', traits: 'traits' };
