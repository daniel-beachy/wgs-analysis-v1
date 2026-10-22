// "Ask AI" (ADR-027): one result, its grades with their reasons and its sources, turned into a prompt to paste into
// any assistant. Nothing is sent anywhere; the person chooses where to paste it. The rules at the end condense
// AGENTS.md §5 (how to talk about results), so a general-purpose assistant answers the way the project asks.
import type { Manifest } from './data';
import { query, sqlString, type Row } from './db';
import { parseSources } from './evidence';

export interface AskItem {
  title: string;
  kind?: string;
  statement?: string;
  why?: string;
  overall?: string;
  evidence?: string;
  evidenceReasons?: string[];
  call?: string;
  callReasons?: string[];
  sources?: { source: string; version?: string; url?: string }[];
  ref?: string;
}

// The table and key column a highlight points at, mirroring take.CITABLE in the pipeline.
const KEY: Record<string, string> = {
  claims: 'claim_id', pgx_genes: 'gene', pgs_traits: 'trait_id', trait_snps: 'rsid', hla_associations: 'association_id',
  special_genes: 'gene', repeats: 'locus_id', cnv_events: 'event_id',
};

const RULES = [
  'Explain what this means in plain language: I am technical but not a geneticist.',
  'Keep the evidence grade (how strong the science is) separate from the call confidence (how sure we are I really have this genotype), and never present the result as stronger than its overall grade.',
  'Give risks in absolute terms ("from about 3 in 100 to about 4 in 100") only where a source gives the numbers; otherwise say the absolute effect is not established.',
  'One copy for a recessive condition means a carrier, not affected. Limited results and computer predictions are hypotheses, and variants of uncertain significance are not answers.',
  'This is a research-grade analysis, not a diagnosis. Anything that could change medical care needs a clinical (CLIA/CAP or ISO 15189) test and a doctor or genetic counsellor. Never suggest changing a medication: say to discuss it with a prescriber or pharmacist.',
  'Cite a source for anything you add, and say when you are unsure.',
];

const list = (xs: string[] | undefined) => (xs ?? []).filter(Boolean).map((x) => `  - ${x}`);

export function promptFor(it: AskItem, m: Manifest): string {
  const whose = m.truth_set
    ? `a public reference genome (${m.truth_set.label}), not mine`
    : 'my own genome (a 30× consumer whole genome)';
  const result = [
    `Result: ${it.title}${it.kind ? ` (${it.kind})` : ''}`,
    it.statement && `What it says: ${it.statement}`,
    it.why && `Why it was highlighted: ${it.why}`,
    it.overall && `Overall grade: ${it.overall} (the weaker of the two grades below; Strong / Moderate / Limited)`,
    it.evidence && `Evidence grade: ${it.evidence}`, ...list(it.evidenceReasons),
    it.call && `Call confidence: ${it.call}`, ...list(it.callReasons),
    it.sources?.length && 'Sources:',
    ...list(it.sources?.map((s) => `${s.source}${s.version ? ` ${s.version}` : ''}${s.url ? `: ${s.url}` : ''}`)),
    it.ref && `(In the analysis: ${it.ref}.)`,
  ].filter(Boolean);
  return [
    `I'm exploring a research-grade analysis of ${whose}, release ${m.id}. Please help me understand this one result.`,
    result.join('\n'),
    ['How to answer:', ...RULES.map((r) => `- ${r}`)].join('\n'),
  ].join('\n\n');
}

const arr = (v: any): string[] | undefined => (v == null ? undefined : Array.isArray(v) ? v : Array.from(v));

export function claimItem(r: Row): AskItem {
  return {
    title: [r.gene ?? r.subject, r.rsid].filter(Boolean).join(' '), kind: r.category_label, statement: r.statement,
    overall: r.overall, evidence: r.evidence_level, evidenceReasons: arr(r.evidence_reasons),
    call: r.call_level, callReasons: arr(r.call_reasons), sources: parseSources(r), ref: `claims, claim_id ${r.claim_id}`,
  };
}

/** A highlight carries only its grades; the reasons live on the row it came from, so look that row up. */
export async function highlightItem(r: Row): Promise<AskItem> {
  const it: AskItem = {
    title: r.title, kind: r.kind_label, statement: r.statement, why: r.notable_because, overall: r.overall,
    evidence: r.evidence_level, call: r.call_level, sources: parseSources(r), ref: `${r.ref_table}, ${KEY[r.ref_table] ?? 'key'} ${r.ref_key}`,
  };
  const key = KEY[r.ref_table];
  if (!key) return it;
  try {
    const src = (await query(`SELECT * FROM ${r.ref_table} WHERE CAST(${key} AS VARCHAR) = ${sqlString(String(r.ref_key))} LIMIT 1`)).rows[0];
    if (src) { it.evidenceReasons = arr(src.evidence_reasons); it.callReasons = arr(src.call_reasons); }
  } catch { /* the grades alone still make a useful prompt */ }
  return it;
}

/** Copy text; falls back to a hidden textarea where the async clipboard API is unavailable. */
export async function copyText(text: string): Promise<boolean> {
  try { await navigator.clipboard.writeText(text); return true; } catch { /* fall back below */ }
  const ta = document.createElement('textarea');
  ta.value = text; ta.setAttribute('readonly', ''); ta.style.position = 'fixed'; ta.style.opacity = '0';
  document.body.appendChild(ta); ta.select();
  const ok = document.execCommand('copy');
  ta.remove();
  return ok;
}
