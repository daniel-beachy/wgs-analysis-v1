<script lang="ts">
  // Graded claims for one dashboard section, each showing what is claimed, how strong the evidence is,
  // how confident your genotype call is, and the sources behind it. Sensitive findings stay hidden until revealed.
  import type { Manifest } from '../lib/data';
  import { query, useRelease, sqlString, type Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { OVERALL_PILL, CALL_PILL, CLASS_PILL, CLASS_LABEL, stars, STAR_WORD, pct, parseSources, consequenceLabel, aaShort, ROLE, INH } from '../lib/evidence';
  import Term from './Term.svelte';

  let { m, section, intro = '', groups = [], empty = 'No graded findings for this section in this release.', filters = true, limit = 0, onrows }:
    { m: Manifest; section: string; intro?: string; groups?: string[]; empty?: string; filters?: boolean; limit?: number; onrows?: (rows: Row[]) => void } = $props();

  let rows = $state<Row[]>([]);
  let loaded = $state(false);
  let error = $state('');
  let cat = $state('');
  let level = $state('');
  let open = $state<string | null>(null);
  let revealAll = $state(sessionStorage.getItem('reveal-sensitive') === '1');
  let revealed = $state<Set<string>>(new Set());

  const RANK = "CASE overall WHEN 'Strong' THEN 0 WHEN 'Moderate' THEN 1 WHEN 'Limited' THEN 2 ELSE 3 END";
  $effect(() => {
    const sec = section;
    const gs = groups;
    loaded = false;
    (async () => {
      try {
        await useRelease(m);
        const r = await query(`SELECT * FROM claims WHERE section = ${sqlString(sec)} ORDER BY ${RANK}, clinvar_stars DESC, subject`);
        // older releases have no "group" column: keep everything
        rows = gs.length ? r.rows.filter((x) => x.group == null || gs.includes(x.group)) : r.rows;
        onrows?.(rows);
      } catch (e) { error = String(e); }
      loaded = true;
    })();
  });

  const cats = $derived([...new Set(rows.map((r) => r.category_label))]);
  const shown = $derived(rows.filter((r) => (!cat || r.category_label === cat) && (!level || r.overall === level)));
  const counts = $derived(Object.fromEntries(['Strong', 'Moderate', 'Limited', 'Not callable'].map((l) => [l, rows.filter((r) => r.overall === l).length])));
  const nSensitive = $derived(rows.filter((r) => r.sensitive).length);
  const hidden = (r: Row) => r.sensitive && !revealAll && !revealed.has(r.claim_id);
  function reveal(r: Row) { revealed = new Set([...revealed, r.claim_id]); }
  function toggleAll() { revealAll = !revealAll; sessionStorage.setItem('reveal-sensitive', revealAll ? '1' : '0'); }
  const json = (s: string | null | undefined, d: any = []) => { try { return s ? JSON.parse(s) : d; } catch { return d; } };
  let allLabs = $state(false);
  let showAll = $state(false);
  const visible = $derived(limit && !showAll ? shown.slice(0, limit) : shown);
  const copiesWord = (r: Row) => (r.zygosity === 'hom' ? (r.chrom === 'X' || r.chrom === 'Y' || r.chrom === 'MT' ? 'your only copy' : 'both copies') : 'one copy');
</script>

{#if !loaded}
  <p class="muted small">Loading findings…</p>
{:else if error}
  <p class="small warn">{error}</p>
{:else if !rows.length}
  <p class="muted small">{empty}</p>
{:else}
  {#if intro}<p class="small muted">{intro}</p>{/if}
  {#if filters}
  <div class="summary row small">
    {#each Object.entries(counts) as [l, n]}
      {#if n}<button class="chip" class:on={level === l} onclick={() => (level = level === l ? '' : l)}><span class="pill {OVERALL_PILL[l]}">{l}</span> {n}</button>{/if}
    {/each}
    <span class="spacer"></span>
    {#if cats.length > 1}
      <select bind:value={cat} aria-label="Category"><option value="">All categories ({rows.length})</option>{#each cats as c}<option value={c}>{c} ({rows.filter((r) => r.category_label === c).length})</option>{/each}</select>
    {/if}
  </div>
  {/if}
  {#if nSensitive}
    <div class="sensitive small">
      <span>🔒 {nSensitive} finding{nSensitive > 1 ? 's' : ''} here relate{nSensitive > 1 ? '' : 's'} to serious conditions and {nSensitive > 1 ? 'are' : 'is'} hidden until you choose to look. Most such findings in healthy people turn out to be carrier status, low-penetrance, or misclassified — and none is a diagnosis.</span>
      <button onclick={toggleAll}>{revealAll ? 'Hide again' : 'Reveal all'}</button>
    </div>
  {/if}
  <div class="tablewrap">
    <table class="small">
      <thead><tr><th>Gene / variant</th><th>What the evidence says</th><th>You carry</th><th><Term t="Evidence grade" label="Evidence" /></th><th><Term t="Call confidence" label="Your call" /></th><th>Overall</th></tr></thead>
      <tbody>
        {#each visible as r (r.claim_id)}
          {#if hidden(r)}
            <tr class="masked"><td colspan="6"><button class="link" onclick={() => reveal(r)}>🔒 Sensitive finding — click to reveal</button></td></tr>
          {:else}
            <tr class:open={open === r.claim_id} onclick={() => (open = open === r.claim_id ? null : r.claim_id)}>
              <td><strong>{r.gene ?? r.subject}</strong><div class="mono faint">{r.rsid ?? `${r.chrom}:${r.pos}`}</div></td>
              <td><span class="pill {CLASS_PILL[r.category] ?? ''}">{r.category_label}</span>{#if r.category !== 'predicted'} <span class="stars" title={STAR_WORD[r.clinvar_stars]}>{stars(r.clinvar_stars)}</span>{/if}
                {#if r.acmg_sf}<span class="pill acmg" title="On the ACMG secondary-findings list: genes where a finding is medically actionable">ACMG SF</span>{/if}
                <div class="cond">{r.conditions || '—'}</div></td>
              <td><span class="mono">{r.genotype}</span><div class="faint">{copiesWord(r)}</div>
                {#if r.role && ROLE[r.role] && r.section !== 'pgx' && r.section !== 'traits'}<span class="pill {ROLE[r.role][1]}" title={r.role_reason ?? ''}>{ROLE[r.role][0]}</span>{/if}
                {#if r.alt_is_major}<span class="pill" title="The reference genome carries the rarer allele here; the listed 'variant' is what most people have. See Reference genome in the primer.">majority allele</span>{/if}</td>
              <td><span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level}</span></td>
              <td><span class="pill {CALL_PILL[r.call_level]}">{r.call_level}</span></td>
              <td><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></td>
            </tr>
            {#if open === r.claim_id}
              <tr class="detail"><td colspan="6">
                <p class="statement">{r.statement}</p>
                {#if json(r.conditions_detail).length}
                  <div class="faint">The condition{json(r.conditions_detail).length > 1 ? 's' : ''}</div>
                  <div class="conds">
                    {#each json(r.conditions_detail) as c}
                      {@const cr = json(r.condition_roles).find((x: any) => x.name === c.name)}
                      <div class="condcard">
                        <div class="row"><strong>{c.name}</strong>
                          {#each c.modes ?? [] as md}{#if INH[md]}<span class="pill"><Term t={INH[md].startsWith("X-linked") ? "X-linked" : INH[md]} label={INH[md]} learn="Inheritance patterns" /></span>{/if}{/each}
                          {#if cr && ROLE[cr.role]}<span class="pill {ROLE[cr.role][1]}">for you: {ROLE[cr.role][0].toLowerCase()}</span>{/if}</div>
                        {#if c.definition}<p class="small">{c.definition.length > 320 ? c.definition.slice(0, 320) + '…' : c.definition}</p>{/if}
                        <div class="small faint">
                          {#if c.onset?.length}Typical onset: {c.onset.join(', ').toLowerCase()} · {/if}
                          {#if c.prevalence}How common: {c.prevalence} · {/if}
                          {#if c.modes_from}inheritance from {c.modes_from === 'name' ? 'the condition name' : c.modes_from}{/if}
                          {#each (c.ids ?? []).filter((x: string) => x.startsWith('MONDO:')).slice(0, 1) as id}
                            · <a href={`https://monarchinitiative.org/${id}`} target="_blank" rel="noreferrer">{id}</a>{/each}
                        </div>
                      </div>
                    {/each}
                  </div>
                  {#if r.conditions_also_listed}<p class="small faint">Also named on the ClinVar record (by fewer labs, or only in submissions listing every disease of the gene): {r.conditions_also_listed}</p>{/if}
                {/if}
                {#if r.acmg_sf || r.actionability}
                  {@const act = json(r.actionability, null)}
                  <div class="callout small">
                    {#if r.acmg_sf}<p><strong><Term t="ACMG SF" label="ACMG secondary-findings gene" /></strong> ({r.acmg_sf}). Labs that sequence genomes report P/LP variants in this gene even when not asked{r.acmg_rule && r.acmg_rule !== 'Any P/LP' ? ` — rule for this gene: ${r.acmg_rule}` : ''}.
                      {r.acmg_reportable ? 'This finding meets that rule, so a clinical lab would report it.' : 'This finding does not meet that rule, so a clinical lab would not report it as a secondary finding.'}</p>{/if}
                    {#if act}<p><strong><Term t="Actionability" label="ClinGen actionability" /></strong>: {act.score}/12 ({act.context.toLowerCase()} context) — {act.intervention?.toLowerCase()} to address {act.outcome?.toLowerCase()}. <a href={act.url} target="_blank" rel="noreferrer">Report</a></p>{/if}
                  </div>
                {/if}
                {#if r.am_score != null || r.revel != null || r.loeuf != null}
                  <div class="scores small">
                    <span class="faint">Computer predictions:</span>
                    {#if r.am_score != null}<span><Term t="AlphaMissense" />: <strong>{r.am_score.toFixed(2)}</strong> <span class="faint">({r.am_score > 0.564 ? 'likely pathogenic' : r.am_score < 0.34 ? 'likely benign' : 'ambiguous'})</span></span>{/if}
                    {#if r.revel != null}<span><Term t="REVEL" />: <strong>{r.revel.toFixed(2)}</strong> <span class="faint">({r.revel >= 0.932 ? 'strong' : r.revel >= 0.773 ? 'moderate' : r.revel >= 0.644 ? 'supporting' : r.revel <= 0.29 ? 'leans benign' : 'not informative'})</span></span>{/if}
                    {#if r.loeuf != null}<span><Term t="LOEUF" />: <strong>{r.loeuf.toFixed(2)}</strong> <span class="faint">({r.loeuf < 0.6 ? 'gene can’t tolerate broken copies' : 'gene tolerates broken copies'})</span></span>{/if}
                  </div>
                {/if}
                <div class="grid cols-3">
                  <div><div class="faint">Why this evidence grade</div><ul>{#each r.evidence_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
                  <div><div class="faint">Why this call confidence</div><ul>{#each r.call_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
                  <div><div class="faint">About the variant</div><ul>
                    <li>{consequenceLabel(r.consequence)}{r.aa_change ? ` (${aaShort(r.aa_change)})` : ''}{r.impact ? ` · ${r.impact.toLowerCase()} impact` : ''}</li>
                    {#if r.inheritance}<li>Inheritance: {r.inheritance.split('/').map((x: string) => INH[x] ?? x).join(', ')}{r.inheritance_basis === 'gene' ? ' (from the gene’s known diseases)' : ''}{r.gene_validity ? ` · gene–disease link: ${r.gene_validity}` : ''}</li>{/if}
                    <li>Frequency: {r.af_1kg != null ? `${pct(r.af_1kg)} of chromosomes (1000 Genomes)` : 'not in 1000 Genomes'}{r.gnomad_popmax_af != null ? ` · up to ${pct(r.gnomad_popmax_af)} in a gnomAD population` : ''}</li>
                  </ul></div>
                </div>
                {#if json(r.submitters).length}
                  {@const subs = json(r.submitters)}
                  <details class="labs small" bind:open={allLabs}>
                    <summary>Who says so: {subs.length}{subs.length >= 25 ? '+' : ''} lab submission{subs.length > 1 ? 's' : ''} to ClinVar
                      ({Object.entries(subs.reduce((a: any, x: any) => ((a[x.class] = (a[x.class] ?? 0) + 1), a), {})).map(([k, n]) => `${n} ${(CLASS_LABEL[k] ?? k).toLowerCase()}`).join(', ')})</summary>
                    <table><thead><tr><th>Lab</th><th>Verdict</th><th>For</th><th>Evaluated</th></tr></thead><tbody>
                      {#each subs as x}<tr><td>{x.submitter}</td><td><span class="pill {CLASS_PILL[x.class] ?? ''}">{x.significance}</span></td><td class="faint">{x.conditions.join('; ')}</td><td class="mono faint">{x.evaluated ?? '—'}</td></tr>{/each}
                    </tbody></table>
                  </details>
                {/if}
                <div class="row">
                  <span class="faint">Sources:</span>
                  {#each parseSources(r) as s}
                    {#if s.url}<a href={s.url} target="_blank" rel="noreferrer">{s.source} {s.version}</a>{:else}<span>{s.source} {s.version}</span>{/if}
                  {/each}
                  <span class="spacer"></span>
                  <button onclick={() => go('variants', r.rsid ?? `${r.chrom}:${r.pos}`)}>🔎 Open in variant explorer</button>
                </div>
              </td></tr>
            {/if}
          {/if}
        {/each}
      </tbody>
    </table>
  </div>
  {#if visible.length < shown.length}<button class="more" onclick={() => (showAll = true)}>Show all {shown.length} (sorted strongest first) ↓</button>{/if}
{/if}

<style>
  .more { margin-top: 8px; width: 100%; }
  .summary { margin: 6px 0 10px; }
  .chip { display: inline-flex; gap: 6px; align-items: center; padding: 3px 10px; }
  .chip.on { border-color: var(--accent); background: color-mix(in srgb, var(--accent) 12%, transparent); }
  .sensitive { display: flex; gap: 12px; align-items: center; background: color-mix(in srgb, var(--warn) 9%, transparent); border: 1px solid color-mix(in srgb, var(--warn) 35%, transparent); border-radius: 10px; padding: 8px 12px; margin-bottom: 10px; }
  .tablewrap { overflow-x: auto; }
  tbody tr:not(.detail):not(.masked) { cursor: pointer; }
  tbody tr:not(.detail):not(.masked):hover { background: var(--panel-2); }
  tr.open { background: color-mix(in srgb, var(--accent) 10%, transparent); }
  tr.detail td { background: var(--panel-2); padding: 12px 16px; }
  tr.detail ul { margin: 4px 0 8px; padding-left: 18px; }
  .statement { font-size: 0.95rem; max-width: 90ch; }
  .cond { max-width: 52ch; color: var(--muted); overflow: hidden; text-overflow: ellipsis; display: -webkit-box; -webkit-line-clamp: 2; line-clamp: 2; -webkit-box-orient: vertical; }
  .stars { color: var(--warn); letter-spacing: 1px; font-size: 0.8rem; }
  .masked td { color: var(--muted); font-style: italic; }
  .link { background: none; border: 0; padding: 0; color: var(--accent); }
  .warn { color: var(--warn); }
  .conds { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 8px; margin: 4px 0 8px; }
  .condcard { border: 1px solid var(--border); border-radius: 10px; padding: 8px 10px; background: var(--panel); }
  .condcard p { margin: 4px 0; }
  .callout { border-left: 3px solid var(--accent); padding: 4px 10px; margin: 8px 0; background: color-mix(in srgb, var(--accent) 6%, transparent); }
  .callout p { margin: 4px 0; }
  .scores { display: flex; flex-wrap: wrap; gap: 14px; margin: 6px 0 8px; }
  .labs { margin: 6px 0 10px; }
  .labs summary { cursor: pointer; color: var(--accent); }
  .labs table { margin-top: 6px; }
  .labs td { vertical-align: top; }
  .pill.acmg { color: var(--accent); background: color-mix(in srgb, var(--accent) 14%, transparent); }
</style>
