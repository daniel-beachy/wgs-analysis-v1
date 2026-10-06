<script lang="ts">
  // Graded claims for one dashboard section, each showing what is claimed, how strong the evidence is,
  // how confident your genotype call is, and the sources behind it. Sensitive findings stay hidden until revealed.
  import type { Manifest } from '../lib/data';
  import { query, useRelease, sqlString, type Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { OVERALL_PILL, CALL_PILL, CLASS_PILL, stars, STAR_WORD, pct, parseSources, consequenceLabel, aaShort } from '../lib/evidence';
  import Term from './Term.svelte';

  let { m, section, intro = '' }: { m: Manifest; section: string; intro?: string } = $props();

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
    loaded = false;
    (async () => {
      try {
        await useRelease(m);
        const r = await query(`SELECT * FROM claims WHERE section = ${sqlString(sec)} ORDER BY ${RANK}, clinvar_stars DESC, subject`);
        rows = r.rows;
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
  const copiesWord = (r: Row) => (r.zygosity === 'hom' ? (r.chrom === 'X' || r.chrom === 'Y' || r.chrom === 'MT' ? 'your only copy' : 'both copies') : 'one copy');
</script>

{#if !loaded}
  <p class="muted small">Loading findings…</p>
{:else if error}
  <p class="small warn">{error}</p>
{:else if !rows.length}
  <p class="muted small">No graded findings for this section in this release.</p>
{:else}
  {#if intro}<p class="small muted">{intro}</p>{/if}
  <div class="summary row small">
    {#each Object.entries(counts) as [l, n]}
      {#if n}<button class="chip" class:on={level === l} onclick={() => (level = level === l ? '' : l)}><span class="pill {OVERALL_PILL[l]}">{l}</span> {n}</button>{/if}
    {/each}
    <span class="spacer"></span>
    {#if cats.length > 1}
      <select bind:value={cat} aria-label="Category"><option value="">All categories ({rows.length})</option>{#each cats as c}<option value={c}>{c} ({rows.filter((r) => r.category_label === c).length})</option>{/each}</select>
    {/if}
  </div>
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
        {#each shown as r (r.claim_id)}
          {#if hidden(r)}
            <tr class="masked"><td colspan="6"><button class="link" onclick={() => reveal(r)}>🔒 Sensitive finding — click to reveal</button></td></tr>
          {:else}
            <tr class:open={open === r.claim_id} onclick={() => (open = open === r.claim_id ? null : r.claim_id)}>
              <td><strong>{r.gene ?? r.subject}</strong><div class="mono faint">{r.rsid ?? `${r.chrom}:${r.pos}`}</div></td>
              <td><span class="pill {CLASS_PILL[r.category] ?? ''}">{r.category_label}</span> <span class="stars" title={STAR_WORD[r.clinvar_stars]}>{stars(r.clinvar_stars)}</span>
                <div class="cond">{r.conditions || '—'}</div></td>
              <td><span class="mono">{r.genotype}</span><div class="faint">{copiesWord(r)}</div>
                {#if r.alt_is_major}<span class="pill" title="The reference genome carries the rarer allele here; the listed 'variant' is what most people have. See Reference genome in the primer.">majority allele</span>{/if}</td>
              <td><span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level}</span></td>
              <td><span class="pill {CALL_PILL[r.call_level]}">{r.call_level}</span></td>
              <td><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></td>
            </tr>
            {#if open === r.claim_id}
              <tr class="detail"><td colspan="6">
                <p class="statement">{r.statement}</p>
                <div class="grid cols-3">
                  <div><div class="faint">Why this evidence grade</div><ul>{#each r.evidence_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
                  <div><div class="faint">Why this call confidence</div><ul>{#each r.call_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
                  <div><div class="faint">About the variant</div><ul>
                    <li>{consequenceLabel(r.consequence)}{r.aa_change ? ` (${aaShort(r.aa_change)})` : ''}{r.impact ? ` · ${r.impact.toLowerCase()} impact` : ''}</li>
                    {#if r.inheritance}<li>Inheritance: {r.inheritance}{r.gene_validity ? ` · gene–disease link: ${r.gene_validity}` : ''}</li>{/if}
                    <li>Frequency: {r.af_1kg != null ? `${pct(r.af_1kg)} of chromosomes (1000 Genomes)` : 'not in 1000 Genomes'}{r.gnomad_popmax_af != null ? ` · up to ${pct(r.gnomad_popmax_af)} in a gnomAD population` : ''}</li>
                  </ul></div>
                </div>
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
{/if}

<style>
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
</style>
