<script lang="ts">
  // Repeat expansions at every STRchive disease locus (ADR-021): your repeat size next to the typical and
  // disease-causing ranges. Adult-onset conditions stay hidden until revealed, like other sensitive findings.
  import type { Manifest } from '../lib/data';
  import { query, useRelease, type Row } from '../lib/db';
  import { OVERALL_PILL, CALL_PILL } from '../lib/evidence';

  let { m }: { m: Manifest } = $props();
  let rows = $state<Row[]>([]);
  let loaded = $state(false);
  let showAll = $state(false);
  let revealAll = $state(sessionStorage.getItem('reveal-sensitive') === '1');

  $effect(() => {
    (async () => {
      await useRelease(m);
      if (m.tables.repeats) rows = (await query(`SELECT * FROM repeats ORDER BY
        CASE status WHEN 'disease_range' THEN 0 WHEN 'intermediate' THEN 1 WHEN 'not_readable' THEN 3 ELSE 2 END,
        CASE call_level WHEN 'Low' THEN 0 WHEN 'Medium' THEN 1 ELSE 2 END, gene`)).rows;
      loaded = true;
    })();
  });

  const STATUS: Record<string, [string, string]> = {
    typical: ['Typical', 'pass'], intermediate: ['Intermediate', 'partial'], disease_range: ['Disease range', 'fail'],
    not_readable: ['Not readable', 'warn'], no_thresholds: ['Not classified', 'soon'], disputed_locus: ['Disputed link', 'soon'],
  };
  const n = (s: string) => rows.filter((r) => r.status === s).length;
  const lowTypical = $derived(rows.filter((r) => r.status === 'typical' && r.call_level !== 'High').length);
  const flagged = $derived(rows.filter((r) => r.status === 'intermediate' || r.status === 'disease_range'));
  const shown = $derived(showAll ? rows : rows.filter((r) => r.status !== 'typical' || r.call_level !== 'High'));
  const nSensitive = $derived(rows.filter((r) => r.sensitive).length);
  function toggle() { revealAll = !revealAll; sessionStorage.setItem('reveal-sensitive', revealAll ? '1' : '0'); }
  const range = (r: Row) => [r.benign_max != null ? `typical ≤ ${r.benign_max}` : '', r.pathogenic_min != null ? `disease ≥ ${r.pathogenic_min}` : ''].filter(Boolean).join(' · ') || 'no published ranges';
</script>

{#if !loaded}
  <p class="muted small">Loading…</p>
{:else if !rows.length}
  <p class="muted small">Not checked in this release (needs the aligned reads and <code>pixi install -e special</code>).</p>
{:else}
  <div class="summary">
    <div class="verdict {flagged.length ? 'warn' : 'good'}">{flagged.length ? '!' : '✓'}</div>
    <div>
      <p class="big"><strong>{rows.length}</strong> repeat-expansion disease loci checked:
        <strong>{n('typical')}</strong> typical{lowTypical ? ` (${lowTypical} with low confidence)` : ''},
        <strong>{flagged.length}</strong> outside the typical range,
        <strong>{n('not_readable')}</strong> not readable{n('no_thresholds') + n('disputed_locus') ? `, ${n('no_thresholds') + n('disputed_locus')} not classified` : ''}.</p>
      <p class="small faint">“Not readable” means too few reads covered the repeat (often very GC-rich stretches such as fragile X, <em>FMR1</em>): unknown, not normal. A clinical repeat test is the way to check those.</p>
    </div>
  </div>
  {#if nSensitive}
    <div class="sensitive small">🔒 {nSensitive} of these loci are for adult-onset conditions (such as Huntington’s disease). Their names and results stay hidden until you choose to see them.
      <button onclick={toggle}>{revealAll ? 'Hide again' : 'Reveal'}</button></div>
  {/if}
  <table class="t small">
    <thead><tr><th>Gene · condition</th><th>Your repeat</th><th>Ranges (STRchive)</th><th>Result</th><th>Call</th><th>Overall</th></tr></thead>
    <tbody>
      {#each shown as r (r.locus_id)}
        {@const st = STATUS[r.status] ?? [r.status, '']}
        {#if r.sensitive && !revealAll}
          <tr class="masked"><td colspan="6">🔒 Adult-onset condition — hidden until you choose to reveal</td></tr>
        {:else}
          <tr>
            <td data-label="Gene"><strong>{r.gene}</strong> <span class="faint">{r.motif?.length <= 6 ? r.motif : 'VNTR'}</span><div class="faint"><a href={r.url} target="_blank" rel="noreferrer">{r.disease}</a></div></td>
            <td data-label="Your repeat">{#if r.status === 'not_readable'}<span class="faint">—</span>{:else}<span class="mono">{r.clinical_estimate ?? r.alleles}</span>{#if r.clinical_estimate && r.clinical_estimate !== r.alleles}<div class="faint">tool count {r.alleles}</div>{/if}<div class="faint">95% range {r.confidence_interval}</div>{/if}</td>
            <td data-label="Ranges" class="faint">{range(r)}</td>
            <td data-label="Result"><span class="pill {st[1]}">{st[0]}</span>{#if r.vcf_agrees}<div class="faint" title="DeepVariant found no length change here either">✓ VCF agrees</div>{/if}</td>
            <td data-label="Call"><span class="pill {CALL_PILL[r.call_level]}">{r.call_level}</span></td>
            <td data-label="Overall"><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></td>
          </tr>
        {/if}
      {/each}
    </tbody>
  </table>
  <button class="link small" onclick={() => (showAll = !showAll)}>{showAll ? 'Show only results needing a note' : `Show all ${rows.length} loci (including ${rows.length - shown.length} clearly typical) →`}</button>
  <p class="small faint">Repeat sizes come from ExpansionHunter (Dolzhenko 2019) on your reads, compared with STRchive’s published ranges. Counts are put on the clinical scale relative to the reference genome, because the tool’s counted region can include interruptions that clinical counts leave out. In a large UK study, this kind of genome screening found 97% of repeat expansions with 99.6% specificity (Ibañez 2022), but sizes are estimates: any result outside the typical range needs a clinical repeat test.</p>
{/if}

<style>
  .summary { display: flex; gap: 14px; align-items: center; margin-bottom: 10px; }
  .summary p { margin: 2px 0; }
  .big { font-size: 1.02rem; }
  .verdict { flex: none; width: 44px; height: 44px; border-radius: 50%; display: grid; place-items: center; font-size: 1.4rem; font-weight: 700; }
  .verdict.good { color: var(--good); background: color-mix(in srgb, var(--good) 16%, transparent); }
  .verdict.warn { color: var(--warn); background: color-mix(in srgb, var(--warn) 16%, transparent); }
  .sensitive { display: flex; gap: 12px; align-items: center; background: color-mix(in srgb, var(--warn) 9%, transparent); border: 1px solid color-mix(in srgb, var(--warn) 35%, transparent); border-radius: 10px; padding: 8px 12px; margin-bottom: 10px; }
  .t { width: 100%; border-collapse: collapse; }
  .t th { text-align: left; color: var(--faint); font-weight: 600; padding: 4px 6px; border-bottom: 1px solid var(--border); }
  .t td { padding: 6px; border-bottom: 1px solid var(--border); vertical-align: top; }
  .masked td { color: var(--faint); }
  .link { background: none; border: 0; padding: 6px 0; color: var(--accent); }
  @media (max-width: 700px) { .t thead { display: none; } .t, .t tbody, .t tr, .t td { display: block; } .t td { border: 0; padding: 2px 6px; } .t tr { border-bottom: 1px solid var(--border); padding: 4px 0; } .t td[data-label]::before { content: attr(data-label) ': '; color: var(--faint); } }
</style>
