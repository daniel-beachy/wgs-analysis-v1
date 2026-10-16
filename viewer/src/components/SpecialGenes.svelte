<script lang="ts">
  // Genes with a near-identical twin (SMN1, GBA1) that standard variant calling can't read, re-read from the aligned
  // reads by specialised callers (ADR-021). "Not determined" is shown as such, never as an all-clear.
  import type { Manifest } from '../lib/data';
  import { query, useRelease, type Row } from '../lib/db';
  import { OVERALL_PILL, CALL_PILL } from '../lib/evidence';
  import Term from './Term.svelte';

  let { m }: { m: Manifest } = $props();
  let rows = $state<Row[]>([]);
  let loaded = $state(false);

  $effect(() => {
    (async () => {
      await useRelease(m);
      if (m.tables.special_genes) rows = (await query('SELECT * FROM special_genes ORDER BY gene')).rows;
      loaded = true;
    })();
  });

  const STATUS: Record<string, [string, string]> = {
    typical: ['Not a carrier', 'pass'], not_determined: ['Not determined', 'warn'], carrier: ['Carrier', 'partial'],
    possible_silent_carrier: ['Possible silent carrier', 'partial'], affected: ['Disease genotype', 'fail'],
  };
  const sources = (r: Row): any[] => { try { return JSON.parse(r.sources ?? '[]'); } catch { return []; } };
</script>

{#if !loaded}
  <p class="muted small">Loading…</p>
{:else if !rows.length}
  <p class="muted small">Not checked in this release (needs the aligned reads and <code>pixi install -e special</code>).</p>
{:else}
  <div class="genes">
    {#each rows as r (r.gene)}
      {@const st = STATUS[r.status] ?? [r.status, '']}
      <div class="gene">
        <div class="head"><strong>{r.gene}</strong> <span class="faint small">{r.condition}</span>
          <span class="pills"><span class="pill {st[1]}">{st[0]}</span><span class="pill {OVERALL_PILL[r.overall]}" title="Overall: the weaker of evidence and call">{r.overall}</span></span></div>
        <p>{r.statement}.</p>
        {#if r.next_step}<p class="small callout"><strong>To find out:</strong> {r.next_step}.</p>{/if}
        <div class="grid cols-2 small">
          <div><div class="faint"><Term t="Evidence grade" label="Evidence" />: <span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level}</span></div><ul>{#each r.evidence_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
          <div><div class="faint"><Term t="Call confidence" label="Genotype call" />: <span class="pill {CALL_PILL[r.call_level]}">{r.call_level}</span></div><ul>{#each r.call_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
        </div>
        <div class="small faint">Read by {r.tool} from your aligned reads ·
          {#each sources(r) as s, i}{#if i}{' · '}{/if}<a href={s.url} target="_blank" rel="noreferrer">{s.source}</a>{/each}</div>
      </div>
    {/each}
  </div>
{/if}

<style>
  .genes { display: grid; gap: 10px; }
  .gene { border: 1px solid var(--border); border-radius: 10px; padding: 10px 12px; background: var(--panel-2, var(--panel)); }
  .head { display: flex; gap: 8px; align-items: baseline; flex-wrap: wrap; }
  .pills { margin-left: auto; display: flex; gap: 6px; }
  p { margin: 6px 0; }
  .callout { background: color-mix(in srgb, var(--warn) 9%, transparent); border-radius: 8px; padding: 6px 10px; }
  ul { margin: 4px 0; padding-left: 18px; }
</style>
