<script lang="ts">
  // Blood group (ABO + RhD) read from the genome (ADR-021), with the deciding evidence and quoted sources.
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
      if (m.tables.blood_group) rows = (await query("SELECT * FROM blood_group ORDER BY system = 'RhD'")).rows;
      loaded = true;
    })();
  });
  const abo = $derived(rows.find((r) => r.system === 'ABO'));
  const rh = $derived(rows.find((r) => r.system === 'RhD'));
  const sign = (r?: Row) => (!r?.result ? '?' : r.result === 'RhD-negative' ? '−' : '+');
  const quotes = (r: Row): any[] => { try { return JSON.parse(r.quotes ?? '[]'); } catch { return []; } };
</script>

{#if !loaded}
  <p class="muted small">Loading…</p>
{:else if !rows.length}
  <p class="muted small">Not typed in this release.</p>
{:else}
  <div class="top">
    <div class="badge">{abo?.result ?? '?'}<sup>{sign(rh)}</sup></div>
    <div><p class="big">Blood group <strong>{abo?.result ?? 'not determined'} {rh?.result ?? ''}</strong></p>
      <p class="small faint">A fun one to check: if you have ever given blood or had a blood test card, it should say the same.</p></div>
  </div>
  <div class="grid cols-2">
    {#each rows as r (r.system)}
      <div class="sys">
        <div class="head"><strong>{r.system}</strong><span class="pills"><span class="pill {CALL_PILL[r.call_level]}">Call {r.call_level}</span><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></span></div>
        <p>{r.statement}.</p>
        <ul class="small">{#each r.call_reasons ?? [] as x}<li>{x}</li>{/each}</ul>
        {#each quotes(r) as q}<blockquote class="small">“{q.quote}” — <a href={q.url} target="_blank" rel="noreferrer">{q.source}</a></blockquote>{/each}
      </div>
    {/each}
  </div>
  <p class="small faint"><Term t="Reference genome" label="Note" />: the reference genome itself carries an O version of ABO, so “no change from the reference” at that position means O. Genome-based typing matched lab typing in over 99% of cases in the study quoted above, but transfusions always use a lab test.</p>
{/if}

<style>
  .top { display: flex; gap: 14px; align-items: center; margin-bottom: 10px; }
  .badge { width: 64px; height: 64px; border-radius: 50%; display: grid; place-items: center; font-size: 1.6rem; font-weight: 800; background: color-mix(in srgb, var(--bad) 14%, transparent); color: var(--text); }
  .badge sup { font-size: 1rem; }
  .big { font-size: 1.1rem; margin: 0; }
  .sys { border: 1px solid var(--border); border-radius: 10px; padding: 10px 12px; }
  .head { display: flex; gap: 8px; align-items: baseline; }
  .pills { margin-left: auto; display: flex; gap: 6px; }
  p { margin: 6px 0; }
  ul { margin: 4px 0; padding-left: 18px; }
  blockquote { margin: 4px 0; padding-left: 8px; border-left: 2px solid var(--border); font-style: italic; color: var(--faint); }
</style>
