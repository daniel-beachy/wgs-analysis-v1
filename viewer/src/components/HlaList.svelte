<script lang="ts">
  // HLA types with well-established disease associations (ADR-021): what you carry, which way it points, and the
  // exact sentence from each study. Rows come from the curated, quoted table; statements are built by rule.
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
      if (m.tables.hla_associations) rows = (await query('SELECT * FROM hla_associations ORDER BY notable DESC, condition, direction')).rows;
      loaded = true;
    })();
  });

  const DIR: Record<string, [string, string]> = { protective: ['Lower risk', 'pass'], risk: ['Higher risk', 'partial'] };
  const quotes = (r: Row): any[] => { try { return JSON.parse(r.quotes ?? '[]'); } catch { return []; } };
</script>

{#if !loaded}
  <p class="muted small">Loading…</p>
{:else if !rows.length}
  <p class="muted small">No HLA types in this release (they come from T1K on the aligned reads).</p>
{:else}
  <div class="list">
    {#each rows as r (r.association_id)}
      {@const d = r.kind === 'rule_out' && r.notable ? ['Rules it out (mostly)', 'pass'] : DIR[r.direction] ?? ['', '']}
      <div class="item" class:dim={!r.notable}>
        <div class="head"><strong>{r.condition}</strong>
          <span class="pills">{#if r.notable}<span class="pill {d[1]}">{d[0]}</span>{:else}<span class="pill soon">Not carried</span>{/if}
            <span class="pill {OVERALL_PILL[r.overall]}" title="Overall: the weaker of evidence and call">{r.overall}</span></span></div>
        <p>{r.statement}.</p>
        {#each quotes(r) as q}
          <blockquote>“{q.quote}” <span class="cite">— <a href={q.url} target="_blank" rel="noreferrer">{q.source}</a></span></blockquote>
        {/each}
        <div class="small faint"><Term t="Evidence grade" label="Evidence" /> {r.evidence_level} ({(r.evidence_reasons ?? []).join('; ')}) · <Term t="Call confidence" label="Genotype call" /> <span class="pill {CALL_PILL[r.call_level]}">{r.call_level}</span> {(r.call_reasons ?? []).join('; ')}</div>
      </div>
    {/each}
  </div>
  <p class="small faint">An association is about groups of people, not a forecast. Odds ratios compare carriers with non-carriers; when a condition is rare, even a large odds ratio leaves most carriers unaffected (each card quotes how common the condition is where the study gives it). HLA types are read with T1K from your reads; a clinical HLA test is the way to confirm one before any medical decision.</p>
{/if}

<style>
  .list { display: grid; gap: 10px; }
  .item { border: 1px solid var(--border); border-radius: 10px; padding: 10px 12px; background: var(--panel-2, var(--panel)); }
  .item.dim { opacity: 0.7; }
  .head { display: flex; gap: 8px; align-items: baseline; flex-wrap: wrap; }
  .pills { margin-left: auto; display: flex; gap: 6px; }
  p { margin: 6px 0; }
  blockquote { margin: 6px 0 6px 4px; padding-left: 10px; border-left: 2px solid var(--border); font-style: italic; color: var(--muted, var(--faint)); font-size: 0.9rem; }
  .cite { font-style: normal; white-space: nowrap; }
</style>
