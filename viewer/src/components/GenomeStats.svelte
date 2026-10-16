<script lang="ts">
  // "How your genome compares": counts set against the 1000 Genomes Project's typical genome, only where the
  // definitions are close enough to compare (ADR-021). Each row quotes the paper.
  import type { Manifest } from '../lib/data';
  import { loadDoc } from '../lib/data';

  let { m }: { m: Manifest } = $props();
  const load = $derived(loadDoc(m, 'genome_stats'));
  const fmt = (n: number) => (n >= 1e6 ? `${(n / 1e6).toFixed(2)} million` : n.toLocaleString());
</script>

{#await load then st}
  {#if st?.rows?.length}
    <div class="rows">
      {#each st.rows as r (r.id)}
        <div class="r">
          <div class="v">{fmt(r.value)}</div>
          <div><strong>{r.label}</strong>{#if r.homozygous != null} <span class="faint">({r.homozygous} on both copies)</span>{/if}
            <div class="small">Typical genome: <strong>{r.typical}</strong>{#if r.in_range === true} <span class="pill pass">within the typical range</span>{:else if r.in_range === false} <span class="pill partial">outside the typical range</span>{/if}</div>
            <div class="small faint">{r.definition}.</div>
            <blockquote class="small">“{r.quote}” — <a href={r.source.fulltext ?? r.source.url} target="_blank" rel="noreferrer">{r.source.source}</a></blockquote>
          </div>
        </div>
      {/each}
    </div>
    <p class="small faint">{st.caveat}</p>
  {:else}
    <p class="muted small">Not computed for this release.</p>
  {/if}
{/await}

<style>
  .rows { display: grid; gap: 10px; }
  .r { display: grid; grid-template-columns: 130px 1fr; gap: 12px; align-items: start; }
  .v { font-size: 1.35rem; font-weight: 700; text-align: right; }
  blockquote { margin: 4px 0; padding-left: 8px; border-left: 2px solid var(--border); font-style: italic; color: var(--faint); }
  @media (max-width: 700px) { .r { grid-template-columns: 1fr; } .v { text-align: left; } }
</style>
