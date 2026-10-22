<script lang="ts">
  // Claude's take (ADR-026): an AI-written reading of the highlights. The pipeline has already split it into blocks
  // and checked every citation; this only renders text segments and citation chips (never raw HTML).
  import type { Manifest } from '../lib/data';
  import { loadDoc } from '../lib/data';
  import { go } from '../lib/router.svelte';
  import { OVERALL_PILL } from '../lib/evidence';
  import Card from './Card.svelte';

  let { m }: { m: Manifest } = $props();
  const load = $derived(loadDoc<any>(m, 'take'));
  const cites = (doc: any) => Object.fromEntries((doc.citations ?? []).map((c: any) => [c.ref, c]));
  function open(c: any) {
    if (!c?.found || !c.page) return;
    if (c.table === 'highlights' && m.tables.highlights) {
      document.getElementById(`hl-${c.key}`)?.scrollIntoView({ behavior: 'smooth', block: 'center' });
    } else go(c.page);
  }
</script>

{#snippet segs(list: any[], cm: Record<string, any>)}
  {#each list as s}
    {#if s.t === 'b'}<strong>{s.v}</strong>{:else if s.t === 'i'}<em>{s.v}</em>{:else if s.t === 'cite'}
      {@const c = cm[s.ref]}
      <button class="cite" class:missing={!c?.found} class:changed={c?.changed} onclick={() => open(c)}
        title={c?.found ? `${c.title} — ${c.overall ?? 'not graded'}${c.changed ? ' (changed since this was written)' : ''}` : `${s.ref} is not in this release`}>
        {#if c?.found}{c.title}{#if c.overall} <span class="pill {OVERALL_PILL[c.overall]}">{c.overall}</span>{/if}{:else}⚠ {s.ref}{/if}</button>
    {:else}{s.v}{/if}
  {/each}
{/snippet}

{#await load then doc}
  <Card title="Claude’s take" subtitle={doc ? `Written by ${doc.author}${doc.written ? ` on ${doc.written}` : ''}${doc.written_for ? ` for release ${doc.written_for}` : ''}.` : ''}>
    {#snippet aside()}<span class="pill opinion" title="Not produced by a rule: an AI's reading of the results">Opinion</span>{/snippet}
    {#if !doc}
      <p class="small muted">No take has been written for this release. A take is an AI’s reading of the highlights, written in a Claude Code session and checked by the pipeline: every paragraph must cite results from this release.</p>
    {:else}
      {@const cm = cites(doc)}
      <p class="small note">This is an AI’s reading of your highlights, not a rule-based result. Every claim in it cites a result (the chips, with that result’s current grade); trust the result and its grade over the wording here.</p>
      {#if doc.stale}
        <p class="small warnbox">⚠ Written for an earlier release: {doc.citations.filter((c: any) => c.changed).length || 'some'} cited result{doc.citations.filter((c: any) => c.changed).length === 1 ? ' has' : 's have'} changed since. Treat it as out of date until it is rewritten.</p>
      {/if}
      {#if doc.problems?.length}
        <details class="small warnbox"><summary>⚠ {doc.problems.length} problem{doc.problems.length > 1 ? 's' : ''} found when checking the citations</summary>
          <ul>{#each doc.problems as p}<li>{p}</li>{/each}</ul></details>
      {/if}
      <div class="take">
        {#each doc.blocks as b}
          {#if b.type === 'h'}<h4>{b.text}</h4>
          {:else if b.type === 'ul'}<ul>{#each b.items as it}<li>{@render segs(it, cm)}</li>{/each}</ul>
          {:else}<p>{@render segs(b.segments, cm)}</p>{/if}
        {/each}
      </div>
    {/if}
  </Card>
{/await}

<style>
  .opinion { color: var(--accent-2); background: color-mix(in srgb, var(--accent-2) 14%, transparent); }
  .note { margin: 0 0 10px; color: var(--muted); border-left: 3px solid var(--accent-2); padding-left: 10px; }
  .warnbox { background: color-mix(in srgb, var(--warn) 9%, transparent); border: 1px solid color-mix(in srgb, var(--warn) 35%, transparent); border-radius: 10px; padding: 8px 12px; }
  .warnbox ul { margin: 6px 0 0; padding-left: 18px; }
  .take { max-width: 80ch; }
  .take h4 { margin: 14px 0 4px; }
  .take p, .take li { line-height: 1.6; }
  .take ul { padding-left: 20px; }
  .cite { display: inline; padding: 0 6px; margin: 0 1px; border-radius: 8px; font-size: 0.85em; line-height: 1.5;
    background: var(--chip); border: 1px solid var(--border); color: var(--text); white-space: normal; }
  .cite:hover { border-color: var(--accent); }
  .cite .pill { font-size: 0.72rem; padding: 0 6px; }
  .cite.changed { border-color: var(--warn); }
  .cite.missing { color: var(--warn); }
</style>
