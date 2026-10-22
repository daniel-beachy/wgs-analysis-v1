<script lang="ts">
  // Highlights: the notable results from every section, strongest evidence first (ADR-025), then Claude's take,
  // an AI-written reading of them that is labelled as opinion (ADR-026).
  import type { Manifest } from '../lib/data';
  import type { SectionDef } from '../lib/sections';
  import { SECTION_DEFS } from '../lib/sections';
  import { query, useRelease, type Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { sectionFor } from '../lib/learn';
  import { getReveal, setReveal } from '../lib/reveal';
  import { OVERALL_PILL, CALL_PILL, parseSources } from '../lib/evidence';
  import Card from '../components/Card.svelte';
  import Term from '../components/Term.svelte';
  import AbsRisk from '../components/AbsRisk.svelte';
  import Take from '../components/Take.svelte';

  let { m, def }: { m: Manifest; def: SectionDef } = $props();
  let rows = $state<Row[]>([]);
  let loaded = $state(false);
  let error = $state('');
  let level = $state('');
  let revealAll = $state(getReveal());
  let revealed = $state<Set<string>>(new Set());
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);

  $effect(() => {
    loaded = false;
    (async () => {
      try {
        await useRelease(m);
        rows = (await query('SELECT * FROM highlights ORDER BY rank')).rows;
      } catch (e) { error = String(e); }
      loaded = true;
    })();
  });

  const counts = $derived(['Strong', 'Moderate', 'Limited'].map((l) => [l, rows.filter((r) => r.overall === l).length] as const));
  const shown = $derived(rows.filter((r) => !level || r.overall === level));
  const nSensitive = $derived(rows.filter((r) => r.sensitive).length);
  const hidden = (r: Row) => r.sensitive && !revealAll && !revealed.has(r.highlight_id);
  const page = (id: string) => SECTION_DEFS.find((s) => s.id === id);
  function toggleAll() { revealAll = !revealAll; setReveal(revealAll); }
</script>

<div class="grid">
  <div class="banner small"><span class="bi">{def.icon}</span>
    <div><strong>{def.title}</strong> — the results most worth knowing about, gathered from every section. Nothing here is new: each card is a result another section already graded, picked by a fixed rule, with its grade unchanged.</div>
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}
  </div>

  <Card title="Your highlights" subtitle="Strongest evidence first. Each card says why it is here; open its section for the full evidence.">
    {#if !loaded}
      <p class="muted small">Loading highlights…</p>
    {:else if error}
      <p class="small warn">{error}</p>
    {:else if !rows.length}
      <p class="muted small">Nothing in this release met a highlight rule. That usually means no single-gene, carrier, medicine or strongly supported trait result stood out.</p>
    {:else}
      <div class="row small chips">
        {#each counts as [l, n]}
          {#if n}<button class="chip" class:on={level === l} onclick={() => (level = level === l ? '' : l)}><span class="pill {OVERALL_PILL[l]}">{l}</span> {n}</button>{/if}
        {/each}
        <span class="spacer"></span>
        <span class="faint">{rows.length} highlights</span>
      </div>
      {#if nSensitive}
        <div class="sensitive small">
          <span>🔒 {nSensitive} highlight{nSensitive > 1 ? 's' : ''} relate{nSensitive > 1 ? '' : 's'} to serious conditions and {nSensitive > 1 ? 'are' : 'is'} hidden until you choose to look. None is a diagnosis.</span>
          <button onclick={toggleAll}>{revealAll ? 'Hide again' : 'Reveal all'}</button>
        </div>
      {/if}
      <div class="cards">
        {#each shown as r (r.highlight_id)}
          {#if hidden(r)}
            <button class="hl masked" onclick={() => (revealed = new Set([...revealed, r.highlight_id]))}>
              <span class="pill">{r.kind_label}</span>
              <span>🔒 Sensitive finding — click to reveal</span>
              <span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span>
            </button>
          {:else}
            <article class="hl" id={`hl-${r.highlight_id}`}>
              <div class="head">
                <span class="pill kind">{r.kind_label}</span>
                <span class="spacer"></span>
                <span class="pill {OVERALL_PILL[r.overall]}" title="Overall grade: the weaker of the two grades below">{r.overall}</span>
              </div>
              <h4>{r.title}</h4>
              <p class="statement">{r.statement}</p>
              <p class="small why"><span class="faint">Why it’s here:</span> {r.notable_because}</p>
              {#if r.kind === 'pgs' && r.extra}<AbsRisk a={r.extra} compact />{/if}
              <div class="grades small">
                <span><Term t="Evidence grade" label="Evidence" /> <span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level ?? '—'}</span></span>
                <span><Term t="Call confidence" label="Your call" /> <span class="pill {CALL_PILL[r.call_level]}">{r.call_level ?? '—'}</span></span>
              </div>
              <div class="foot small">
                <span class="srcs">
                  {#each parseSources(r).slice(0, 3) as s}
                    {#if s.url}<a href={s.url} target="_blank" rel="noreferrer">{s.source}{s.version ? ` ${s.version}` : ''}</a>{:else}<span>{s.source}</span>{/if}
                  {/each}
                </span>
                {#if page(r.section)}<button class="link" onclick={() => go(r.section)}>Open in {page(r.section)?.title} →</button>{/if}
              </div>
            </article>
          {/if}
        {/each}
      </div>
    {/if}
  </Card>

  <Take {m} />

  <Card title="How highlights are picked" subtitle="Fixed rules over the graded results, so the same data always gives the same highlights.">
    <ul class="small rules">
      <li><strong>Single-gene findings and carrier results</strong>: every one the Health and Carrier sections list.</li>
      <li><strong>Computer predictions</strong>: shown as hypotheses, always Limited.</li>
      <li><strong>Risk factors</strong>: only those graded at least Moderate; most common-variant risks are too small to stand out.</li>
      <li><strong>Medicine genes</strong>: any gene whose predicted effect isn’t normal, with how many medicines’ guideline advice differs.</li>
      <li><strong>HLA types</strong>: associations the HLA section marks as notable for you.</li>
      <li><strong>Polygenic scores</strong>: a featured trait where you are in the top or bottom 5% and the studies don’t disagree.</li>
      <li><strong>Traits</strong>: a well-supported single-variant trait where fewer than 1 in 4 people share your genotype.</li>
      <li><strong>Hard-to-read genes, repeats and deletions</strong>: anything those specialised scans flagged.</li>
    </ul>
    <p class="small faint">Variants of uncertain significance are left out: they are not answers yet. A highlight is never graded higher than the result it came from; a self-check enforces that.</p>
  </Card>
</div>

<style>
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
  @media (max-width: 600px) { .banner { flex-wrap: wrap; } }
  .chips { margin: 0 0 10px; }
  .chip { display: inline-flex; gap: 6px; align-items: center; padding: 3px 10px; }
  .chip.on { border-color: var(--accent); background: color-mix(in srgb, var(--accent) 12%, transparent); }
  .sensitive { display: flex; gap: 12px; align-items: center; background: color-mix(in srgb, var(--warn) 9%, transparent); border: 1px solid color-mix(in srgb, var(--warn) 35%, transparent); border-radius: 10px; padding: 8px 12px; margin-bottom: 10px; }
  .cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(min(100%, 320px), 1fr)); gap: 12px; }
  .hl { display: flex; flex-direction: column; gap: 6px; background: var(--panel-2); border: 1px solid var(--border); border-radius: 12px; padding: 12px 14px; min-width: 0; text-align: left; }
  .hl.masked { flex-direction: row; align-items: center; justify-content: space-between; color: var(--muted); font-style: italic; cursor: pointer; }
  .head { display: flex; align-items: center; gap: 6px; }
  .kind { color: var(--accent); background: color-mix(in srgb, var(--accent) 12%, transparent); }
  h4 { margin: 0; font-size: 1rem; overflow-wrap: anywhere; }
  .statement { margin: 0; font-size: 0.92rem; overflow-wrap: anywhere; }
  .why { margin: 0; }
  .grades { display: flex; gap: 14px; flex-wrap: wrap; }
  .foot { display: flex; gap: 8px; align-items: center; flex-wrap: wrap; margin-top: auto; padding-top: 4px; }
  .srcs { display: flex; gap: 8px; flex-wrap: wrap; flex: 1; min-width: 0; }
  .srcs a, .srcs span { overflow-wrap: anywhere; }
  .link { background: none; border: 0; padding: 0; color: var(--accent); }
  .rules { padding-left: 18px; margin: 0 0 8px; }
  .rules li { margin-bottom: 4px; }
</style>
