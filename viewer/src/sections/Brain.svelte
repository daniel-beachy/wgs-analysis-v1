<script lang="ts">
  // Brain & mind: neurological and psychiatric polygenic scores (with the caveats they need) and rare
  // protein-changing variants in genes SFARI links to brain development.
  import type { Manifest } from '../lib/data';
  import type { SectionDef } from '../lib/sections';
  import { query, useRelease, type Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { sectionFor } from '../lib/learn';
  import { OVERALL_PILL, consequenceLabel, aaShort, pct } from '../lib/evidence';
  import Card from '../components/Card.svelte';
  import Explain from '../components/Explain.svelte';
  import Term from '../components/Term.svelte';
  import PgsList from '../components/PgsList.svelte';

  let { m, def, sex }: { m: Manifest; def: SectionDef; sex?: string } = $props();
  let brain = $state<Row[]>([]);
  let showAll = $state(false);
  let open = $state<string | null>(null);
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);

  $effect(() => {
    (async () => {
      await useRelease(m);
      if (m.tables.brain_variants) brain = (await query(`SELECT * FROM brain_variants WHERE gene IS NOT NULL
        ORDER BY flagged DESC, sfari_score NULLS LAST, gene`)).rows;
    })();
  });
  // SFARI Gene's own scoring categories.
  const SFARI: Record<string, string> = { '1': 'high confidence', '2': 'strong candidate', '3': 'suggestive evidence' };
  const flagged = $derived(brain.filter((b) => b.flagged));
  const shown = $derived(showAll ? brain : flagged);
  const key = (b: Row) => `${b.chrom}:${b.pos}:${b.alt}`;
</script>

<div class="grid">
  <div class="banner small"><span class="bi">{def.icon}</span>
    <div><strong>{def.title}</strong> — {def.blurb}</div>
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}
  </div>

  <Explain title="Read this first: what these scores can and cannot say" learn="Polygenic scores">
    <p>Scores for brain-related conditions and for things like educational attainment are among the most debated in genetics. They are shown because you asked to see everything, with these caveats:</p>
    <ul>
      <li><strong>They predict weakly.</strong> Even the best explain a few percent of the differences between people; most people with a high score never develop the condition.</li>
      <li><strong>Environment is tangled in.</strong> Educational-attainment and cognition scores partly capture family background and opportunity, not just biology: comparisons between siblings, which remove family background, find much smaller effects (<a href="https://doi.org/10.1038/s41588-022-01062-7" target="_blank" rel="noreferrer">Howe et al., 2022</a>).</li>
      <li><strong>They were built mostly from people of European ancestry</strong> and work less well for others.</li>
      <li><strong>They are not a diagnosis</strong> and say nothing about who you are.</li>
    </ul>
  </Explain>

  <Card title="1 · Polygenic scores" subtitle="Your percentile among the 1000 Genomes reference people most similar to you. Click a score for what it predicts and how well.">
    <PgsList {m} section="brain" {sex} />
  </Card>

  <Card title="2 · Rare variants in brain-development genes" subtitle={`Rare, protein-changing variants in genes the SFARI Gene database links to autism and neurodevelopment${m.knowledge?.versions?.sfari ? ` (${m.knowledge.versions.sfari})` : ''}`}>
    <Explain title="What does this mean?" learn="From variant to meaning">
      <p>Everyone carries rare variants in these genes. Only those that prediction tools flag as likely damaging are listed by default, and they are always graded <strong>Limited</strong>: no lab has classified them, and many such variants turn out harmless. Variants already classified by labs appear in <button class="link" onclick={() => go('health')}>Health</button>.</p>
    </Explain>
    {#if !brain.length}<p class="muted small">No rare protein-changing variants in these genes in this release.</p>
    {:else}
      <p class="small"><strong>{flagged.length}</strong> flagged by prediction tools · {brain.length - flagged.length} more rare protein-changing variant{brain.length - flagged.length === 1 ? '' : 's'} that the tools consider unlikely to matter.</p>
      {#if shown.length}
      <div class="tablewrap">
        <table class="small">
          <thead><tr><th>Gene</th><th>Change</th><th>SFARI gene score</th><th>How rare</th><th>Evidence</th></tr></thead>
          <tbody>
            {#each shown as b (key(b))}
              <tr onclick={() => (open = open === key(b) ? null : key(b))}>
                <td><strong>{b.gene}</strong><div class="mono faint">{b.rsid ?? `${b.chrom}:${b.pos}`} · {b.zygosity === 'hom' ? 'both copies' : 'one copy'}</div></td>
                <td>{consequenceLabel(b.consequence)}{b.aa_change ? ` (${aaShort(b.aa_change)})` : ''}</td>
                <td>{b.sfari_score ? `${b.sfari_score} — ${SFARI[String(b.sfari_score)] ?? ''}` : '—'}{b.syndromic ? ' · syndromic' : ''}</td>
                <td title="Highest frequency in any gnomAD population">{b.gnomad_popmax_af != null ? pct(b.gnomad_popmax_af) : b.af_1kg != null ? pct(b.af_1kg) : 'not seen in reference data'}</td>
                <td>{#if b.evidence_level}<span class="pill {OVERALL_PILL[b.evidence_level]}">{b.evidence_level}</span>{:else}<span class="faint">not flagged</span>{/if}</td>
              </tr>
              {#if open === key(b)}
                <tr class="detail"><td colspan="5"><ul>
                  {#each b.evidence_reasons ?? [] as x}<li>{x}</li>{/each}
                  {#if b.revel != null}<li><Term t="REVEL" />: {b.revel.toFixed(2)}</li>{/if}
                  {#if b.am_score != null}<li><Term t="AlphaMissense" />: {b.am_score.toFixed(2)}</li>{/if}
                  {#if b.loeuf != null}<li><Term t="LOEUF" />: {b.loeuf.toFixed(2)}</li>{/if}
                  <li>SFARI: {b.category} · {b.reports} report{b.reports === 1 ? '' : 's'} · <a href={`https://gene.sfari.org/database/human-gene/${b.gene}`} target="_blank" rel="noreferrer">gene page</a></li>
                </ul></td></tr>
              {/if}
            {/each}
          </tbody>
        </table>
      </div>
      {/if}
      {#if brain.length > flagged.length}<button class="link small" onclick={() => (showAll = !showAll)}>{showAll ? 'Show flagged only' : 'Show all rare protein-changing variants →'}</button>{/if}
    {/if}
  </Card>
</div>

<style>
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
  .link { background: none; border: 0; padding: 4px 0; color: var(--accent); }
  tbody tr:not(.detail) { cursor: pointer; }
  .detail ul { margin: 4px 0; padding-left: 18px; }
</style>
