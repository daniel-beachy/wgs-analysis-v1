<script lang="ts">
  // Traits: well-replicated single variants first (easy to check against yourself), then polygenic scores.
  import type { Manifest } from '../lib/data';
  import type { SectionDef } from '../lib/sections';
  import { go } from '../lib/router.svelte';
  import { sectionFor } from '../lib/learn';
  import Card from '../components/Card.svelte';
  import Explain from '../components/Explain.svelte';
  import Term from '../components/Term.svelte';
  import TraitSnps from '../components/TraitSnps.svelte';
  import PgsList from '../components/PgsList.svelte';

  let { m, def, sex }: { m: Manifest; def: SectionDef; sex?: string } = $props();
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);
</script>

<div class="grid">
  <div class="banner small"><span class="bi">{def.icon}</span>
    <div><strong>{def.title}</strong> — {def.blurb}</div>
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}
  </div>

  <Explain title="How to read this page (start here)" learn="Polygenic scores">
    <p>Most traits are shaped by thousands of variants plus environment, so genetics gives a <em>tendency</em>, not an answer. This page has two kinds of result:</p>
    <ol>
      <li><strong>Single variants</strong> with a large, well-replicated effect (eye colour, milk digestion, alcohol flush). You can often check these against yourself — a fun test of the data.</li>
      <li><strong><Term t="PRS / PGS (Polygenic risk score)" label="Polygenic scores" /></strong>: one number adding up thousands to millions of small effects, shown as your <Term t="Percentile" label="percentile" /> among reference people whose DNA is most like yours. The 50th percentile is typical; most people fall between the 10th and 90th.</li>
    </ol>
    <p>The pills read as elsewhere: <strong>Evidence</strong> (how well the result predicts the trait in published studies) · <strong>Your call</strong> (how reliably it was measured in your data) · <strong>Overall</strong> (the weaker of the two). Polygenic scores are never graded Strong: even good ones explain only a slice of the variation between people.</p>
  </Explain>

  <Card title="1 · Single-variant traits" subtitle="Variants the GWAS Catalog links to a trait in several studies. Click a row for the study and grades.">
    <TraitSnps {m} />
  </Card>

  <Card title="2 · Polygenic scores" subtitle="Your percentile for body traits, against the 1000 Genomes reference people most similar to you. Click a score for details.">
    <PgsList {m} section="traits" {sex} />
  </Card>

  <Card title="3 · Browse more scores" subtitle="Every other trait in the PGS Catalog with a score that independent researchers have tested — picked by the same rule, not by hand.">
    <details><summary class="small">Show the list</summary>
      <PgsList {m} section="" {sex} featured={false} empty="No other independently tested scores in this release." />
    </details>
  </Card>
</div>

<style>
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
</style>
