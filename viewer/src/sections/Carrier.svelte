<script lang="ts">
  // Carrier: recessive conditions where you have one altered copy. Explains inheritance with a Punnett square,
  // then lists findings grouped by how solid they are.
  import type { Manifest } from '../lib/data';
  import type { SectionDef } from '../lib/sections';
  import type { Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { sectionFor } from '../lib/learn';
  import Card from '../components/Card.svelte';
  import ClaimList from '../components/ClaimList.svelte';
  import Explain from '../components/Explain.svelte';
  import Term from '../components/Term.svelte';

  let { m, def }: { m: Manifest; def: SectionDef } = $props();
  let rows = $state<Row[]>([]);
  const solid = $derived(rows.filter((r) => ['pathogenic', 'likely_pathogenic'].includes(r.category)));
  const weak = $derived(rows.filter((r) => !['pathogenic', 'likely_pathogenic'].includes(r.category)));
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);
  // Punnett square cells: [father allele, mother allele]
  const CELLS: [string, string][] = [['N', 'N'], ['c', 'N'], ['N', 'c'], ['c', 'c']];
  const OUTCOME = (a: string, b: string) => (a === 'c' && b === 'c' ? ['Affected', 'bad'] : a === 'c' || b === 'c' ? ['Carrier', 'mid'] : ['Unaffected', 'ok']);
</script>

<div class="grid">
  <div class="banner small"><span class="bi">{def.icon}</span>
    <div><strong>{def.title}</strong> — conditions that need <em>two</em> altered copies of a gene, where you have one. Carriers are almost always healthy; this matters mainly for family planning.</div>
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}
  </div>

  <div class="two">
    <Card title="What being a carrier means" subtitle="Everyone carries a few recessive variants — most people never find out">
      <div class="punnett">
        <svg viewBox="0 0 300 220" role="img" aria-label="Punnett square: two carrier parents have a 25% chance per child of an affected child">
          <text x="190" y="18" class="lab" text-anchor="middle">Parent 1 (carrier)</text>
          <text x="145" y="44" class="al" text-anchor="middle">N</text><text x="235" y="44" class="al c" text-anchor="middle">c</text>
          <text x="22" y="140" class="lab" text-anchor="middle" transform="rotate(-90 22 140)">Parent 2 (carrier)</text>
          <text x="66" y="100" class="al" text-anchor="middle">N</text><text x="66" y="190" class="al c" text-anchor="middle">c</text>
          {#each CELLS as [a, b], i}
            {@const x = 100 + (i % 2) * 90}
            {@const y = 56 + Math.floor(i / 2) * 82}
            {@const o = OUTCOME(a, b)}
            <rect {x} {y} width="86" height="78" rx="10" class="cell {o[1]}" />
            <text x={x + 43} y={y + 34} text-anchor="middle" class="gt">{[a, b].sort().join(' ')}</text>
            <text x={x + 43} y={y + 54} text-anchor="middle" class="out">{o[0]}</text>
            <text x={x + 43} y={y + 68} text-anchor="middle" class="out">1 in 4</text>
          {/each}
        </svg>
        <div class="small">
          <p><strong>N</strong> = working copy, <strong>c</strong> = altered copy. Each parent passes on one of their two copies at random.</p>
          <p>If <em>both</em> parents carry a variant in the <em>same</em> gene, each child has a <strong>1 in 4</strong> chance of the condition, 1 in 2 of being a carrier, 1 in 4 of neither.</p>
          <p>If only one parent is a carrier, children can be carriers but (for <Term t="Autosomal recessive" label="autosomal recessive" /> conditions) not affected.</p>
        </div>
      </div>
      <Explain title="X-linked conditions are different" learn="Inheritance patterns">
        <p>Genes on the X chromosome follow a different pattern. People with XY chromosomes have only one X, so a single altered copy can cause an <Term t="X-linked" label="X-linked" /> recessive condition (e.g. some colour blindness, haemophilia). People with XX are usually carriers. This dashboard handles that automatically, using the sex inferred from your coverage.</p>
      </Explain>
    </Card>

    <Card title="Your carrier summary">
      <div class="summary">
        <div class="stat"><div class="n">{solid.length}</div><div class="small muted">classified carrier results<br /><span class="faint">(lab-classified pathogenic / likely pathogenic)</span></div></div>
        <div class="stat"><div class="n faint">{weak.length}</div><div class="small muted">weaker possibilities<br /><span class="faint">(labs disagree, or computer prediction)</span></div></div>
      </div>
      <p class="small muted">The average person carries one to three recessive disease variants that a carrier screen would detect, so a handful here is expected.</p>
      <div class="small partner">
        <strong>👫 Planning a family?</strong> What matters is whether a partner carries a variant in the <em>same gene</em>. Clinical carrier screening for both partners (ACOG recommends it for anyone considering pregnancy) checks hundreds of genes properly, including ones short-read sequencing struggles with — such as <em>SMN1</em> (spinal muscular atrophy), <em>HBA1/2</em> (alpha thalassaemia) and <em>FMR1</em> (fragile X), which this dashboard does not assess yet.
      </div>
    </Card>
  </div>

  <Card title="Carrier findings" subtitle="One copy of a variant in a recessive gene. Click a row for the condition, who says so, and the reasons for each grade.">
    <ClaimList {m} section="carrier" onrows={(r) => (rows = r)} empty="No carrier findings: none of your variants is a known pathogenic variant in a recessive gene." />
  </Card>
</div>

<style>
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
  .two { display: grid; grid-template-columns: 3fr 2fr; gap: 16px; }
  @media (max-width: 900px) { .two { grid-template-columns: 1fr; } }
  .punnett { display: grid; grid-template-columns: minmax(220px, 300px) 1fr; gap: 16px; align-items: center; }
  @media (max-width: 600px) { .punnett { grid-template-columns: 1fr; } }
  .punnett p { margin: 0 0 8px; }
  svg { width: 100%; height: auto; }
  .lab { fill: var(--muted); font-size: 11px; }
  .al { fill: var(--text); font-size: 18px; font-weight: 700; }
  .al.c, .gt { fill: var(--text); }
  .gt { font-size: 18px; font-weight: 700; font-family: ui-monospace, monospace; }
  .out { fill: var(--muted); font-size: 11px; }
  .cell { stroke: var(--border); }
  .cell.ok { fill: color-mix(in srgb, var(--good) 14%, var(--panel)); }
  .cell.mid { fill: color-mix(in srgb, var(--warn) 16%, var(--panel)); }
  .cell.bad { fill: color-mix(in srgb, var(--bad, #e5484d) 22%, var(--panel)); }
  .summary { display: flex; gap: 24px; margin-bottom: 8px; }
  .stat { display: flex; gap: 10px; align-items: center; }
  .n { font-size: 2.2rem; font-weight: 700; }
  .partner { background: var(--panel-2); border: 1px solid var(--border); border-radius: 10px; padding: 10px 12px; }
</style>
