<script lang="ts">
  import type { SectionDef } from '../lib/sections';
  import { sectionFor } from '../lib/learn';
  import { go } from '../lib/router.svelte';
  import type { Manifest } from '../lib/data';
  import Card from '../components/Card.svelte';
  import ClaimList from '../components/ClaimList.svelte';
  let { def, status, m }: { def: SectionDef; status?: { status: string; detail: string; note: string }; m?: Manifest } = $props();
  const CLAIM_INTRO: Record<string, string> = {
    health: 'Your variants that ClinVar links to a disease or risk, graded by how strong that link is and how confidently your genotype was called. Most are risk factors with small effects; pathogenic entries are rarer and deserve clinical confirmation.',
    carrier: 'Variants in genes where disease needs two altered copies, and you have one. Usually no effect on you; relevant for family planning.',
    pgx: 'Variants ClinVar records as affecting a drug response. Full medicine-by-medicine guidance (CPIC/DPWG, star alleles) arrives in Step 5; until then, ClinVar entries are capped at Moderate.',
    traits: 'Variants ClinVar records as associated with a trait. Curated trait panels and polygenic scores arrive in Step 6; until then these are mostly Limited.',
  };
  const early = $derived(!!m?.tables.claims && def.module != null && def.module in CLAIM_INTRO);
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);
  // detail format from the capability model: "+ extra a; extra b | - missing kind (why); kind (why)"
  const parts = $derived.by(() => {
    const out = { plus: [] as string[], missing: [] as string[] };
    for (const seg of (status?.detail ?? '').split(' | ')) {
      const body = seg.slice(2).trim();
      if (seg.startsWith('+ ')) out.plus.push(...body.split(/;\s*/));
      else if (seg.startsWith('- missing ')) out.missing.push(...seg.slice(10).split(/;\s*(?![^(]*\))/).map((x) => x.replace(/^external:/, '').replace(/\(ingest an? /, '(').replace(/^(\S+) \((.+)\)$/, '$2 ($1)')));
    }
    return out;
  });
</script>

{#if early && m}
  <div class="grid">
    <div class="banner small"><span class="bi">{def.icon}</span><div><strong>{def.title}</strong> — {def.blurb}<br /><span class="faint">Early view: graded ClinVar findings. The full section arrives in step {def.step}{def.evidence ? `, built on ${def.evidence}` : ''}</span></div>
      {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}</div>
    <Card title="Graded findings" subtitle="From ClinVar, joined to your variants locally. Click a row for the reasons and sources.">
      <ClaimList {m} section={def.module!} intro={CLAIM_INTRO[def.module!]} />
    </Card>
  </div>
{:else}
<Card>
  <div class="up">
    <div class="icon">{def.icon}</div>
    <h2>{def.title}</h2>
    <p class="lead">{def.blurb}</p>
    {#if def.evidence}<p class="small"><strong>Evidence it will rest on:</strong> <span class="muted">{def.evidence}</span></p>{/if}
    <div class="row center">
      <span class="pill soon">Arrives in step {def.step} of the build plan</span>
      {#if status}<span class="pill {status.status}" title={status.detail}>your data: {status.status}</span>{/if}
    </div>
    {#if parts.plus.length || parts.missing.length || status?.note}
      <div class="detail small">
        {#if parts.plus.length}<p><strong>Your files also make possible:</strong></p><ul>{#each parts.plus as p}<li>{p}</li>{/each}</ul>{/if}
        {#if parts.missing.length}<p><strong>Optional extras not found</strong> <span class="faint">(results from Linux-only tools; ingested automatically if present):</span></p><ul class="faint">{#each parts.missing as p}<li>{p}</li>{/each}</ul>{/if}
        {#if status?.note}<p class="faint">{status.note}</p>{/if}
      </div>
    {/if}
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Read the background while you wait</button>{/if}
  </div>
</Card>
{/if}

<style>
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
  .up { text-align: center; padding: 30px 10px; max-width: 720px; margin: 0 auto; }
  .icon { font-size: 3rem; }
  .lead { margin: 0 auto 1em; }
  .center { justify-content: center; margin: 14px 0 8px; }
  .detail { max-width: 600px; margin: 6px auto 16px; text-align: left; background: var(--panel-2); border: 1px solid var(--border); border-radius: 10px; padding: 4px 16px; }
  .detail ul { margin: 0 0 8px; padding-left: 20px; }
  .detail p { margin: 8px 0 4px; }
</style>
