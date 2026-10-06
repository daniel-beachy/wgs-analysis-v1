<script lang="ts">
  import type { SectionDef } from '../lib/sections';
  import { sectionFor } from '../lib/learn';
  import { go } from '../lib/router.svelte';
  import Card from '../components/Card.svelte';
  let { def, status }: { def: SectionDef; status?: { status: string; detail: string; note: string } } = $props();
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

<style>
  .up { text-align: center; padding: 30px 10px; max-width: 720px; margin: 0 auto; }
  .icon { font-size: 3rem; }
  .lead { margin: 0 auto 1em; }
  .center { justify-content: center; margin: 14px 0 8px; }
  .detail { max-width: 600px; margin: 6px auto 16px; text-align: left; background: var(--panel-2); border: 1px solid var(--border); border-radius: 10px; padding: 4px 16px; }
  .detail ul { margin: 0 0 8px; padding-left: 20px; }
  .detail p { margin: 8px 0 4px; }
</style>
