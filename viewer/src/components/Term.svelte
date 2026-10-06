<script lang="ts">
  // A glossary term: dotted underline; hover or focus shows the plain-language definition.
  import { lookup, sectionFor } from '../lib/learn';
  import { go } from '../lib/router.svelte';
  let { t, label = '', learn = '' }: { t: string; label?: string; learn?: string } = $props();
  const g = $derived(lookup(t));
  const target = $derived(sectionFor(learn || t));
</script>

<span class="term" tabindex="0" role="button" aria-label={g ? `${g.term}: ${g.meaning}` : t}>
  {label || t}
  {#if g || target}
    <span class="tip" role="tooltip">
      {#if g}<strong>{g.term}</strong><span>{g.meaning}</span>{/if}
      {#if target}<button class="link" onclick={() => go('learn', target)}>Learn more →</button>{/if}
    </span>
  {/if}
</span>

<style>
  .term { position: relative; border-bottom: 1px dotted var(--faint); cursor: help; }
  .tip {
    display: none; position: absolute; z-index: 50; left: 0; top: calc(100% + 6px); width: max-content; max-width: 300px;
    background: var(--panel-2); border: 1px solid var(--border); border-radius: 10px; padding: 10px 12px;
    box-shadow: var(--shadow); font-size: 0.85rem; font-weight: 400; color: var(--text); line-height: 1.45;
    text-transform: none; letter-spacing: 0; white-space: normal;
  }
  .tip::before { content: ''; position: absolute; left: 0; right: 0; top: -8px; height: 8px; }
  .tip strong { display: block; margin-bottom: 2px; }
  .term:hover .tip, .term:focus-within .tip, .term:focus .tip { display: block; }
  .link { background: none; border: 0; padding: 4px 0 0; color: var(--accent); display: block; }
</style>
