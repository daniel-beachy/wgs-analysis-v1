<script lang="ts">
  // "What does this mean?" panel: collapsed by default, links to the matching primer section.
  import type { Snippet } from 'svelte';
  import { sectionFor } from '../lib/learn';
  import { go } from '../lib/router.svelte';
  let { title = 'What does this mean?', learn = '', open = false, children }: { title?: string; learn?: string; open?: boolean; children: Snippet } = $props();
  const target = $derived(learn ? sectionFor(learn) : undefined);
</script>

<details class="explain" {open}>
  <summary>💡 {title}</summary>
  <div class="body">
    {@render children()}
    {#if target}<button class="link" onclick={() => go('learn', target)}>Read more in the primer →</button>{/if}
  </div>
</details>

<style>
  .explain { border: 1px dashed var(--border); border-radius: 10px; padding: 8px 12px; margin: 10px 0; background: color-mix(in srgb, var(--accent) 5%, transparent); }
  summary { cursor: pointer; color: var(--accent); font-weight: 600; font-size: 0.9rem; }
  .body { padding-top: 8px; font-size: 0.92rem; }
  .body :global(p:last-child) { margin-bottom: 0.4em; }
  .link { background: none; border: 0; padding: 0; color: var(--accent); }
</style>
