<script lang="ts">
  import type { Snippet } from 'svelte';
  let { title = '', subtitle = '', id = '', children, aside }: { title?: string; subtitle?: string; id?: string; children: Snippet; aside?: Snippet } = $props();
</script>

<section class="card" {id}>
  {#if title || aside}
    <header>
      <div>
        {#if title}<h3>{title}</h3>{/if}
        {#if subtitle}<p class="muted small">{subtitle}</p>{/if}
      </div>
      {#if aside}<div>{@render aside()}</div>{/if}
    </header>
  {/if}
  {@render children()}
</section>

<style>
  .card { background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 18px 20px; min-width: 0; }
  header { display: flex; justify-content: space-between; gap: 12px; align-items: flex-start; margin-bottom: 10px; }
  header p { margin: 0; }
  h3 { margin: 0 0 2px; }
</style>
