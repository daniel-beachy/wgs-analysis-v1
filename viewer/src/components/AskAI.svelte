<script lang="ts">
  // "Ask AI" (ADR-027): copies one result as a ready-made prompt. Nothing is sent; the person pastes it where they like.
  import type { Manifest } from '../lib/data';
  import { promptFor, copyText, type AskItem } from '../lib/askai';

  let { m, item }: { m: Manifest; item: () => AskItem | Promise<AskItem> } = $props();
  let state = $state<'' | 'copied' | 'failed'>('');

  async function ask(e: MouseEvent) {
    e.stopPropagation();
    state = (await copyText(promptFor(await item(), m))) ? 'copied' : 'failed';
    setTimeout(() => (state = ''), 2500);
  }
</script>

<button class="ask" onclick={ask}
  title="Copy this result, its grades, reasons and sources as a prompt for any AI assistant. Nothing is sent: pasting it into a cloud assistant shares this one result with that service.">
  {state === 'copied' ? '✓ Copied — paste into your assistant' : state === 'failed' ? 'Couldn’t copy' : '✨ Ask AI'}
</button>

<style>
  .ask { font-size: 0.85rem; padding: 3px 10px; white-space: nowrap; }
</style>
