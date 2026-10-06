<script lang="ts">
  // The genetics primer, rendered inside the dashboard (works offline). Mermaid diagrams load on demand.
  import { marked } from 'marked';
  import { PRIMER, SECTIONS, slug } from '../lib/learn';
  import { route, go } from '../lib/router.svelte';

  let { theme = 'dark' }: { theme?: string } = $props();
  let el: HTMLElement;

  const renderer = new marked.Renderer();
  renderer.heading = ({ tokens, depth }) => {
    const text = marked.Parser.parseInline(tokens);
    return `<h${depth} id="learn-${slug(text)}">${text}</h${depth}>`;
  };
  renderer.code = ({ text, lang }) =>
    lang === 'mermaid' ? `<pre class="mermaid">${text.replace(/</g, '&lt;')}</pre>` : `<pre><code>${text.replace(/&/g, '&amp;').replace(/</g, '&lt;')}</code></pre>`;
  // In-document links like (#8-reading-the-qc-report) become dashboard routes.
  renderer.link = ({ href, tokens }) => {
    const text = marked.Parser.parseInline(tokens);
    if (href.startsWith('#')) return `<a href="#/learn/${href.slice(1)}">${text}</a>`;
    if (/^https?:/.test(href)) return `<a href="${href}" target="_blank" rel="noreferrer">${text}</a>`;
    return text;
  };
  const html = marked.parse(PRIMER, { renderer, async: false }) as string;

  $effect(() => {
    const target = route.rest[0];
    if (target) requestAnimationFrame(() => document.getElementById(`learn-${target}`)?.scrollIntoView({ behavior: 'smooth', block: 'start' }));
  });

  $effect(() => {
    const t = theme;
    const blocks = el?.querySelectorAll<HTMLElement>('pre.mermaid');
    if (!blocks?.length) return;
    import('mermaid').then(({ default: mermaid }) => {
      mermaid.initialize({ startOnLoad: false, theme: t === 'dark' ? 'dark' : 'default', securityLevel: 'strict' });
      blocks.forEach((b) => { b.removeAttribute('data-processed'); if (b.dataset.src) b.innerHTML = b.dataset.src; else b.dataset.src = b.innerHTML; });
      mermaid.run({ nodes: Array.from(blocks) });
    });
  });
</script>

<div class="layout">
  <nav class="toc">
    <div class="faint small">Contents</div>
    {#each SECTIONS as s}
      <a class:sub={s.level === 3} href={`#/learn/${s.slug}`} onclick={(e) => { e.preventDefault(); go('learn', s.slug); }}>{s.title}</a>
    {/each}
  </nav>
  <article class="prose" bind:this={el}>{@html html}</article>
</div>

<style>
  .layout { display: grid; grid-template-columns: 240px minmax(0, 1fr); gap: 24px; }
  @media (max-width: 900px) { .layout { grid-template-columns: 1fr; } .toc { position: static !important; max-height: none !important; } }
  .toc { position: sticky; top: 76px; align-self: start; max-height: calc(100vh - 96px); overflow: auto; display: flex; flex-direction: column; gap: 2px; font-size: 0.86rem; }
  .toc a { color: var(--muted); padding: 3px 8px; border-radius: 6px; }
  .toc a:hover { background: var(--panel); color: var(--text); text-decoration: none; }
  .toc a.sub { padding-left: 20px; font-size: 0.8rem; }
  .prose { background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 26px 34px; max-width: 900px; }
  .prose :global(h1) { font-size: 1.8rem; }
  .prose :global(h2) { margin-top: 1.8em; padding-top: 0.6em; border-top: 1px solid var(--border); scroll-margin-top: 80px; }
  .prose :global(h3) { margin-top: 1.4em; scroll-margin-top: 80px; }
  .prose :global(img) { max-width: 100%; height: auto; border-radius: 10px; }
  .prose :global(table) { margin: 12px 0; font-size: 0.9rem; display: block; overflow-x: auto; }
  .prose :global(pre) { background: var(--code); padding: 12px 14px; border-radius: 8px; overflow-x: auto; font-size: 0.82rem; line-height: 1.4; }
  .prose :global(pre code) { background: none; padding: 0; }
  .prose :global(pre.mermaid) { background: transparent; text-align: center; }
  .prose :global(blockquote) { margin: 12px 0; padding: 8px 14px; border-left: 3px solid var(--accent); background: var(--panel-2); border-radius: 0 8px 8px 0; color: var(--muted); }
  .prose :global(li) { margin: 3px 0; }
</style>
