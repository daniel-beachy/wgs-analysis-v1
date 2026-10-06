<script lang="ts">
  import { loadIndex, loadManifest, loadDoc, type Index, type Manifest } from './lib/data';
  import { route, go } from './lib/router.svelte';
  import { SECTION_DEFS } from './lib/sections';
  import Overview from './sections/Overview.svelte';
  import QC from './sections/QC.svelte';
  import Variants from './sections/Variants.svelte';
  import Sql from './sections/Sql.svelte';
  import Learn from './sections/Learn.svelte';
  import Upcoming from './sections/Upcoming.svelte';
  import Changes from './sections/Changes.svelte';

  const PORTFOLIO = 'https://daniel-beachy.github.io';
  let theme = $state(localStorage.getItem('theme') ?? 'dark');
  $effect(() => { document.documentElement.dataset.theme = theme; localStorage.setItem('theme', theme); });

  let index = $state<Index | null>(null);
  let releaseId = $state<string>('');
  let m = $state<Manifest | null>(null);
  let error = $state('');
  let sections = $state<any[]>([]);
  let sex = $state<string | undefined>();

  loadIndex().then((i) => { index = i; releaseId = i.latest ?? i.releases[0]?.id; })
    .catch((e) => (error = `Could not find any analysis results (${e}). Run the pipeline, then reopen the dashboard.`));

  $effect(() => {
    const entry = index?.releases.find((r) => r.id === releaseId);
    if (!entry) return;
    loadManifest(entry).then(async (mm) => {
      m = mm;
      sections = (await loadDoc(mm, 'sections')) ?? [];
      sex = (await loadDoc<any>(mm, 'qc'))?.inferred_sex?.call;
    }).catch((e) => (error = String(e)));
  });

  const def = $derived(SECTION_DEFS.find((s) => s.id === route.tab) ?? SECTION_DEFS[0]);
  $effect(() => { document.title = `${def.title} · Genome Dashboard`; });
</script>

<header>
  <div class="bar">
    <a class="brand" href="#/overview"><img src="./favicon.svg" alt="" width="26" height="26" /><span>Genome Dashboard</span></a>
    {#if m}<span class="faint small mono sample">{m.sample}</span>{/if}
    <span class="spacer"></span>
    {#if index && index.releases.length > 1}
      <label class="small muted">Release <select bind:value={releaseId}>{#each index.releases as r}<option value={r.id}>{r.id}</option>{/each}</select></label>
    {:else if m}<span class="small faint" title="Each pipeline run creates a dated release">release {m.id}</span>{/if}
    <button class="icon" title="Toggle light / dark" onclick={() => (theme = theme === 'dark' ? 'light' : 'dark')}>{theme === 'dark' ? '☀️' : '🌙'}</button>
    <a class="small" href={PORTFOLIO} target="_blank" rel="noreferrer" title="Daniel Beachy — portfolio">Portfolio ↗</a>
  </div>
  <nav>
    {#each SECTION_DEFS as s}
      <a href={`#/${s.id}`} class:active={def.id === s.id} class:soon={!s.ready} onclick={(e) => { e.preventDefault(); go(s.id); }}>{s.icon} {s.title}</a>
    {/each}
  </nav>
</header>

<main>
  {#if error}
    <div class="card"><h2>No data yet</h2><p>{error}</p><p class="small muted">The dashboard looks for <code>wgs-data/releases/index.json</code> next to (or above) the launcher.</p></div>
  {:else if !m}
    <p class="muted">Loading your genome…</p>
  {:else if def.id === 'overview'}<Overview {m} index={index!} />
  {:else if def.id === 'qc'}<QC {m} />
  {:else if def.id === 'changes'}<Changes {m} index={index!} />
  {:else if def.id === 'variants'}<Variants {m} {sex} />
  {:else if def.id === 'sql'}<Sql {m} />
  {:else if def.id === 'learn'}<Learn {theme} />
  {:else}<Upcoming {def} {m} status={sections.find((x) => x.id === (def.module ?? def.id))} />
  {/if}
</main>

<footer class="small faint">
  <p>For learning and curiosity — not a medical test or diagnosis. Discuss anything health-related with a clinician and confirm with a clinical-grade test.
    Runs entirely on this computer; your genome never leaves it. · <a href={PORTFOLIO} target="_blank" rel="noreferrer">daniel-beachy.github.io</a></p>
</footer>

<style>
  header { position: sticky; top: 0; z-index: 20; background: color-mix(in srgb, var(--bg) 88%, transparent); backdrop-filter: blur(10px); border-bottom: 1px solid var(--border); }
  .bar { display: flex; align-items: center; gap: 14px; padding: 10px 24px 4px; max-width: 1400px; margin: 0 auto; }
  .brand { display: flex; align-items: center; gap: 10px; font-weight: 700; font-size: 1.05rem; color: var(--text); }
  .brand:hover { text-decoration: none; }
  .icon { padding: 4px 8px; }
  nav { display: flex; gap: 2px; padding: 4px 18px 0; overflow-x: auto; max-width: 1400px; margin: 0 auto; scrollbar-width: none; }
  nav a { white-space: nowrap; padding: 8px 12px; color: var(--muted); border-bottom: 2px solid transparent; font-size: 0.9rem; }
  nav a:hover { color: var(--text); text-decoration: none; }
  nav a.active { color: var(--text); border-bottom-color: var(--accent); }
  nav a.soon { opacity: 0.65; }
  main { max-width: 1400px; margin: 0 auto; padding: 22px 24px 40px; }
  footer { max-width: 1400px; margin: 0 auto; padding: 0 24px 30px; }
  .sample { padding-top: 2px; }
  @media (max-width: 700px) { .sample { display: none; } }
</style>
