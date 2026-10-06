<script lang="ts">
  import type { Index, Manifest } from '../lib/data';
  import { loadDoc } from '../lib/data';
  import { compact, dateTime, int } from '../lib/format';
  import { SECTION_DEFS } from '../lib/sections';
  import { go } from '../lib/router.svelte';
  import Card from '../components/Card.svelte';
  import Term from '../components/Term.svelte';

  let { m, index }: { m: Manifest; index: Index } = $props();
  const load = $derived(Promise.all([loadDoc(m, 'qc'), loadDoc(m, 'sections'), loadDoc(m, 'coverage'), loadDoc(m, 'knowledge'), loadDoc(m, 'annotate')]));
  const tiles = SECTION_DEFS.filter((s) => s.id !== 'overview');
</script>

{#await load then [qc, sections, cov, know, ann]}
  {@const checks = qc?.assessment ?? []}
  {@const passed = checks.filter((c: any) => c.status === 'pass').length}
  <div class="grid">
    <section class="hero">
      <div>
        <p class="eyebrow">Your genome · release {m.id}</p>
        <h1>About 3 billion letters, read {cov ? `${Math.round(cov.autosomal_mean_depth)} times over` : 'many times over'}.</h1>
        <p class="lead">This dashboard turns your whole-genome sequencing data into something you can explore and keep learning from. Every claim it makes will show the evidence behind it and how confident it is — and when the science changes, a refresh shows exactly what changed for you.</p>
        <div class="row">
          <button class="primary" onclick={() => go('learn')}>📖 New to genetics? Start here</button>
          <button onclick={() => go('qc')}>✅ Check the data quality</button>
          <button onclick={() => go('variants')}>🔎 Explore variants</button>
        </div>
      </div>
      <div class="facts">
        {#if qc}
          <div><span class="kpi">{compact(qc.metrics.snv_count + qc.metrics.indel_count)}</span><span class="muted small">confident <Term t="Variant" label="variants" /></span></div>
          <div><span class="kpi">{qc.metrics.autosomal_mean_depth}×</span><span class="muted small">average <Term t="Coverage" label="depth" /></span></div>
          <div><span class="kpi">{qc.inferred_sex?.call}</span><span class="muted small">chromosomal sex</span></div>
          <div><span class="kpi">{passed}/{checks.length}</span><span class="muted small">quality checks pass</span></div>
        {/if}
      </div>
    </section>

    <div class="grid cols-3">
      {#each tiles as s}
        {@const st = s.module ? sections?.find((x: any) => x.id === s.module) : undefined}
        <button class="tile" onclick={() => go(s.id)}>
          <div class="row"><span class="ti">{s.icon}</span><strong>{s.title}</strong><span class="spacer"></span>
            {#if s.ready}<span class="pill pass">ready</span>{:else}<span class="pill soon">step {s.step}</span>{/if}</div>
          <p class="small muted">{s.blurb}</p>
          {#if st && !s.ready}<p class="small faint">Your data: <span class="pill {st.status}">{st.status}</span></p>{/if}
        </button>
      {/each}
    </div>

    <div class="grid cols-2">
      <Card title="How confident is each claim?" subtitle="Every finding carries two grades — how strong the science is, and how sure we are of your genotype. The overall grade is the weaker of the two.">
        <table class="small">
          <tbody>
            <tr><td><span class="pill pass">Strong</span></td><td>Expert-reviewed, replicated, and your genotype is confidently called.</td></tr>
            <tr><td><span class="pill partial">Moderate</span></td><td>Solid evidence with caveats, e.g. a single study or a modest effect.</td></tr>
            <tr><td><span class="pill">Limited</span></td><td>Preliminary, conflicting or tiny effect — shown for curiosity and learning.</td></tr>
            <tr><td><span class="pill fail">Not callable</span></td><td>Your data can’t answer this position (low coverage), so no claim is made.</td></tr>
          </tbody>
        </table>
        <button class="link small" onclick={() => go('learn', 'how-this-project-grades-its-confidence')}>How the grading works →</button>
        <p class="small faint">This is a personal learning tool, not a medical test. Anything health-related should be confirmed by a clinician with a clinical-grade test.</p>
      </Card>
      <Card title="Releases" subtitle="Each analysis run is saved as a dated, read-only release, so results can be compared over time.">
        <table class="small">
          <thead><tr><th>Release</th><th>Created</th><th>What changed</th></tr></thead>
          <tbody>
            {#each index.releases.filter((r) => r.sample === m.sample) as r}
              <tr><td class="mono">{r.id}{r.id === m.id ? ' ●' : ''}</td><td>{dateTime(r.created)}</td><td class="small muted">{r.summary ?? ''}</td></tr>
            {/each}
          </tbody>
        </table>
        <p class="small faint">Pipeline v{m.pipeline_version}. <button class="link" onclick={() => go('changes')}>See what changed →</button></p>
      </Card>
    </div>

    {#if know?.sources}
      <Card title="Where the knowledge comes from" subtitle="Public databases, downloaded and joined to your variants on this computer — your genome is never sent anywhere. Refresh with `wgs knowledge refresh`.">
        {#if ann}
          <div class="row small muted annstats">
            <span><strong>{int(ann.alleles)}</strong> variant alleles annotated</span>
            <span><strong>{int(ann.high_impact)}</strong> high-impact</span>
            <span><strong>{int(ann.in_clinvar)}</strong> in ClinVar</span>
            <span><strong>{int(ann.claims)}</strong> graded findings</span>
          </div>
        {/if}
        <table class="small">
          <thead><tr><th>Source</th><th>Version used</th><th>What it adds</th><th>Updates</th><th>Licence</th></tr></thead>
          <tbody>
            {#each Object.entries(know.sources) as [id, s]}
              {@const src = s as any}
              <tr><td><a href={src.homepage} target="_blank" rel="noreferrer">{src.title ?? id}</a></td><td class="mono">{src.version}</td>
                <td class="muted">{src.description}</td><td class="muted">{src.cadence}</td><td class="faint">{src.licence}</td></tr>
            {/each}
            <tr><td>Evidence model</td><td class="mono">v{know.evidence_model}</td><td class="muted">The grading rules themselves (ADR-013). Changing them is tracked like a data update.</td><td></td><td></td></tr>
          </tbody>
        </table>
      </Card>
    {/if}
  </div>
{/await}

<style>
  .hero { display: grid; grid-template-columns: minmax(0, 1.6fr) minmax(220px, 1fr); gap: 28px; padding: 30px; border-radius: 18px; border: 1px solid var(--border);
    background: radial-gradient(circle at 85% 20%, color-mix(in srgb, var(--accent-2) 18%, transparent), transparent 45%), radial-gradient(circle at 10% 90%, color-mix(in srgb, var(--accent) 16%, transparent), transparent 50%), var(--panel); }
  @media (max-width: 860px) { .hero { grid-template-columns: 1fr; } }
  .eyebrow { text-transform: uppercase; letter-spacing: 0.08em; font-size: 0.75rem; color: var(--accent); font-weight: 700; margin-bottom: 6px; }
  h1 { font-size: 2rem; max-width: 22ch; }
  .facts { display: grid; grid-template-columns: 1fr 1fr; gap: 14px; align-content: center; }
  .facts > div { display: flex; flex-direction: column; background: color-mix(in srgb, var(--panel-2) 80%, transparent); border: 1px solid var(--border); border-radius: 12px; padding: 12px 14px; }
  .tile { display: flex; flex-direction: column; justify-content: flex-start; text-align: left; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 14px 16px; color: var(--text); }
  .tile:hover { border-color: var(--accent); transform: translateY(-1px); }
  .tile p { margin: 6px 0 0; }
  .ti { font-size: 1.2rem; }
  .link { background: none; border: 0; padding: 0; color: var(--accent); }
  .annstats { gap: 18px; margin-bottom: 8px; }
</style>
