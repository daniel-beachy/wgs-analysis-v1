<script lang="ts">
  // "What changed": each release stores a diff against the one before it (src/wgs/diff.py).
  // This shows why your results moved: a knowledge source updated, the evidence model changed, or your data did.
  import type { Index, Manifest, ReleaseEntry } from '../lib/data';
  import { loadManifest, loadDoc } from '../lib/data';
  import { dateTime, int } from '../lib/format';
  import { go } from '../lib/router.svelte';
  import { CLASS_LABEL, OVERALL_PILL } from '../lib/evidence';
  import Card from '../components/Card.svelte';
  import Explain from '../components/Explain.svelte';
  import Term from '../components/Term.svelte';

  let { m, index }: { m: Manifest; index: Index } = $props();

  const releases = $derived(index.releases.filter((r) => r.sample === m.sample));
  let pickedId = $state('');
  const picked = $derived(releases.find((r) => r.id === (pickedId || m.id)) ?? releases[0]);
  let revealed = $state(sessionStorage.getItem('reveal-sensitive') === '1');

  const load = $derived(picked ? loadManifest(picked).then(async (mm) => ({ mm, ch: await loadDoc<any>(mm, 'changes'), k: await loadDoc<any>(mm, 'knowledge') })) : Promise.resolve(null));

  const FIELD: Record<string, string> = {
    category: 'ClinVar classification', clinvar_stars: 'ClinVar review stars', evidence_level: 'Evidence grade',
    call_level: 'Call confidence', overall: 'Overall grade', gene_validity: 'Gene–disease validity (ClinGen)', section: 'Dashboard section',
  };
  const SECTION: Record<string, string> = { health: 'Health', carrier: 'Carrier', pgx: 'Medicines', traits: 'Traits' };
  const cls = (c: string | null) => (c ? CLASS_LABEL[c] ?? c : 'not in ClinVar');
  const fmt = (f: string, v: any) => (v == null ? '—' : f === 'clinvar_stars' ? `${v}★` : f === 'category' ? cls(v) : f === 'section' ? SECTION[v] ?? v : String(v));
  const where = (r: any) => r.rsid ?? `${r.chrom}:${r.pos}`;
  const masked = (r: any) => r.sensitive && !revealed;
  let showReworded = $state(false);

  // Word-level diff (LCS) so a wording update highlights only the words that changed.
  function wordDiff(a: string, b: string): { t: string; k: 'same' | 'del' | 'add' }[] {
    const x = (a ?? '').split(/(\s+)/), y = (b ?? '').split(/(\s+)/);
    const L = Array.from({ length: x.length + 1 }, () => new Uint16Array(y.length + 1));
    for (let i = x.length - 1; i >= 0; i--)
      for (let j = y.length - 1; j >= 0; j--) L[i][j] = x[i] === y[j] ? L[i + 1][j + 1] + 1 : Math.max(L[i + 1][j], L[i][j + 1]);
    const out: { t: string; k: 'same' | 'del' | 'add' }[] = [];
    const push = (t: string, k: 'same' | 'del' | 'add') => (out.at(-1)?.k === k ? (out.at(-1)!.t += t) : out.push({ t, k }));
    let i = 0, j = 0;
    while (i < x.length && j < y.length) {
      if (x[i] === y[j]) { push(x[i], 'same'); i++; j++; }
      else if (L[i + 1][j] >= L[i][j + 1]) push(x[i++], 'del');
      else push(y[j++], 'add');
    }
    while (i < x.length) push(x[i++], 'del');
    while (j < y.length) push(y[j++], 'add');
    return out;
  }
</script>

<div class="grid">
  <Card title="What changed" subtitle="Genomic knowledge keeps moving. Each release records exactly how your results differ from the release before it — and why.">
    <Explain learn="Keeping results current">
      <p>Your DNA does not change, but what science knows about it does. When you run <code>wgs knowledge refresh</code>, the tool downloads the newest ClinVar, gene and population data, then re-grades every finding. A new release is created only if something actually changed, and this page shows the difference: which knowledge sources moved to a new version, which findings appeared or disappeared, and which of your variants were <strong>reclassified</strong> (for example from “uncertain” to “benign”).</p>
    </Explain>
    <div class="timeline">
      {#each releases as r (r.id)}
        <button class="rel" class:on={picked?.id === r.id} onclick={() => (pickedId = r.id)}>
          <span class="dot"></span>
          <span class="mono small">{r.id}{r.id === index.latest ? ' · latest' : ''}</span>
          <span class="small muted">{(r as any).summary ?? (r.previous ? '' : 'First release.')}</span>
        </button>
      {/each}
    </div>
  </Card>

  {#await load}
    <p class="muted">Loading…</p>
  {:then data}
    {#if data}
      {@const ch = data.ch}
      {#if !ch}
        <Card><p class="muted">This release was made before change-tracking existed, so there is nothing to compare.</p></Card>
      {:else if ch.first}
        <Card title="First release"><p class="muted">Nothing to compare against yet. The next time knowledge is refreshed (or the pipeline changes), this page will show what moved.</p></Card>
      {:else}
        {@const cc = ch.claims.counts ?? {}}
        <div class="grid cols-4">
          <div class="card stat"><span class="kpi">{Object.keys(ch.knowledge).length}</span><span class="small muted">knowledge sources updated</span></div>
          <div class="card stat"><span class="kpi good">+{cc.added ?? 0}</span><span class="small muted">new findings</span></div>
          <div class="card stat"><span class="kpi bad">−{cc.removed ?? 0}</span><span class="small muted">findings no longer reported</span></div>
          <div class="card stat"><span class="kpi warn">{cc.regraded ?? cc.changed ?? 0}</span><span class="small muted">{`findings regraded${cc.regraded != null ? ` · ${int(cc.changed - cc.regraded)} wording updates` : ''} · ${int(ch.reclassified.count)} ClinVar reclassifications`}</span></div>
        </div>

        <Card title="Why things moved" subtitle={`Compared with release ${ch.previous}`}>
          <table class="small">
            <thead><tr><th>Source</th><th>Before</th><th>After</th></tr></thead>
            <tbody>
              {#each Object.entries(ch.knowledge as Record<string, [string, string]>) as [src, [a, b]]}
                <tr><td>{data.k?.sources?.[src]?.title ?? src}</td><td class="mono">{a ?? 'not used'}</td><td class="mono">{b ?? 'removed'}</td></tr>
              {:else}<tr><td colspan="3" class="muted">No knowledge source changed version.</td></tr>{/each}
              {#if ch.evidence_model}<tr><td>Evidence model (grading rules)</td><td class="mono">v{ch.evidence_model[0] ?? '—'}</td><td class="mono">v{ch.evidence_model[1]}</td></tr>{/if}
            </tbody>
          </table>
          {#if ch.claims.note}<p class="small faint">{ch.claims.note}</p>{/if}
        </Card>

        {@const regraded = ch.claims.changed.filter((r: any) => (r.kind ?? 'regraded') === 'regraded')}
        {@const reworded = ch.claims.changed.filter((r: any) => r.kind === 'reworded')}
        {#if ch.claims.added.length || ch.claims.removed.length || ch.claims.changed.length}
          <Card title="Your findings" subtitle="Graded claims that were added, removed or regraded. Wording-only updates (new condition names, updated submission counts) are folded away below.">
            {#if [...ch.claims.added, ...ch.claims.removed, ...ch.claims.changed].some((r: any) => r.sensitive)}
              <p class="small sens">🔒 Some changes involve sensitive findings and are masked. <button onclick={() => { revealed = !revealed; sessionStorage.setItem('reveal-sensitive', revealed ? '1' : '0'); }}>{revealed ? 'Hide' : 'Reveal'}</button></p>
            {/if}
            <table class="small">
              <thead><tr><th></th><th>Finding</th><th>Section</th><th>What changed</th><th>Now</th></tr></thead>
              <tbody>
                {#each regraded as r}
                  <tr><td><span class="pill partial">regraded</span></td>
                    {#if masked(r)}<td colspan="4" class="faint">🔒 sensitive finding</td>{:else}
                    <td><strong>{r.gene ?? r.subject}</strong> <span class="mono faint">{where(r)}</span><div class="faint">{r.category_label}</div></td>
                    <td>{SECTION[r.section] ?? r.section}</td>
                    <td>{#each Object.entries(r.fields as Record<string, [any, any]>) as [f, [a, b]]}{#if f === 'statement'}<div class="muted wd">{#each wordDiff(a, b) as seg}{#if seg.k === 'del'}<del>{seg.t}</del>{:else if seg.k === 'add'}<ins>{seg.t}</ins>{:else}{seg.t}{/if}{/each}</div>{:else}<div>{FIELD[f] ?? f}: <span class="old">{fmt(f, a)}</span> → <strong>{fmt(f, b)}</strong></div>{/if}{/each}</td>
                    <td><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></td>{/if}</tr>
                {/each}
                {#each ch.claims.added as r}
                  <tr><td><span class="pill pass">new</span></td>
                    {#if masked(r)}<td colspan="4" class="faint">🔒 sensitive finding</td>{:else}
                    <td><strong>{r.gene ?? r.subject}</strong> <span class="mono faint">{r.rsid ?? ''}</span><div class="faint">{r.category_label}</div></td>
                    <td>{SECTION[r.section] ?? r.section}</td><td class="muted">{r.statement}</td>
                    <td><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></td>{/if}</tr>
                {/each}
                {#each ch.claims.removed as r}
                  <tr><td><span class="pill fail">gone</span></td>
                    {#if masked(r)}<td colspan="4" class="faint">🔒 sensitive finding</td>{:else}
                    <td><strong>{r.gene ?? r.subject}</strong> <span class="mono faint">{r.rsid ?? ''}</span><div class="faint">{r.category_label}</div></td>
                    <td>{SECTION[r.section] ?? r.section}</td><td class="muted">No longer reported (e.g. reclassified as benign or withdrawn).</td><td></td>{/if}</tr>
                {/each}
                {#if reworded.length}
                  <tr><td colspan="5"><button onclick={() => (showReworded = !showReworded)}>{showReworded ? 'Hide' : 'Show'} {int(cc.changed - (cc.regraded ?? 0))} wording-only updates</button>
                    <span class="faint"> — same grade, but ClinVar's text changed (e.g. a condition renamed or a new lab submission).</span></td></tr>
                  {#if showReworded}
                    {#each reworded as r}
                      <tr><td><span class="pill">reworded</span></td>
                        {#if masked(r)}<td colspan="4" class="faint">🔒 sensitive finding</td>{:else}
                        <td><strong>{r.gene ?? r.subject}</strong> <span class="mono faint">{where(r)}</span><div class="faint">{r.category_label}</div></td>
                        <td>{SECTION[r.section] ?? r.section}</td>
                        <td class="muted">{#each wordDiff(r.fields.statement?.[0], r.fields.statement?.[1]) as seg}{#if seg.k === 'del'}<del>{seg.t}</del>{:else if seg.k === 'add'}<ins>{seg.t}</ins>{:else}{seg.t}{/if}{/each}</td>
                        <td><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></td>{/if}</tr>
                    {/each}
                  {/if}
                {/if}
              </tbody>
            </table>
            {#if (cc.added ?? 0) + (cc.removed ?? 0) + (cc.changed ?? 0) > ch.claims.added.length + ch.claims.removed.length + ch.claims.changed.length}
              <p class="small faint">Showing the first 500 of each kind.</p>{/if}
          </Card>
        {/if}

        {#if ch.reclassified.count}
          <Card title="ClinVar reclassifications among your variants" subtitle="Every variant you carry whose ClinVar class changed — including the many benign ones that never become findings.">
            <div class="grid cols-2">
              <table class="small">
                <thead><tr><th>Before</th><th>After</th><th class="num">Variants</th></tr></thead>
                <tbody>{#each ch.reclassified.by_change as b}<tr><td>{cls(b.old === 'not in ClinVar' ? null : b.old)}</td><td>{cls(b.new === 'not in ClinVar' ? null : b.new)}</td><td class="num">{int(b.n)}</td></tr>{/each}</tbody>
              </table>
              <div class="tablewrap">
                <table class="small">
                  <thead><tr><th>Variant</th><th>Gene</th><th>Before</th><th>After</th></tr></thead>
                  <tbody>
                    {#each ch.reclassified.items.slice(0, 200) as r}
                      <tr class="click" onclick={() => go('variants', where(r))}><td class="mono">{where(r)}</td><td>{r.gene ?? ''}</td>
                        <td>{cls(r.old_class)} {r.old_stars != null ? `${r.old_stars}★` : ''}</td><td>{cls(r.new_class)} {r.new_stars != null ? `${r.new_stars}★` : ''}</td></tr>
                    {/each}
                  </tbody>
                </table>
              </div>
            </div>
            {#if ch.reclassified.review_changed || ch.reclassified.detail_changed}<p class="small muted">Also: {int(ch.reclassified.review_changed ?? 0)} variants kept their class but changed review status (<Term t="Review stars" label="stars" />), and {int(ch.reclassified.detail_changed ?? 0)} changed only in wording (for example “Benign” → “Benign/Likely benign”). These are not counted as reclassifications.</p>{/if}
            {#if ch.reclassified.count > 200}<p class="small faint">Showing 200 of {int(ch.reclassified.count)}, pathogenic-related first. Query the full set in the SQL tab by comparing two releases’ annotations tables.</p>{/if}
          </Card>
        {/if}
      {/if}
      <p class="small faint">Release created {dateTime(data.mm.created)} · pipeline v{data.mm.pipeline_version} · evidence model v{data.mm.knowledge?.evidence_model ?? '—'}</p>
    {/if}
  {/await}
</div>

<style>
  .timeline { display: flex; flex-direction: column; gap: 2px; margin-top: 8px; border-left: 2px solid var(--border); padding-left: 0; }
  .rel { display: grid; grid-template-columns: 14px 170px 1fr; gap: 10px; align-items: center; text-align: left; background: none; border: 0; border-radius: 0 8px 8px 0; padding: 6px 10px 6px 0; margin-left: -7px; }
  .rel:hover, .rel.on { background: var(--panel-2); }
  .dot { width: 12px; height: 12px; border-radius: 50%; background: var(--border); border: 2px solid var(--panel); }
  .rel.on .dot { background: var(--accent); }
  .stat { display: flex; flex-direction: column; }
  .kpi.good { color: var(--good); } .kpi.bad { color: var(--bad); } .kpi.warn { color: var(--warn); }
  .old { text-decoration: line-through; color: var(--faint); }
  del { color: var(--bad); text-decoration: line-through; background: color-mix(in srgb, var(--bad) 12%, transparent); }
  ins { color: var(--good); text-decoration: none; background: color-mix(in srgb, var(--good) 14%, transparent); }
  .wd { margin-top: 4px; }
  .sens { background: color-mix(in srgb, var(--warn) 9%, transparent); padding: 6px 10px; border-radius: 8px; }
  .tablewrap { max-height: 420px; overflow: auto; }
  .click { cursor: pointer; } .click:hover { background: var(--panel-2); }
</style>
