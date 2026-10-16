<script lang="ts">
  // Ancestry (ADR-020): which present-day reference groups your DNA is most like, measured against the 1000 Genomes
  // panel with the method's own error rates, plus the maternal (mtDNA) and paternal (Y) lineages.
  import type { Manifest } from '../lib/data';
  import { loadDoc } from '../lib/data';
  import type { SectionDef } from '../lib/sections';
  import { query, useRelease, type Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { sectionFor } from '../lib/learn';
  import { OVERALL_PILL, CALL_PILL } from '../lib/evidence';
  import Card from '../components/Card.svelte';
  import Explain from '../components/Explain.svelte';
  import Term from '../components/Term.svelte';

  let { m, def }: { m: Manifest; def: SectionDef; sex?: string } = $props();
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);

  let comps = $state<Row[]>([]);
  let pops = $state<Row[]>([]);
  let haps = $state<Row[]>([]);
  let pts = $state<Row[]>([]);
  let doc = $state<any>(null);
  let loaded = $state(false);
  let error = $state('');
  let hover = $state<{ x: number; y: number; text: string } | null>(null);

  $effect(() => {
    (async () => {
      try {
        await useRelease(m);
        const t = m.tables;
        doc = await loadDoc(m, 'ancestry');
        if (t.ancestry_components) comps = (await query('SELECT * FROM ancestry_components ORDER BY share DESC')).rows;
        if (t.ancestry_populations) pops = (await query('SELECT * FROM ancestry_populations ORDER BY neighbours DESC, loo_recall DESC')).rows;
        if (t.haplogroups) haps = (await query("SELECT * FROM haplogroups ORDER BY line = 'paternal'")).rows;
        if (t.ancestry_pca) pts = (await query('SELECT * FROM ancestry_pca ORDER BY you')).rows;
      } catch (e) { error = String(e); }
      loaded = true;
    })();
  });

  // Categorical slots in fixed order (validated for colour-blind separation on both app surfaces).
  const SLOT = ['var(--s1)', 'var(--s2)', 'var(--s3)', 'var(--s4)', 'var(--s5)', 'var(--s6)', 'var(--s7)'];
  const SUPER = ['EUR', 'AFR', 'EAS', 'SAS', 'AMR'];
  const NAME: Record<string, string> = { AFR: 'African', EUR: 'European', EAS: 'East Asian', SAS: 'South Asian', AMR: 'Admixed American' };
  const superColour = (s: string) => SLOT[SUPER.indexOf(s)] ?? 'var(--faint)';
  const popCodes = $derived(pops.map((p) => p.code).sort());
  const popColour = (c: string) => SLOT[popCodes.indexOf(c)] ?? 'var(--faint)';
  const pct = (x: number) => (x >= 0.995 ? '100%' : x < 0.005 ? '0%' : `${(x * 100).toFixed(x < 0.1 ? 1 : 0)}%`);
  const detected = $derived(comps.filter((c) => c.detected));
  const main = $derived(doc?.main as string | undefined);
  const near = $derived(doc?.nearest);
  const zoom = $derived<number[]>(doc?.zoom_pcs ?? [3, 4]);
  const popName = (c: string) => pops.find((p) => p.code === c)?.description ?? c;
  const parse = (s: any) => { try { return typeof s === 'string' ? JSON.parse(s) : s ?? []; } catch { return []; } };

  // Scatter geometry
  const W = 420, H = 300, PAD = 28;
  function scale(rows: Row[], a: string, b: string) {
    const xs = rows.map((r) => r[a]), ys = rows.map((r) => r[b]);
    const [x0, x1] = [Math.min(...xs), Math.max(...xs)], [y0, y1] = [Math.min(...ys), Math.max(...ys)];
    return { x: (v: number) => PAD + ((v - x0) / (x1 - x0 || 1)) * (W - 2 * PAD), y: (v: number) => H - PAD - ((v - y0) / (y1 - y0 || 1)) * (H - 2 * PAD) };
  }
  const global = $derived(pts.length ? scale(pts, 'pc1', 'pc2') : null);
  const inGroup = $derived(pts.filter((p) => p.you || p.superpop === main));
  const local = $derived(inGroup.length ? scale(inGroup, `pc${zoom[0]}`, `pc${zoom[1]}`) : null);
  const you = $derived(pts.find((p) => p.you));
  const show = (e: MouseEvent, text: string) => { const r = (e.currentTarget as SVGElement).ownerSVGElement!.getBoundingClientRect(); hover = { x: e.clientX - r.left, y: e.clientY - r.top, text }; };
</script>

<div class="grid">
  <div class="banner small"><span class="bi">{def.icon}</span>
    <div><strong>{def.title}</strong> — {def.blurb}</div>
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}
  </div>

  <Explain title="Read this first: what genetic ancestry can and cannot say">
    <ul>
      <li><strong>It is similarity, not a family tree.</strong> Your DNA is compared with 2,504 present-day people from 26 populations sequenced by the 1000 Genomes Project. “European 100%” means your DNA sits with the European reference people, not that every ancestor lived in Europe.</li>
      <li><strong>Reference groups are broad.</strong> Five continental groups, and inside each a handful of sampled populations. Groups that were never sampled can’t show up, and neighbouring populations often can’t be told apart. The page says when that happens, using the method’s measured error rate.</li>
      <li><strong>Lineages are single lines.</strong> Your mother’s mother’s mother… line (mitochondrial DNA) and, for men, your father’s father’s father… line (Y chromosome). Each is one ancestor out of thousands in any generation, and the haplogroup names describe migrations thousands of years ago, not recent family.</li>
      <li>Commercial tests report finer regions using private reference panels with tens of thousands of people; this open method does not try to match them.</li>
    </ul>
  </Explain>

  {#if !loaded}<p class="muted small">Loading…</p>
  {:else if error}<p class="small warn">{error}</p>
  {:else if !comps.length && !haps.length}<p class="muted small">Ancestry was not computed for this release (needs the polygenic-score stage’s PCA or a VCF).</p>
  {:else}
    {#if comps.length}
      <Card title="1 · Continental make-up" subtitle="Your DNA written as the closest mix of the five 1000 Genomes continental groups. Shares below a group’s noise floor can’t be told apart from zero.">
        <div class="stack" role="img" aria-label={detected.map((c) => `${NAME[c.superpop]} ${pct(c.share)}`).join(', ')}>
          {#each detected as c (c.superpop)}<span style:width="{c.share * 100}%" style:background={superColour(c.superpop)} title="{NAME[c.superpop]} {pct(c.share)}"></span>{/each}
        </div>
        <div class="legend small">{#each detected as c (c.superpop)}<span><i style:background={superColour(c.superpop)}></i>{NAME[c.superpop]} <strong>{pct(c.share)}</strong></span>{/each}</div>
        <table class="t small">
          <thead><tr><th>Group</th><th>Your share</th><th><Term t="Noise floor" label="Noise floor" /></th><th>Grade</th><th>Why</th></tr></thead>
          <tbody>
            {#each comps as c (c.superpop)}
              <tr class:off={!c.detected}>
                <td><i class="dot" style:background={superColour(c.superpop)}></i>{NAME[c.superpop] ?? c.superpop}</td>
                <td data-label="Your share">{#if c.detected}<strong>{pct(c.share)}</strong>{:else}<span class="faint">not detected</span>{/if}</td>
                <td data-label="Noise floor">{pct(c.noise_floor)}</td>
                <td data-label="Grade"><span class="pill {OVERALL_PILL[c.overall]}">{c.overall}</span></td>
                <td class="faint" data-label="Why">{(c.evidence_reasons ?? []).join('. ')}.</td>
              </tr>
            {/each}
          </tbody>
        </table>
        {#if comps[0]?.call_reasons?.length}<p class="small">Your call: <span class="pill {CALL_PILL[comps[0].call_level]}">{comps[0].call_level}</span> {comps[0].call_reasons.join('. ')}.</p>{/if}
        <p class="faint small">How it works: the polygenic-score step places you and every reference person on 10 ancestry axes (<Term t="Principal components (PCA)" label="principal components" />, from {doc?.method?.pca_variants?.toLocaleString() ?? 'many'} positions spread across the genome). Your point is written as the closest mix of the five group centres (all shares ≥ 0, adding to 100%; the approach of rye, Conley et al. 2023). Each reference person is scored the same way with themselves left out: the noise floor is the share a group gets by chance among 95% of reference people who have none of it.</p>
      </Card>
    {/if}

    {#if pts.length && global && you}
      <Card title="2 · You among the reference people" subtitle="Each dot is one 1000 Genomes person; the ring is you. Left: the two main ancestry axes. Right: the people in your main group, on the two axes that best separate its populations.">
        <div class="plots">
          <figure>
            <svg viewBox="0 0 {W} {H}" role="img" aria-label="Ancestry PCA, axes 1 and 2">
              <text x={W / 2} y={H - 6} class="ax">ancestry axis 1</text>
              <text x="10" y={H / 2} class="ax" transform="rotate(-90 10 {H / 2})">axis 2</text>
              {#each pts.filter((p) => !p.you) as p (p.iid)}
                <circle cx={global.x(p.pc1)} cy={global.y(p.pc2)} r="2.6" fill={superColour(p.superpop)} opacity="0.7" role="presentation"
                  onmouseenter={(e) => show(e, `${p.iid} · ${NAME[p.superpop]} · ${p.population}`)} onmouseleave={() => (hover = null)} />
              {/each}
              <circle cx={global.x(you.pc1)} cy={global.y(you.pc2)} r="7" class="you" />
            </svg>
            <figcaption class="legend small">{#each SUPER as s}<span><i style:background={superColour(s)}></i>{NAME[s]}</span>{/each}<span><i class="ring"></i>You</span></figcaption>
          </figure>
          {#if local && main}
            <figure>
              <svg viewBox="0 0 {W} {H}" role="img" aria-label="Your main group, axes {zoom[0]} and {zoom[1]}">
                <text x={W / 2} y={H - 6} class="ax">ancestry axis {zoom[0]}</text>
                <text x="10" y={H / 2} class="ax" transform="rotate(-90 10 {H / 2})">axis {zoom[1]}</text>
                {#each inGroup.filter((p) => !p.you) as p (p.iid)}
                  <circle cx={local.x(p[`pc${zoom[0]}`])} cy={local.y(p[`pc${zoom[1]}`])} r="3" fill={popColour(p.population)} opacity="0.75" role="presentation"
                    onmouseenter={(e) => show(e, `${p.iid} · ${popName(p.population)}`)} onmouseleave={() => (hover = null)} />
                {/each}
                <circle cx={local.x(you[`pc${zoom[0]}`])} cy={local.y(you[`pc${zoom[1]}`])} r="7" class="you" />
              </svg>
              <figcaption class="legend small">{#each popCodes as c}<span><i style:background={popColour(c)}></i>{c}</span>{/each}<span><i class="ring"></i>You</span></figcaption>
            </figure>
          {/if}
          {#if hover}<div class="tip small" style:left="{hover.x + 12}px" style:top="{hover.y + 12}px">{hover.text}</div>{/if}
        </div>
      </Card>
    {/if}

    {#if pops.length && near}
      <Card title="3 · Nearest reference population" subtitle="Which {NAME[main ?? ''] ?? ''} populations your 10 nearest reference people come from, and how well the method separates them.">
        <p><strong>Closest:</strong> {near.group.map(popName).join(' or ')}{' '}<span class="pill {OVERALL_PILL[near.overall]}">{near.overall}</span></p>
        <ul class="small">{#each near.evidence_reasons ?? [] as r}<li>{r}.</li>{/each}</ul>
        <table class="t small">
          <thead><tr><th>Population</th><th>Your 10 nearest</th><th>Method accuracy</th><th>Often confused with</th></tr></thead>
          <tbody>
            {#each pops as p (p.code)}
              {@const cw = parse(p.confused_with)}
              <tr class:hl={p.in_nearest_group}>
                <td><i class="dot" style:background={popColour(p.code)}></i>{#if p.url}<a href={p.url} target="_blank" rel="noreferrer">{p.description ?? p.code}</a>{:else}{p.description ?? p.code}{/if} <span class="faint">({p.code}{p.samples ? `, ${p.samples} people` : ''})</span></td>
                <td data-label="Your 10 nearest"><span class="nb"><span style:width="{p.neighbours * 10}%"></span></span> {p.neighbours}</td>
                <td data-label="Method accuracy">{Math.round(p.loo_recall * 100)}%</td>
                <td class="faint" data-label="Often confused with">{Object.entries(cw).filter(([, v]) => (v as number) >= 0.05).map(([k, v]) => `${k} ${Math.round((v as number) * 100)}%`).join(', ') || '—'}</td>
              </tr>
            {/each}
          </tbody>
        </table>
        <p class="faint small">“Method accuracy” is the share of that population’s reference people placed back in it when each is left out and classified by its 10 nearest neighbours. Population names and sampling details come from the International Genome Sample Resource.</p>
      </Card>
    {/if}

    <Card title="4 · Maternal and paternal lineages" subtitle="Haplogroups: branches of the human family tree read from DNA passed down one line, almost unchanged.">
      <div class="lines">
        {#each haps as h (h.line)}
          {@const path = parse(h.path).filter((p: any) => p.state === 'D')}
          <div class="line">
            <div class="faint small">{h.line === 'maternal' ? 'Mother’s line · mitochondrial DNA' : 'Father’s line · Y chromosome'}</div>
            <div class="hg">{h.haplogroup}</div>
            <div class="pills"><span class="pill {OVERALL_PILL[h.overall]}">{h.overall}</span><span class="pill {CALL_PILL[h.call_level]}">Call {h.call_level}</span></div>
            {#if path.length}
              <div class="path small">{#each path.slice(-10) as p, i}{#if i}<span class="faint">{' › '}</span>{/if}<span title={p.marker ? `defined by ${p.marker}` : ''}>{p.haplogroup}</span>{/each}</div>
            {:else if parse(h.found).length}
              <div class="small"><span class="faint">Defining mutations found:</span> {parse(h.found).join(' ')}</div>
            {/if}
            <ul class="small">{#each h.call_reasons ?? [] as r}<li>{r}</li>{/each}{#each h.evidence_reasons ?? [] as r}<li>{r}</li>{/each}</ul>
            {#if parse(h.alternatives).length}<div class="small faint">Next-best fits: {parse(h.alternatives).map((a: any) => `${a.haplogroup} (${a.quality.toFixed(2)})`).join(', ')}</div>{/if}
            <div class="small">{h.tool} · {h.tree} · from your {h.input}{' · '}<a href={h.url} target="_blank" rel="noreferrer">This branch on YFull (age estimates, related samples)</a></div>
          </div>
        {/each}
        {#if doc?.paternal_note}<div class="line"><div class="faint small">Father’s line · Y chromosome</div><p class="small">{doc.paternal_note}.</p></div>{/if}
        {#if !haps.length}<p class="muted small">No haplogroups in this release (needs the haplo tools: <code>pixi install -e haplo</code>).</p>{/if}
      </div>
    </Card>
  {/if}
</div>

<style>
  /* Dark is the app default; light steps of the same hues under [data-theme='light']. */
  .grid { --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --s5: #d55181; --s6: #008300; --s7: #9085e9; }
  :global([data-theme='light']) .grid { --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --s4: #eda100; --s5: #e87ba4; --s6: #008300; --s7: #6250d6; }
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
  .stack { display: flex; gap: 2px; height: 18px; border-radius: 4px; overflow: hidden; margin: 6px 0; }
  .stack span { display: block; height: 100%; }
  .legend { display: flex; flex-wrap: wrap; gap: 4px 14px; margin: 4px 0 8px; }
  .legend i, .dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 5px; vertical-align: -1px; }
  .legend i.ring { background: none; border: 2px solid var(--text); width: 8px; height: 8px; }
  .t { width: 100%; border-collapse: collapse; margin: 6px 0; }
  .t th { text-align: left; color: var(--faint); font-weight: 600; padding: 4px 6px; border-bottom: 1px solid var(--border); }
  .t td { padding: 6px; border-bottom: 1px solid var(--border); vertical-align: top; }
  .t tr.off td { opacity: 0.65; }
  .t tr.hl td { background: color-mix(in srgb, var(--accent) 8%, transparent); }
  .nb { display: inline-block; width: 60px; height: 6px; border-radius: 3px; background: var(--chip); vertical-align: middle; margin-right: 6px; overflow: hidden; }
  .nb span { display: block; height: 100%; background: var(--text); }
  .plots { position: relative; display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 16px; }
  figure { margin: 0; }
  svg { width: 100%; height: auto; display: block; background: var(--panel-2, var(--panel)); border-radius: 10px; }
  .ax { fill: var(--faint); font-size: 11px; text-anchor: middle; }
  circle.you { fill: none; stroke: var(--text); stroke-width: 2.5; }
  .tip { position: absolute; pointer-events: none; background: var(--panel); border: 1px solid var(--border); border-radius: 6px; padding: 3px 8px; box-shadow: var(--shadow); white-space: nowrap; }
  .lines { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 12px; }
  .line { border: 1px solid var(--border); border-radius: 10px; padding: 10px 12px; }
  .hg { font-size: 1.6rem; font-weight: 700; margin: 2px 0 4px; }
  .pills { display: flex; gap: 6px; margin-bottom: 6px; }
  .path { margin: 4px 0; line-height: 1.6; }
  ul { margin: 4px 0; padding-left: 18px; }
  @media (max-width: 700px) { .t thead { display: none; } .t, .t tbody, .t tr, .t td { display: block; } .t td { border: 0; padding: 2px 6px; } .t td[data-label]::before { content: attr(data-label) ': '; color: var(--faint); } .t tr { border-bottom: 1px solid var(--border); padding: 4px 0; } }
</style>
