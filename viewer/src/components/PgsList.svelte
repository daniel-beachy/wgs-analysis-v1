<script lang="ts">
  // Polygenic scores for one dashboard section: your percentile among genetically similar reference people,
  // with the Model 3 grades (how well the score predicts in independent studies × how well it matched your data).
  import type { Manifest } from '../lib/data';
  import { query, useRelease, sqlString, type Row } from '../lib/db';
  import { OVERALL_PILL, CALL_PILL } from '../lib/evidence';
  import Term from './Term.svelte';

  let { m, section, sex, featured = true, empty = 'No polygenic scores for this section in this release.' }:
    { m: Manifest; section: string; sex?: string; featured?: boolean; empty?: string } = $props();

  let rows = $state<Row[]>([]);
  let loaded = $state(false);
  let error = $state('');
  let open = $state<string | null>(null);
  let browse = $state(false);

  $effect(() => {
    const sec = section;
    loaded = false;
    (async () => {
      try {
        await useRelease(m);
        if (!m.tables.pgs_scores) { rows = []; loaded = true; return; }
        rows = (await query(`SELECT * FROM pgs_scores WHERE section = ${sqlString(sec)}
          ORDER BY featured DESC, CASE overall WHEN 'Strong' THEN 0 WHEN 'Moderate' THEN 1 WHEN 'Limited' THEN 2 ELSE 3 END, label`)).rows;
      } catch (e) { error = String(e); }
      loaded = true;
    })();
  });

  const sexOk = (r: Row) => !r.sex || !sex || (sex === 'XY' ? r.sex === 'M' : sex === 'XX' ? r.sex === 'F' : true);
  const main = $derived(rows.filter((r) => r.featured && sexOk(r)));
  const other = $derived(rows.filter((r) => !r.featured || !sexOk(r)));
  const shown = $derived(featured ? (browse ? [...main, ...other] : main) : rows);
  const nth = (p: number) => { const n = Math.round(p); const s = n % 100 >= 11 && n % 100 <= 13 ? 'th' : ({ 1: 'st', 2: 'nd', 3: 'rd' } as any)[n % 10] ?? 'th'; return `${n}${s}`; };
  const band = (p: number) => (p >= 95 ? 'top 5%' : p >= 80 ? 'upper range' : p <= 5 ? 'bottom 5%' : p <= 20 ? 'lower range' : 'middle range');
  const POP: Record<string, string> = { EUR: 'European', AFR: 'African', EAS: 'East Asian', SAS: 'South Asian', AMR: 'admixed American' };
  const doi = (d: string | null) => (d ? `https://doi.org/${d}` : null);
  // PGS Catalog writes ancestry mixes as "European:72.4|African:19.3"; show "European 72%, African 19%".
  const ancestryMix = (a: string) =>
    a.split('|').map((x) => { const [n, v] = x.split(':'); return v ? `${n} ${Math.round(+v)}%` : n; }).join(', ');
</script>

{#if !loaded}
  <p class="muted small">Loading scores…</p>
{:else if error}
  <p class="small warn">{error}</p>
{:else if !rows.length}
  <p class="muted small">{empty}</p>
{:else}
  <div class="list">
    {#each shown as r (r.pgs_id)}
      {@const p = r.percentile}
      <div class="score" class:dim={!sexOk(r)}>
        <button class="head" onclick={() => (open = open === r.pgs_id ? null : r.pgs_id)} aria-expanded={open === r.pgs_id}>
          <div class="name"><strong>{r.label}</strong>
            <div class="faint small">{r.pgs_id} · {r.n_variants?.toLocaleString()} variants{!sexOk(r) ? ` · built for ${r.sex === 'F' ? 'women' : 'men'}` : ''}</div></div>
          <div class="bar" title={p != null ? `${nth(p)} percentile` : 'no percentile'}>
            {#if p != null}
              <div class="track"><span class="mid"></span><span class="mark" style:left="{p}%"></span></div>
              <div class="small"><strong>{nth(p)}</strong> <span class="faint">percentile · {band(p)}</span></div>
            {:else}<div class="small faint">not scored</div>{/if}
          </div>
          <div class="pills">
            <span class="pill {OVERALL_PILL[r.overall]}" title="Overall: the weaker of evidence and call">{r.overall}</span>
          </div>
        </button>
        {#if open === r.pgs_id}
          <div class="detail small">
            {#if r.trait_description}<p>{r.trait_description}{#if r.trait_url}{' '}<a href={r.trait_url} target="_blank" rel="noreferrer">ontology</a>{/if}</p>{/if}
            {#if p != null}<p>Your score is higher than about <strong>{Math.round(p)}%</strong> of the {POP[r.compared_with] ?? r.compared_with ?? 'reference'} reference group (1000 Genomes) — people whose DNA is most like yours. A higher percentile means more of the variants that push “{r.trait_reported ?? r.label}” up.</p>{/if}
            {#if r.eval_effect}<p><strong>How much the score explains:</strong> {r.eval_effect}{r.eval_n ? ` in ${r.eval_n.toLocaleString()} people` : ''}{r.eval_ancestry ? ` of ${r.eval_ancestry.replace(/^./, (c: string) => c.toUpperCase())} ancestry` : ''}{r.eval_author ? ` (${r.eval_author}${r.eval_published ? `, ${String(r.eval_published).slice(0, 4)}` : ''}${r.eval_independent ? ', independent study' : ', score authors’ own study'})` : ''}.</p>{/if}
            <div class="grid cols-3">
              <div><div class="faint"><Term t="Evidence grade" label="Evidence" />: <span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level}</span></div><ul>{#each r.evidence_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
              <div><div class="faint"><Term t="Call confidence" label="Your call" />: <span class="pill {CALL_PILL[r.call_confidence]}">{r.call_confidence}</span></div><ul>{#each r.call_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
              <div><div class="faint">About the score</div><ul>
                {#if r.name}<li>{r.name}{r.method ? ` — ${r.method}` : ''}</li>{/if}
                {#if r.ancestry_gwas}<li>Built from: {ancestryMix(r.ancestry_gwas)}</li>{/if}
                {#if r.matched != null}<li>{r.matched.toLocaleString()} of {r.total?.toLocaleString()} variants found in your data</li>{/if}
                {#if r.why}<li>Chosen because: {r.why}</li>{/if}
                <li>{#if r.first_author}{r.first_author}{r.published ? ` (${String(r.published).slice(0, 4)})` : ''}{r.journal ? `, ${r.journal}` : ''} · {/if}<a href={`https://www.pgscatalog.org/score/${r.pgs_id}/`} target="_blank" rel="noreferrer">PGS Catalog</a>{#if doi(r.doi)} · <a href={doi(r.doi)} target="_blank" rel="noreferrer">paper</a>{/if}</li>
              </ul></div>
            </div>
          </div>
        {/if}
      </div>
    {/each}
  </div>
  {#if featured && other.length}
    <button class="link small" onclick={() => (browse = !browse)}>{browse ? 'Show featured only' : `Browse ${other.length} more score${other.length > 1 ? 's' : ''} →`}</button>
  {/if}
{/if}

<style>
  .list { display: flex; flex-direction: column; gap: 6px; }
  .score { border: 1px solid var(--border); border-radius: 10px; background: var(--panel-2, var(--panel)); }
  .score.dim { opacity: 0.6; }
  .head { all: unset; box-sizing: border-box; cursor: pointer; display: grid; grid-template-columns: minmax(160px, 1.3fr) 2fr 92px; gap: 14px; align-items: center; width: 100%; padding: 9px 12px; }
  .head:focus-visible { outline: 2px solid var(--accent); border-radius: 10px; }
  .pills { justify-self: end; }
  @media (max-width: 700px) { .head { grid-template-columns: 1fr; } .pills { justify-self: start; } }
  .track { position: relative; height: 8px; border-radius: 4px; margin: 4px 0 3px;
    background: linear-gradient(90deg, color-mix(in srgb, var(--accent) 30%, transparent), var(--chip) 30%, var(--chip) 70%, color-mix(in srgb, var(--warn) 40%, transparent)); }
  .mid { position: absolute; left: 50%; top: -2px; width: 1px; height: 12px; background: var(--faint); }
  .mark { position: absolute; top: -4px; width: 4px; height: 16px; margin-left: -2px; border-radius: 2px; background: var(--text); box-shadow: 0 0 0 2px var(--panel); }
  .detail { padding: 0 12px 10px; border-top: 1px dashed var(--border); }
  .detail p { margin: 8px 0; }
  .detail ul { margin: 4px 0; padding-left: 18px; }
  .link { background: none; border: 0; padding: 6px 0; color: var(--accent); }
</style>
