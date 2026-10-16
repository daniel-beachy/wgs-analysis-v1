<script lang="ts">
  // Polygenic scores for one dashboard section, one row per trait (ADR-019): the headline is the median across the
  // studies behind the trait's scores, the band is the range of every score, and each score is listed with its paper,
  // how well it has been tested and how well it matched your data (Model 3 grades).
  import type { Manifest } from '../lib/data';
  import { query, useRelease, sqlString, type Row } from '../lib/db';
  import { OVERALL_PILL, CALL_PILL } from '../lib/evidence';
  import Term from './Term.svelte';
  import AbsRisk from './AbsRisk.svelte';

  let { m, section, sex, featured = true, empty = 'No polygenic scores for this section in this release.' }:
    { m: Manifest; section: string; sex?: string; featured?: boolean; empty?: string } = $props();

  let traits = $state<Row[]>([]);
  let scores = $state<Record<string, Row[]>>({});
  let loaded = $state(false);
  let error = $state('');
  let open = $state<string | null>(null);
  let openScore = $state<string | null>(null);
  let browse = $state(false);

  $effect(() => {
    const sec = section;
    loaded = false;
    (async () => {
      try {
        await useRelease(m);
        if (!m.tables.pgs_traits || !m.tables.pgs_scores) { traits = []; loaded = true; return; }
        traits = (await query(`SELECT * FROM pgs_traits WHERE section = ${sqlString(sec)}
          ORDER BY featured DESC, CASE overall WHEN 'Strong' THEN 0 WHEN 'Moderate' THEN 1 WHEN 'Limited' THEN 2 ELSE 3 END, label`)).rows;
        const rows = (await query(`SELECT * FROM pgs_scores WHERE trait_group IN
          (SELECT trait_id FROM pgs_traits WHERE section = ${sqlString(sec)})
          ORDER BY trait_group, independent_pubs DESC, independent_n DESC, evaluated_n DESC`)).rows;
        const by: Record<string, Row[]> = {};
        for (const r of rows) (by[r.trait_group] ??= []).push(r);
        scores = by;
      } catch (e) { error = String(e); }
      loaded = true;
    })();
  });

  const sexOk = (r: Row) => !r.sex || !sex || (sex === 'XY' ? r.sex === 'M' : sex === 'XX' ? r.sex === 'F' : true);
  const main = $derived(traits.filter((r) => r.featured && sexOk(r)));
  const other = $derived(traits.filter((r) => !r.featured || !sexOk(r)));
  const shown = $derived(featured ? (browse ? [...main, ...other] : main) : traits);
  const nth = (p: number) => { if (p < 0.5) return 'below 1st'; if (p > 99.5) return 'above 99th'; const n = Math.round(p); const s = n % 100 >= 11 && n % 100 <= 13 ? 'th' : ({ 1: 'st', 2: 'nd', 3: 'rd' } as any)[n % 10] ?? 'th'; return `${n}${s}`; };
  const band = (p: number) => (p >= 95 ? 'top 5%' : p >= 80 ? 'upper range' : p <= 5 ? 'bottom 5%' : p <= 20 ? 'lower range' : 'middle range');
  const POP: Record<string, string> = { EUR: 'European', AFR: 'African', EAS: 'East Asian', SAS: 'South Asian', AMR: 'admixed American' };
  const AGREE: Record<string, [string, string]> = {
    consistent: ['Scores agree', 'pass'], same_direction: ['Agree on direction', 'partial'],
    disagree: ['Scores disagree', 'warn'], single: ['One score', 'soon'], none: ['Not scored', 'fail'],
  };
  const doi = (d: string | null) => (d ? `https://doi.org/${d}` : null);
  const year = (d: any) => (d ? String(d).slice(0, 4) : '');
  const used = (r: Row) => r.passed && r.percentile != null;
  // PGS Catalog writes ancestry mixes as "European:72.4|African:19.3"; show "European 72%, African 19%".
  const ancestryMix = (a: string) =>
    a.split('|').map((x) => { const [n, v] = x.split(':'); return v ? `${n} ${Math.round(+v)}%` : n; }).join(', ');
  const studies = (t: Row): any[] => { try { return JSON.parse(t.studies ?? '[]'); } catch { return []; } };
</script>

{#snippet scoreDetail(r: Row)}
  <div class="sdetail">
    {#if used(r)}<p>This score puts you {#if r.percentile < 0.5}lower than <strong>over 99%</strong>{:else}higher than about <strong>{Math.round(r.percentile)}%</strong>{/if} of the {POP[r.compared_with] ?? r.compared_with ?? 'reference'} reference group (1000 Genomes) — people whose DNA is most like yours. A higher percentile means more of the variants that push “{r.trait_reported ?? r.label}” up.</p>{/if}
    <AbsRisk a={r.abs_risk} />
    {#if r.eval_effect}<p><strong>How much the score explains:</strong> {r.eval_effect}{r.eval_n ? ` in ${r.eval_n.toLocaleString()} people` : ''}{r.eval_ancestry ? ` of ${r.eval_ancestry.replace(/^./, (c: string) => c.toUpperCase())} ancestry` : ''}{r.eval_author ? ` (${r.eval_author}${r.eval_published ? `, ${year(r.eval_published)}` : ''}${r.eval_independent ? ', independent study' : ', score authors’ own study'})` : ''}{#if r.eval_pmid}{' · '}<a href={`https://pubmed.ncbi.nlm.nih.gov/${r.eval_pmid}/`} target="_blank" rel="noreferrer">PubMed</a>{/if}.</p>{/if}
    <div class="grid cols-3">
      <div><div class="faint"><Term t="Evidence grade" label="Evidence" />: <span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level}</span></div><ul>{#each r.evidence_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
      <div><div class="faint"><Term t="Call confidence" label="Your call" />: <span class="pill {CALL_PILL[r.call_confidence]}">{r.call_confidence}</span></div><ul>{#each r.call_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
      <div><div class="faint">About the score</div><ul>
        {#if r.name}<li>{r.name}{r.method ? ` — ${r.method}` : ''}</li>{/if}
        {#if r.ancestry_gwas}<li>Built from: {ancestryMix(r.ancestry_gwas)}</li>{/if}
        {#if r.matched != null}<li>{r.matched.toLocaleString()} of {r.total?.toLocaleString()} score positions readable in your data (your genotype was measured there; most are the common version, which is expected; coverage this high is what makes the score reliable)</li>{/if}
        {#if r.why}<li>{r.why}</li>{/if}
        <li>{r.rank_in_trait === 1 ? 'The selection rule’s pick for this trait' : `Ranked #${r.rank_in_trait} for this trait by the selection rule`} (most independent tests first)</li>
      </ul></div>
    </div>
  </div>
{/snippet}

{#if !loaded}
  <p class="muted small">Loading scores…</p>
{:else if error}
  <p class="small warn">{error}</p>
{:else if !traits.length}
  <p class="muted small">{empty}</p>
{:else}
  <div class="list">
    {#each shown as t (t.trait_id)}
      {@const p = t.percentile}
      {@const ss = scores[t.trait_id] ?? []}
      {@const ag = AGREE[t.agreement] ?? ['', '']}
      <div class="score" class:dim={!sexOk(t)}>
        <button class="head" onclick={() => (open = open === t.trait_id ? null : t.trait_id)} aria-expanded={open === t.trait_id}>
          <div class="name"><strong>{t.label}</strong>
            <div class="faint small">{t.n_scored} of {t.n_scores} score{t.n_scores === 1 ? '' : 's'} computed · {t.n_studies} stud{t.n_studies === 1 ? 'y' : 'ies'}{!sexOk(t) ? ` · built for ${t.sex === 'F' ? 'women' : 'men'}` : ''}</div></div>
          <div class="bar" title={p != null ? `median ${nth(p)} percentile; range ${nth(t.low)}–${nth(t.high)}` : 'no percentile'}>
            {#if p != null}
              <div class="track"><span class="mid"></span>
                {#if t.high > t.low}<span class="range" style:left="{t.low}%" style:width="{t.high - t.low}%"></span>{/if}
                {#each ss.filter(used) as s (s.pgs_id)}<span class="tick" style:left="{s.percentile}%"></span>{/each}
                <span class="mark" style:left="{p}%"></span></div>
              <div class="small"><strong>{nth(p)}</strong> <span class="faint">percentile{t.n_scored > 1 ? ` · range ${nth(t.low)}–${nth(t.high)}` : ''} · {band(p)}</span></div>
              <AbsRisk a={t.abs_risk} compact />
            {:else}<div class="small faint">not scored</div>{/if}
          </div>
          <div class="pills">
            {#if t.n_scored > 1}<span class="pill {ag[1]}" title={t.agreement_text}>{ag[0]}</span>{/if}
            <span class="pill {OVERALL_PILL[t.overall]}" title="Grade of the best-tested score ({t.best_pgs_id}): the weaker of its evidence and your call">{t.overall}</span>
          </div>
        </button>
        {#if open === t.trait_id}
          <div class="detail small">
            {#if t.trait_description}<p>{t.trait_description}{#if t.trait_url}{' '}<a href={t.trait_url} target="_blank" rel="noreferrer">ontology</a>{/if}</p>{/if}
            {#if p != null}
              <p><strong>Your result: {nth(p)} percentile</strong> among the {POP[t.compared_with] ?? t.compared_with ?? 'reference'} reference group (1000 Genomes).
                {t.n_studies > 1 ? `This is the middle (median) of ${t.n_studies} studies’ results` : 'This comes from one study'}{t.n_scored > t.n_studies ? `; several scores from one study count once` : ''}.
                {#if t.n_scored > 1}Every score is shown below; the shaded band is their range ({nth(t.low)} to {nth(t.high)}).{/if}</p>
              {#if t.n_scored > 1}<p><span class="pill {ag[1]}">{ag[0]}</span>{' '}{t.agreement_text}.</p>{/if}
              <div class="studies">{#each studies(t) as st}<span class="chip">{st.first_author ?? st.pgp_id}{st.published ? ` ${year(st.published)}` : ''}: <strong>{nth(st.percentile)}</strong>{st.pgs_ids.length > 1 ? ` (median of ${st.pgs_ids.length})` : ''}</span>{/each}</div>
            {/if}
            {#if t.risk_pgs_id}
              <AbsRisk a={t.abs_risk} />
              {#if JSON.parse(t.abs_risk ?? '{}').status === 'ok'}<p class="faint">Risk in numbers uses {t.risk_pgs_id}, the most-tested of these scores with a published effect that can be converted; it is that one score’s estimate, not the median’s.</p>{/if}
            {/if}
            <table class="scores">
              <thead><tr><th>Score</th><th>Study</th><th>You</th><th><Term t="Evidence grade" label="How well tested" /></th><th><Term t="Call confidence" label="Your call" /></th></tr></thead>
              <tbody>
                {#each ss as r (r.pgs_id)}
                  <tr class:off={!used(r)} class="srow" onclick={() => (openScore = openScore === r.pgs_id ? null : r.pgs_id)}>
                    <td><a href={`https://www.pgscatalog.org/score/${r.pgs_id}/`} target="_blank" rel="noreferrer" onclick={(e) => e.stopPropagation()}>{r.pgs_id}</a><div class="faint">{r.n_variants?.toLocaleString()} variants</div></td>
                    <td>{r.first_author ?? ''}{r.published ? ` ${year(r.published)}` : ''}{r.journal ? `, ${r.journal}` : ''}<div>{#if r.pmid}<a href={`https://pubmed.ncbi.nlm.nih.gov/${r.pmid}/`} target="_blank" rel="noreferrer" onclick={(e) => e.stopPropagation()}>PubMed</a>{/if}{#if r.pmid && doi(r.doi)}{' · '}{/if}{#if doi(r.doi)}<a href={doi(r.doi)} target="_blank" rel="noreferrer" onclick={(e) => e.stopPropagation()}>paper</a>{/if}</div></td>
                    <td>{#if used(r)}<strong>{nth(r.percentile)}</strong>{:else}<span class="faint">not scored</span>{/if}</td>
                    <td><span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level}</span><div class="faint">{r.independent_pubs ? `${r.independent_pubs} independent group${r.independent_pubs === 1 ? '' : 's'}, ${r.independent_n?.toLocaleString()} people` : `only its authors tested it${r.evaluated_n ? `, ${r.evaluated_n.toLocaleString()} people` : ''}`}</div></td>
                    <td><span class="pill {CALL_PILL[r.call_confidence]}">{r.call_confidence}</span>{#if r.match_rate != null}<div class="faint">{(r.match_rate * 100).toFixed(1)}% matched</div>{/if}</td>
                  </tr>
                  {#if openScore === r.pgs_id}<tr><td colspan="5">{@render scoreDetail(r)}</td></tr>{/if}
                {/each}
              </tbody>
            </table>
            <p class="faint">Click a score for its full detail. “How well tested” counts research groups other than the authors who checked the score in new people (Moderate needs at least 2 groups and 20,000 people; polygenic scores are never graded Strong).</p>
          </div>
        {/if}
      </div>
    {/each}
  </div>
  {#if featured && other.length}
    <button class="link small" onclick={() => (browse = !browse)}>{browse ? 'Show featured only' : `Browse ${other.length} more trait${other.length > 1 ? 's' : ''} →`}</button>
  {/if}
{/if}

<style>
  .list { display: flex; flex-direction: column; gap: 6px; }
  .score { border: 1px solid var(--border); border-radius: 10px; background: var(--panel-2, var(--panel)); }
  .score.dim { opacity: 0.6; }
  .head { all: unset; box-sizing: border-box; cursor: pointer; display: grid; grid-template-columns: minmax(160px, 1.3fr) 2fr 150px; gap: 14px; align-items: center; width: 100%; padding: 9px 12px; }
  .head:focus-visible { outline: 2px solid var(--accent); border-radius: 10px; }
  .pills { justify-self: end; display: flex; flex-direction: column; align-items: flex-end; gap: 4px; }
  @media (max-width: 700px) { .head { grid-template-columns: 1fr; } .pills { justify-self: start; flex-direction: row; } }
  .track { position: relative; height: 8px; border-radius: 4px; margin: 4px 0 3px;
    background: linear-gradient(90deg, color-mix(in srgb, var(--accent) 30%, transparent), var(--chip) 30%, var(--chip) 70%, color-mix(in srgb, var(--warn) 40%, transparent)); }
  .mid { position: absolute; left: 50%; top: -2px; width: 1px; height: 12px; background: var(--faint); }
  .range { position: absolute; top: 0; height: 8px; border-radius: 4px; background: color-mix(in srgb, var(--text) 22%, transparent); }
  .tick { position: absolute; top: -1px; width: 2px; height: 10px; margin-left: -1px; background: var(--faint); }
  .mark { position: absolute; top: -4px; width: 4px; height: 16px; margin-left: -2px; border-radius: 2px; background: var(--text); box-shadow: 0 0 0 2px var(--panel); }
  .detail { padding: 0 12px 10px; border-top: 1px dashed var(--border); }
  .detail p { margin: 8px 0; }
  .detail ul { margin: 4px 0; padding-left: 18px; }
  .studies { display: flex; flex-wrap: wrap; gap: 6px; margin: 4px 0 8px; }
  .chip { background: var(--chip); border-radius: 999px; padding: 1px 9px; }
  .scores { width: 100%; border-collapse: collapse; margin-top: 6px; }
  .scores th { text-align: left; font-weight: 600; color: var(--faint); padding: 4px 6px; border-bottom: 1px solid var(--border); }
  .scores td { padding: 6px; border-bottom: 1px solid var(--border); vertical-align: top; }
  .srow { cursor: pointer; }
  .srow:hover { background: color-mix(in srgb, var(--chip) 60%, transparent); }
  .srow.off { opacity: 0.6; }
  .sdetail { padding: 2px 4px 6px; }
  @media (max-width: 700px) { .scores, .scores tbody, .scores tr, .scores td { display: block; } .scores thead { display: none; } .scores td { border: 0; padding: 2px 6px; } .srow { border-bottom: 1px solid var(--border); padding: 4px 0; } }
  .link { background: none; border: 0; padding: 6px 0; color: var(--accent); }
</style>
