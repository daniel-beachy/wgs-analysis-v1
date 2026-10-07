<script lang="ts">
  // "Your risk vs typical risk" for one polygenic score (ADR-018). Every number names its source.
  let { a, compact = false }: { a: any; compact?: boolean } = $props();
  const r = $derived(typeof a === 'string' ? JSON.parse(a) : a);
  const pct = (x: number | null | undefined) => (x == null ? '' : x * 100 >= 10 ? Math.round(x * 100).toString() : (x * 100).toFixed(1).replace(/\.0$/, ''));
  const METHOD: Record<string, string> = {
    auroc: 'Converted from how well the score separated people who got the disease from those who did not (AUROC), using the standard liability-threshold (binormal) model.',
    or_per_sd: 'Converted from the published odds ratio per standard deviation of the score, calibrated so the average across everyone equals the typical risk.',
    hr_per_sd: 'Converted from the published hazard ratio per standard deviation of the score, calibrated so the average across everyone equals the typical risk.',
  };
</script>

{#if r}
  {#if compact}
    {#if r.status === 'ok'}
      <div class="small abs" title="Estimated absolute risk; open for sources">≈ <strong>{r.you_per_100}</strong> in 100 <span class="faint">vs {r.typical_per_100} typical</span></div>
    {/if}
  {:else if r.status === 'ok'}
    <div class="box">
      <div class="nums">
        <div><div class="big">{r.you_per_100}<span>&nbsp;in 100</span></div><div class="faint">people with your score{r.you_low != null && r.you_high != null && pct(r.you_low) !== pct(r.you_high) ? ` (range ${pct(r.you_low)}–${pct(r.you_high)})` : ''}</div></div>
        <div class="vs faint">vs</div>
        <div><div class="big typ">{r.typical_per_100}<span>&nbsp;in 100</span></div><div class="faint">typical person</div></div>
        <div class="ratio faint">{r.ratio >= 1 ? `${r.ratio.toFixed(1)}× typical` : `${(1 / r.ratio).toFixed(1)}× lower than typical`}</div>
      </div>
      <ul>
        <li><strong>Typical risk:</strong> {r.baseline.text}{' — '}{#if r.baseline.url}<a href={r.baseline.url} target="_blank" rel="noreferrer">{r.baseline.source}</a>{:else}{r.baseline.source}{/if}{#if r.baseline.quote}<blockquote>“{r.baseline.quote}”</blockquote>{/if}</li>
        <li><strong>Score effect:</strong> {r.effect.metric} {r.effect.estimate}{r.effect.ci?.[0] != null ? ` [${r.effect.ci[0]}–${r.effect.ci[1]}]` : ''}{r.effect.trait ? ` for “${r.effect.trait}”` : ''}{r.effect.cohort ? `, ${r.effect.cohort}` : ''}{r.effect.n ? ` (${r.effect.n.toLocaleString()} people, ${r.effect.cases?.toLocaleString()} cases)` : ''}{r.effect.ppm_id ? ` · ${r.effect.ppm_id}` : ''}{' — '}{#if r.effect.url}<a href={r.effect.url} target="_blank" rel="noreferrer">{r.effect.source}</a>{:else}{r.effect.source}{/if}{#if r.effect.pmid}{' · '}<a href={`https://pubmed.ncbi.nlm.nih.gov/${r.effect.pmid}/`} target="_blank" rel="noreferrer">PubMed {r.effect.pmid}</a>{/if}{#if r.effect.quote}<blockquote>“{r.effect.quote}”</blockquote>{/if}</li>
        <li class="faint">{METHOD[r.method] ?? r.method} Genes are only part of risk: lifestyle, family history and age matter as much or more, and a score never says you will or won’t get a disease.</li>
      </ul>
    </div>
  {:else}
    <p class="faint"><strong>Your risk in numbers:</strong> {r.status === 'not_applicable' ? 'not applicable' : 'not established'} — {r.reason}.</p>
  {/if}
{/if}

<style>
  .abs { margin-top: 2px; }
  .box { border: 1px solid var(--border); border-radius: 8px; padding: 8px 10px; margin: 8px 0; }
  .nums { display: flex; align-items: end; gap: 18px; flex-wrap: wrap; }
  .big { font-size: 1.5rem; font-weight: 700; line-height: 1.1; }
  .big span { font-size: 0.85rem; font-weight: 400; }
  .big.typ { color: var(--muted, var(--faint)); }
  .vs { padding-bottom: 18px; }
  .ratio { padding-bottom: 18px; margin-left: auto; }
  ul { margin: 6px 0 0; padding-left: 18px; }
  blockquote { margin: 3px 0 3px 4px; padding-left: 8px; border-left: 2px solid var(--border); font-style: italic; color: var(--faint); }
</style>
