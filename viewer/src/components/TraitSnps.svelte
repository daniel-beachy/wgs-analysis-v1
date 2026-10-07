<script lang="ts">
  // Well-replicated single-variant traits (GWAS Catalog): your genotype, what the study found per copy of the
  // effect allele, and how common your genotype is. All wording comes from the pipeline (traits stage).
  import type { Manifest } from '../lib/data';
  import { query, useRelease, sqlString, type Row } from '../lib/db';
  import { OVERALL_PILL, CALL_PILL, pct } from '../lib/evidence';
  import Term from './Term.svelte';

  let { m, section = 'traits' }: { m: Manifest; section?: string } = $props();
  let rows = $state<Row[]>([]);
  let loaded = $state(false);
  let open = $state<string | null>(null);

  $effect(() => {
    const sec = section;
    (async () => {
      await useRelease(m);
      if (m.tables.trait_snps) rows = (await query(`SELECT * FROM trait_snps WHERE section = ${sqlString(sec)} ORDER BY body_system, label`)).rows;
      loaded = true;
    })();
  });
  const key = (r: Row) => `${r.rsid}-${r.label}`;
  const copies = (n: number | null) => (n == null ? 'unknown' : n === 0 ? 'no copies' : n === 1 ? 'one copy' : 'two copies');
</script>

{#if !loaded}<p class="muted small">Loading…</p>
{:else if !rows.length}<p class="muted small">No single-variant traits in this release.</p>
{:else}
  <div class="tablewrap">
    <table class="small">
      <thead><tr><th>Trait</th><th>You carry</th><th>What the study found</th><th><Term t="Percentile" label="Your genotype is shared by" /></th><th>Overall</th></tr></thead>
      <tbody>
        {#each rows as r (key(r))}
          <tr class:open={open === key(r)} onclick={() => (open = open === key(r) ? null : key(r))}>
            <td><strong>{r.label}</strong><div class="mono faint">{r.gene ?? ''} · {r.rsid}</div></td>
            <td><span class="mono">{r.genotype ?? '—'}</span>{#if r.effect_allele}<div class="faint">{copies(r.effect_copies)} of {r.effect_allele}</div>{/if}</td>
            <td>{r.effect_text ?? 'The study does not say which allele has the effect, so your genotype cannot be read in either direction.'}</td>
            <td>{#if r.genotype_share != null}<strong>{pct(r.genotype_share)}</strong> <span class="faint">of people</span>{#if r.genotype_share_eur != null}<div class="faint">{pct(r.genotype_share_eur)} of Europeans</div>{/if}{:else}<span class="faint">—</span>{/if}</td>
            <td><span class="pill {OVERALL_PILL[r.overall]}">{r.overall}</span></td>
          </tr>
          {#if open === key(r)}
            <tr class="detail"><td colspan="5">
              <div class="grid cols-3">
                <div><div class="faint"><Term t="Evidence grade" label="Evidence" />: <span class="pill {OVERALL_PILL[r.evidence_level]}">{r.evidence_level}</span></div><ul>{#each r.evidence_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
                <div><div class="faint"><Term t="Call confidence" label="Your call" />: <span class="pill {CALL_PILL[r.call_confidence]}">{r.call_confidence}</span></div><ul>{#each r.call_reasons ?? [] as x}<li>{x}</li>{/each}</ul></div>
                <div><div class="faint">The study</div><ul>
                  <li>“{r.trait_reported}”{r.sample ? ` — ${r.sample}` : ''}</li>
                  <li>{r.first_author ?? ''}{r.published ? ` (${String(r.published).slice(0, 4)})` : ''}{#if r.pmid} · <a href={`https://pubmed.ncbi.nlm.nih.gov/${r.pmid}/`} target="_blank" rel="noreferrer">PubMed</a>{/if}{#if r.study} · <a href={`https://www.ebi.ac.uk/gwas/studies/${r.study}`} target="_blank" rel="noreferrer">{r.study}</a>{/if}</li>
                  <li>{r.publications} publication{r.publications === 1 ? '' : 's'} report this variant for the trait{r.ci_text ? ` · 95% CI ${r.ci_text}` : ''}</li>
                  {#if r.effect_af != null}<li>{r.effect_allele} is on {pct(r.effect_af)} of chromosomes worldwide (1000 Genomes)</li>{/if}
                </ul></div>
              </div>
            </td></tr>
          {/if}
        {/each}
      </tbody>
    </table>
  </div>
{/if}

<style>
  tbody tr:not(.detail) { cursor: pointer; }
  tr.open td { background: color-mix(in srgb, var(--accent) 6%, transparent); }
  .detail ul { margin: 4px 0; padding-left: 18px; }
</style>
