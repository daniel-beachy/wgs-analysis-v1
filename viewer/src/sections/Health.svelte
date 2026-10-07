<script lang="ts">
  // Health: medically actionable genes first (with the reassurance of what was checked), then single-gene findings,
  // computer predictions, disagreements and risk factors. Each list is a graded ClaimList.
  import type { Manifest } from '../lib/data';
  import type { SectionDef } from '../lib/sections';
  import { query, useRelease, type Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { sectionFor } from '../lib/learn';
  import Card from '../components/Card.svelte';
  import ClaimList from '../components/ClaimList.svelte';
  import PgsList from '../components/PgsList.svelte';
  import Explain from '../components/Explain.svelte';
  import Term from '../components/Term.svelte';

  let { m, def, sex }: { m: Manifest; def: SectionDef; sex?: string } = $props();

  let genes = $state<Row[]>([]);
  let acmgClaims = $state<Row[]>([]);
  let pick = $state<Row | null>(null);
  let counts = $state<Record<string, number>>({});
  const hasAcmg = $derived(!!m.tables.acmg_genes);

  $effect(() => {
    (async () => {
      await useRelease(m);
      if (m.tables.acmg_genes) genes = (await query('SELECT * FROM acmg_genes ORDER BY category, gene')).rows;
      const c = await query(`SELECT * FROM claims WHERE acmg_sf IS NOT NULL AND category IN ('pathogenic', 'likely_pathogenic')`).catch(() => ({ rows: [] as Row[] }));
      acmgClaims = c.rows;
      const g = await query(`SELECT "group" AS g, count(*) AS n FROM claims WHERE section = 'health' GROUP BY 1`).catch(() => ({ rows: [] as Row[] }));
      counts = Object.fromEntries(g.rows.map((r) => [r.g, r.n]));
    })();
  });

  const reportable = $derived(acmgClaims.filter((c) => c.acmg_reportable));
  const byCat = $derived(Object.entries(genes.reduce((a: Record<string, Row[]>, g) => ((a[g.category] ??= []).push(g), a), {})));
  const meanCallable = $derived(genes.length ? genes.reduce((s, g) => s + (g.callable_fraction ?? 0), 0) / genes.length : 0);
  const CAT_ICON: Record<string, string> = { Cancer: '🎗️', Cardiovascular: '🫀', 'Inborn error of metabolism': '⚗️', Miscellaneous: '🧩' };
  const flagged = (g: string) => acmgClaims.find((c) => c.gene === g);
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);
</script>

<div class="grid">
  <div class="banner small"><span class="bi">{def.icon}</span>
    <div><strong>{def.title}</strong> — variants with documented links to disease, sorted from “a clinician would act on this” to “interesting but weak”. Every finding shows how strong the science is, how sure we are of your genotype, and who says so.</div>
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}
  </div>

  <Explain title="How to read this page (start here)" learn="From variant to meaning">
    <p>Almost everyone who looks at their genome finds a few “pathogenic” entries. That is normal: each of us carries a handful of disease variants, mostly as <Term t="Carrier" label="carriers" /> of recessive conditions or with <Term t="Penetrance" label="low penetrance" />. A finding here is a <em>prompt to learn and maybe ask a clinician</em>, not a diagnosis.</p>
    <ol>
      <li><strong>Actionable genes</strong> — the 84 genes where professional guidelines say a finding should be reported because something can be done about it.</li>
      <li><strong>Single-gene findings</strong> — variants labs have classified as (likely) pathogenic.</li>
      <li><strong>Computer predictions</strong> — rare variants no lab has classified, that prediction tools flag. Always Limited.</li>
      <li><strong>Labs disagree</strong> — at least one lab says pathogenic, others don’t.</li>
      <li><strong>Risk factors</strong> — common variants that nudge risk up or down a little.</li>
    </ol>
    <p>Findings about serious conditions are hidden until you click. The pills read: <strong>Evidence</strong> (how good the science is) · <strong>Your call</strong> (how sure we are of your genotype) · <strong>Overall</strong> (the weaker of the two).</p>
  </Explain>

  {#if hasAcmg}
    <Card title="1 · Medically actionable genes" subtitle={`ACMG secondary-findings list ${m.knowledge?.versions?.acmg_sf ?? ''}: genes where a pathogenic variant changes medical care`}>
      <div class="acmg-head">
        {#if reportable.length === 0}
          <div class="verdict good">✓</div>
          <div>
            <p class="big"><strong>No reportable findings in {genes.length} actionable genes.</strong></p>
            <p class="small muted">A clinical lab applying the <Term t="ACMG SF" /> rules to your variants would report nothing here. On average {(meanCallable * 100).toFixed(1)}% of each gene region was confidently readable in your data.
              {#if acmgClaims.length}{acmgClaims.length} variant{acmgClaims.length > 1 ? 's' : ''} in these genes {acmgClaims.length > 1 ? 'are' : 'is'} listed in ClinVar but {acmgClaims.length > 1 ? "don't" : "doesn't"} meet the reporting rule (e.g. one copy of a recessive variant).{/if}</p>
          </div>
        {:else}
          <div class="verdict warn">!</div>
          <div>
            <p class="big"><strong>{reportable.length} finding{reportable.length > 1 ? 's' : ''} a clinical lab would report.</strong></p>
            <p class="small muted">See “Single-gene findings” below (click to reveal). These are worth confirming with a clinical test and discussing with a doctor or genetic counsellor — many are manageable with screening.</p>
          </div>
        {/if}
      </div>
      <div class="cats">
        {#each byCat as [cat, gs]}
          <div class="cat">
            <div class="faint small">{CAT_ICON[cat] ?? ''} {cat} · {gs.length}</div>
            <div class="chips">
              {#each gs as g}
                {@const f = flagged(g.gene)}
                <button class="gene" class:flag={f?.acmg_reportable} class:note={f && !f.acmg_reportable} class:low={g.callable_fraction < 0.9}
                  class:sel={pick?.gene === g.gene} onclick={() => (pick = pick?.gene === g.gene ? null : g)}
                  title={`${g.conditions} · ${(g.callable_fraction * 100).toFixed(0)}% callable`}>{g.gene}</button>
              {/each}
            </div>
          </div>
        {/each}
      </div>
      {#if pick}
        <div class="pick small">
          <strong>{pick.gene}</strong> — {pick.conditions}. <span class="faint">Inheritance: {pick.inheritance}. Reported when: {pick.report}.</span><br />
          Readable in your data: <strong>{(pick.callable_fraction * 100).toFixed(1)}%</strong> of the {(pick.gene_bp / 1000).toFixed(0)} kb gene region.
          {#if flagged(pick.gene)}You have a ClinVar-listed variant here — see the lists below.{:else}No pathogenic or likely pathogenic variant found.{/if}
          <button class="link" onclick={() => go('variants', pick!.gene)}>Browse your variants in {pick.gene} →</button>
        </div>
      {/if}
      <p class="legend small faint"><span class="gene demo">GENE</span> nothing found · <span class="gene demo note">GENE</span> listed variant, not reportable · <span class="gene demo flag">GENE</span> reportable · <span class="gene demo low">GENE</span> under 90% readable. Click a gene for details.</p>
      <Explain title="What “no findings” does and doesn’t mean">
        <p>It means none of your small variants (single-letter changes and short insertions/deletions) in these genes is classified pathogenic in ClinVar in a way that meets the reporting rule. It does <em>not</em> rule out: larger deletions or duplications of whole exons, repeat expansions, variants in the few percent of each gene that couldn’t be read confidently, or variants not yet classified. A family history of one of these conditions matters more than this result — tell your doctor about it either way.</p>
      </Explain>
    </Card>
  {/if}

  <Card title="2 · Single-gene findings" subtitle="Variants labs have classified as pathogenic or likely pathogenic for a condition">
    <ClaimList {m} section="health" groups={['monogenic', 'compound']} empty="No pathogenic or likely pathogenic variants outside your carrier results." />
    <Explain title="Why “pathogenic” rarely means “you will get it”" learn="Monogenic vs polygenic">
      <p>“Pathogenic” describes the variant, not you. Whether it matters depends on <Term t="Inheritance" label="inheritance" /> (one copy is usually harmless for recessive conditions), <Term t="Penetrance" label="penetrance" /> (many dominant variants affect only a minority of carriers), and the rest of your genome and life. Each card shows which applies.</p>
    </Explain>
  </Card>

  <Card title="3 · Computer predictions" subtitle="Rare variants nobody has classified yet, flagged by prediction tools. Always Limited evidence.">
    <ClaimList {m} section="health" groups={['predicted']} empty="No rare, unclassified variants that the prediction tools agree are damaging in disease genes." />
    <Explain title="How the predictions work, and why they’re capped" learn="How this project grades">
      <p><Term t="AlphaMissense" /> (Google DeepMind) and <Term t="REVEL" /> score how likely a protein change is to be harmful; <Term t="LOEUF" /> says whether a gene tolerates broken copies. A variant only appears here if it is rare (under 0.1% in every population), sits in a known disease gene, and the tools agree (or REVEL alone reaches the ClinGen “strong” threshold). Even then it is a hypothesis: these tools are roughly 90% accurate on known variants, so across millions of variants they still raise false alarms. Recessive-gene predictions appear under Carrier.</p>
    </Explain>
  </Card>

  <Card title="4 · Labs disagree" subtitle="At least one lab says pathogenic, others say uncertain or benign">
    <ClaimList {m} section="health" groups={['uncertain']} empty="No conflicting classifications with a pathogenic vote." />
  </Card>

  <Card title="5 · Risk factors" subtitle="Common variants with small, documented effects on risk">
    <ClaimList {m} section="health" groups={['risk']} limit={12} empty="No risk-factor variants found." />
    <Explain title="Small effects add up" learn="Risk numbers">
      <p>Each risk factor shifts risk only a little, and most common diseases depend on thousands of variants plus lifestyle. The <Term t="Polygenic risk score" label="polygenic scores" /> below combine many variants into one percentile — a much better summary than any single entry here.</p>
    </Explain>
  </Card>

  <Card title="6 · Polygenic risk scores" subtitle="Your percentile for common diseases and health measures, against the 1000 Genomes reference people most similar to you. Click a score for details.">
    <PgsList {m} section="health" {sex} empty="Polygenic scores were not computed for this release." />
    <p class="small faint">A percentile is a tendency, not a forecast: lifestyle, family history and chance usually matter as much or more. Brain-related scores are on the <button class="link inline" onclick={() => go('brain')}>Brain &amp; mind</button> page.</p>
  </Card>

  {#if Object.keys(counts).length}
    <p class="small faint">This release: {Object.entries(counts).map(([g, n]) => `${n} ${g}`).join(' · ')}. Carrier results are on the <button class="link inline" onclick={() => go('carrier')}>Carrier</button> page.</p>
  {/if}
</div>

<style>
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
  .acmg-head { display: flex; gap: 16px; align-items: center; margin-bottom: 12px; }
  .acmg-head p { margin: 2px 0; }
  .big { font-size: 1.05rem; }
  .verdict { flex: none; width: 48px; height: 48px; border-radius: 50%; display: grid; place-items: center; font-size: 1.5rem; font-weight: 700; }
  .verdict.good { color: var(--good); background: color-mix(in srgb, var(--good) 16%, transparent); }
  .verdict.warn { color: var(--warn); background: color-mix(in srgb, var(--warn) 16%, transparent); }
  .cats { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 12px; }
  .chips { display: flex; flex-wrap: wrap; gap: 4px; margin-top: 4px; }
  .gene { font: 600 0.72rem/1 var(--mono, ui-monospace, monospace); padding: 4px 6px; border-radius: 6px; border: 1px solid var(--border); background: color-mix(in srgb, var(--good) 8%, var(--panel)); color: var(--text); cursor: pointer; }
  .gene.note { background: color-mix(in srgb, var(--warn) 18%, var(--panel)); border-color: color-mix(in srgb, var(--warn) 50%, transparent); }
  .gene.flag { background: color-mix(in srgb, var(--bad, #e5484d) 22%, var(--panel)); border-color: var(--bad, #e5484d); }
  .gene.low { border: 1.5px dashed var(--warn); }
  .gene.sel { outline: 2px solid var(--accent); }
  .gene.demo { cursor: default; padding: 2px 4px; }
  .pick { margin-top: 10px; padding: 8px 12px; border-radius: 10px; background: var(--panel-2); border: 1px solid var(--border); }
  .legend { margin-top: 10px; }
  .link { background: none; border: 0; padding: 0; color: var(--accent); }
  .link.inline { font-size: inherit; }
</style>
