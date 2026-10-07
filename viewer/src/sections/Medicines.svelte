<script lang="ts">
  // Medicines (pharmacogenomics): which versions of the drug-processing genes you carry, what that does to
  // each enzyme's "speed", and which published prescribing guidelines change because of it.
  // Data: pgx_genes / pgx_drugs (PharmCAT + Cyrius + T1K, built by stages/pgx.py), claims (section pgx).
  import type { Manifest } from '../lib/data';
  import type { SectionDef } from '../lib/sections';
  import { DATA } from '../lib/data';
  import { query, useRelease, type Row } from '../lib/db';
  import { go } from '../lib/router.svelte';
  import { sectionFor } from '../lib/learn';
  import { OVERALL_PILL, CALL_PILL } from '../lib/evidence';
  import Card from '../components/Card.svelte';
  import ClaimList from '../components/ClaimList.svelte';
  import Explain from '../components/Explain.svelte';
  import Term from '../components/Term.svelte';

  let { m, def }: { m: Manifest; def: SectionDef } = $props();

  let genes = $state<Row[]>([]);
  let drugs = $state<Row[]>([]);
  let about = $state<Record<string, Row>>({});
  let loaded = $state(false);
  let error = $state('');
  let openGene = $state<string | null>(null);
  let openDrug = $state<string | null>(null);
  let search = $state('');
  let showStandard = $state(false);
  let showFda = $state(false);
  let showNone = $state(false);

  $effect(() => {
    loaded = false;
    (async () => {
      try {
        await useRelease(m);
        genes = (await query('SELECT * FROM pgx_genes ORDER BY gene')).rows;
        drugs = (await query('SELECT * FROM pgx_drugs ORDER BY drug, matched DESC, source')).rows;
        about = {};
        if (m.tables.gene_about) {
          for (const r of (await query('SELECT * FROM gene_about WHERE gene IN (SELECT gene FROM pgx_genes)')).rows) about[r.gene] = r;
        }
      } catch (e) { error = String(e); }
      loaded = true;
    })();
  });

  // Enzyme "speed" scale for the gauge diagram.
  const SPEED = ['poor', 'intermediate', 'normal', 'rapid', 'ultrarapid'];
  const SPEED_LABEL = ['Poor', 'Intermediate', 'Normal', 'Rapid', 'Ultrarapid'];
  const speedIndex = (g: Row) => {
    const p = (g.phenotype ?? '').toLowerCase();
    if (p.includes('ultrarapid')) return 4;
    if (p.includes('rapid')) return 3;
    if (p.includes('intermediate') || p.includes('decreased')) return 1;
    if (p.includes('poor')) return 0;
    if (p.includes('normal')) return 2;
    return -1;
  };
  const gauged = $derived(genes.filter((g) => speedIndex(g) >= 0 && !g.gene.startsWith('HLA') && g.gene !== 'G6PD' && g.gene !== 'MT-RNR1'));
  const notable = $derived(genes.filter((g) => speedIndex(g) >= 0 && speedIndex(g) !== 2));
  const other = $derived(genes.filter((g) => !gauged.includes(g)));
  const untyped = $derived(genes.filter((g) => g.call_level === 'Not callable'));

  // Drug list: group the per-source rows by drug. The bucket is the strongest CPIC/DPWG tier, decided in the
  // pipeline from ClinPGx's curated flags (stages/pgx.py `tier`). FDA label rows are shown but never set the bucket.
  const GUIDE = new Set(['CPIC', 'DPWG']);
  const RANK = ['none', 'standard', 'note', 'change'];
  const byDrug = $derived.by(() => {
    const out = new Map<string, Row[]>();
    for (const d of drugs) (out.get(d.drug) ?? out.set(d.drug, []).get(d.drug)!).push(d);
    return [...out.entries()].map(([drug, rows]) => {
      const guideAll = rows.filter((r) => GUIDE.has(r.source));
      const guide = guideAll.filter((r) => r.matched);
      // No CPIC/DPWG guideline: an FDA label that matched you, else listed as "not for your result" (never hidden).
      const tier = guideAll.length ? guideAll.reduce((t, r) => (RANK.indexOf(r.tier) > RANK.indexOf(t) ? r.tier : t), 'none')
        : rows.some((r) => r.matched) ? 'fda' : 'none';
      const best = guide.find((r) => r.overall === 'Strong') ?? guide.find((r) => r.overall === 'Moderate') ?? guide[0];
      // CPIC and DPWG reach different conclusions for your genotype (e.g. venlafaxine: "no action" vs "avoid").
      const bySrc = new Map<string, string>();
      for (const r of guide) bySrc.set(r.source, RANK.indexOf(r.tier) > RANK.indexOf(bySrc.get(r.source) ?? 'none') ? r.tier : bySrc.get(r.source)!);
      const differ = bySrc.size > 1 && new Set(bySrc.values()).size > 1;
      return { drug, tier, differ, rows: rows.filter((r) => r.matched), guide,
        notes: rows.filter((r) => !r.matched && r.tier === 'note'), unmatched: rows.some((r) => r.matched || r.tier === 'note') ? [] : rows,
        genes: [...new Set(rows.flatMap((r) => r.genes ?? []))].sort(), overall: best?.overall ?? null };
    });
  });
  const q = $derived(search.trim().toLowerCase());
  const match = (d: { drug: string; genes: string[] }) => !q || d.drug.toLowerCase().includes(q) || d.genes.some((g) => g.toLowerCase().includes(q));
  const inTier = (t: string) => byDrug.filter((d) => d.tier === t && match(d));
  const actionDrugs = $derived(inTier('change'));
  const noteDrugs = $derived(inTier('note'));
  const standardDrugs = $derived(inTier('standard'));
  const noneDrugs = $derived(inTier('none'));
  const fdaOnly = $derived(inTier('fda'));
  // Gene description: MedlinePlus (written for the public, expert-reviewed) first, else NCBI Gene's summary.
  const aboutGene = (g: string) => {
    const a = about[g];
    if (a?.medlineplus_text) return { name: a.name, text: a.medlineplus_text, src: 'MedlinePlus Genetics', url: a.medlineplus_url, when: a.medlineplus_reviewed ? `reviewed ${a.medlineplus_reviewed}` : '', machine: false };
    if (a?.ncbi_summary) {
      const machine = !['RefSeq', 'OMIM'].includes(a.ncbi_summary_source);
      return { name: a.name, text: a.ncbi_summary, src: `NCBI Gene (${a.ncbi_summary_source})`, url: a.ncbi_url, machine,
        when: machine ? 'machine-generated from gene annotations; can describe unrelated functions' : `written by ${a.ncbi_summary_source} curators` };
    }
    return a?.name ? { name: a.name, text: '', src: 'NCBI Gene', url: a.ncbi_url, when: '', machine: false } : null;
  };
  const firstPara = (t: string) => t.split(/\n/)[0];

  // Guideline text arrives as HTML fragments (lists, <br>, entities): turn it into plain text with bullets.
  function plain(html: string | null | undefined): string {
    if (!html) return '';
    const s = html.replace(/<li>/gi, '\n• ').replace(/<br\s*\/?>/gi, '\n').replace(/<\/(p|ul|ol)>/gi, '\n');
    const t = new DOMParser().parseFromString(s, 'text/html').body.textContent ?? '';
    return t.replace(/\n{3,}/g, '\n\n').trim();
  }
  const json = (s: string | null | undefined, d: any = null) => { try { return s ? JSON.parse(s) : d; } catch { return d; } };
  const reportUrl = $derived(m.documents.pgx_report ? DATA + m.documents.pgx_report.path : null);
  const learnSlug = $derived(def.learn ? sectionFor(def.learn) : undefined);
  const speedColor = (i: number) => ['var(--bad)', 'var(--warn)', 'var(--good)', 'var(--accent)', 'var(--accent-2)'][i];
  // CPIC population codes, spelled out.
  const POP: Record<string, string> = { CVI: 'heart/vascular use', NVI: 'stroke/brain-vessel use', ACS: 'after a heart attack', PCI: 'after a stent', 'non-ACS': 'not after a heart attack', 'non-PCI': 'no stent' };
  const popLabel = (p: string) => p.split(/\s+/).map((w) => POP[w] ?? w).join(', ');
  const geneRow = (g: string) => genes.find((x) => x.gene === g);
  // Compact diplotype for the gauge: long allele names collapse to their bracketed short name ("…(HapB3)" -> "HapB3").
  const shortDip = (d: string | null) =>
    !d || d.length <= 14 ? (d ?? '') : d.split('/').map((a) => /\(([^()]+)\)\s*$/.exec(a)?.[1] ?? a.slice(0, 8) + '…').join('/');
</script>

<div class="grid">
  <div class="banner small"><span class="bi">{def.icon}</span>
    <div><strong>{def.title}</strong> — how your versions of the drug-processing genes change the way your body handles common medicines, and what the published prescribing guidelines say for someone with your genes.</div>
    {#if learnSlug}<button onclick={() => go('learn', learnSlug)}>📖 Background</button>{/if}
  </div>

  <Explain title="How to read this page (start here)" learn="Pharmacogenomics">
    <p>Most drugs are cleared or activated by a few enzymes, mainly in the liver. Your genes decide which <Term t="Star allele" label="versions" /> of those enzymes you make, and so how <em>fast</em> they work — your <Term t="Metaboliser phenotype" label="metaboliser type" />:</p>
    <ul>
      <li><strong>Poor / intermediate</strong> — the enzyme is slow. A drug it clears can build up (more side effects); a <Term t="Prodrug" label="pro-drug" /> it activates (codeine, clopidogrel) may not work.</li>
      <li><strong>Normal</strong> — what dosing labels assume.</li>
      <li><strong>Rapid / ultrarapid</strong> — the opposite: drugs are cleared too quickly, or pro-drugs over-activated.</li>
    </ul>
    <p>Expert panels (<Term t="CPIC" /> in the US, <Term t="DPWG" /> in the Netherlands) publish what a prescriber should do for each type. This page runs <Term t="PharmCAT" /> — the standard open tool, used by clinical labs — on your data and shows every guideline that matches you.</p>
    <p><strong>This is not a reason to start, stop or change a medicine on your own.</strong> It is something to show a doctor or pharmacist, ideally after confirming with a clinical pharmacogenomic test.</p>
  </Explain>

  {#if !loaded}
    <p class="muted small">Loading your medicine results…</p>
  {:else if error}
    <p class="small warn">{error}</p>
  {:else}
    <Card title="Your enzyme speeds" subtitle="Each row is a gene; the dot shows how fast your version of the enzyme works compared with the typical version">
      <div class="gauge">
        <div class="ghead"><span></span>{#each SPEED_LABEL as l, i}<span style:color={speedColor(i)}>{l}</span>{/each}</div>
        {#each gauged as g}
          {@const i = speedIndex(g)}
          <button class="grow" class:sel={openGene === g.gene} onclick={() => { openGene = openGene === g.gene ? null : g.gene; document.getElementById(`gene-${g.gene}`)?.scrollIntoView({ behavior: 'smooth', block: 'center' }); }}>
            <span class="gname mono">{g.gene}</span>
            {#each SPEED as _, j}
              <span class="cell" class:normal={j === 2}>{#if j === i}<span class="dot" style:background={speedColor(i)} title={`${g.diplotype} · ${g.phenotype}`}></span><span class="dip mono">{shortDip(g.diplotype)}</span>{/if}</span>
            {/each}
          </button>
        {/each}
      </div>
      <p class="small muted">
        {#if notable.length}You differ from “normal” in <strong>{notable.length}</strong> gene{notable.length > 1 ? 's' : ''}: {notable.map((g) => `${g.gene} (${g.phenotype.toLowerCase()})`).join(', ')}. That is typical — almost everyone has at least one.{:else}All typed enzymes are in the normal range.{/if}
        Click a row for details.
      </p>
      <Explain title="Why “intermediate” isn’t a disease">
        <p>These are normal human variations, not illnesses. A slow enzyme only matters when you take a drug that depends on it — then the guideline might suggest a different dose or a different drug. About 95% of people have at least one gene on this page that changes some guideline.</p>
      </Explain>
    </Card>

    <Card title="Your genes, one by one" subtitle="Diplotype (your two versions), what it means, and how sure we are">
      <div class="genes">
        {#each genes as g (g.gene)}
          {@const open = openGene === g.gene}
          {@const det = json(g.outside_detail)}
          {@const voi = json(g.variants_of_interest, [])}
          <div class="gene" id={`gene-${g.gene}`} class:open class:untyped={g.call_level === 'Not callable'}>
            <button class="ghd" onclick={() => (openGene = open ? null : g.gene)}>
              <span class="mono gsym">{g.gene}</span>
              <span class="mono dip2" title={g.diplotype}>{shortDip(g.diplotype) || '—'}</span>
              <span class="ph" title={g.phenotype}>{g.phenotype}</span>
              <span class="pill {OVERALL_PILL[g.overall]}" title="Overall confidence: the weaker of evidence and call">{g.overall}</span>
            </button>
            {#if open}
              <div class="gbody small">
                {#if aboutGene(g.gene)}
                  {@const ab = aboutGene(g.gene)!}
                  {#if ab.name}<p><strong>{g.gene}</strong> = <em>{ab.name}</em></p>{/if}
                  {#if ab.machine}
                    <details><summary class="faint">Technical summary (machine-generated)</summary><p class="faint">{ab.text}</p></details>
                  {:else if ab.text}<p>{firstPara(ab.text)}</p>{/if}
                  {#if !ab.machine && firstPara(ab.text) !== ab.text}<details><summary class="faint">More about {g.gene}</summary><p class="rec">{ab.text.slice(firstPara(ab.text).length).trim()}</p></details>{/if}
                  <p class="faint">From <a href={ab.url} target="_blank" rel="noreferrer">{ab.src} ↗</a>{ab.when ? ` · ${ab.when}` : ''}</p>
                {/if}
                {#if g.allele1_function || g.allele2_function}
                  <p>{#if g.allele2}Your two copies: <span class="mono">{g.allele1}</span> — {g.allele1_function ?? 'unknown function'}; <span class="mono">{g.allele2}</span> — {g.allele2_function ?? 'unknown function'}.{:else}Your version: <span class="mono">{g.allele1}</span> — {g.allele1_function ?? 'unknown function'}.{/if}
                    {#if g.activity_score && g.activity_score !== 'n/a'}<Term t="Activity score" /> <strong>{g.activity_score}</strong> (typical = 2).{/if}</p>
                {/if}
                <div class="row">
                  <span>Evidence <span class="pill {OVERALL_PILL[g.evidence_level]}">{g.evidence_level}</span></span>
                  <span>Your call <span class="pill {CALL_PILL[g.call_level]}">{g.call_level}</span></span>
                  <span class="faint">Called by {g.call_source ?? '—'}</span>
                </div>
                <ul class="reasons">
                  {#each [...(g.evidence_reasons ?? []), ...(g.call_reasons ?? [])] as r}<li>{r}</li>{/each}
                </ul>
                {#if g.positions_total}
                  <p class="faint">{g.positions_total} defining positions checked · {g.positions_variant} differ from the reference · {g.positions_missing} unreadable{#if g.missing_rsids?.length}{' '}({g.missing_rsids.slice(0, 6).join(', ')}{g.missing_rsids.length > 6 ? '…' : ''}){/if}.</p>
                {/if}
                {#if voi.length}
                  <p>Your variant positions: {#each voi as v, k}<span class="mono">{v.rsid ?? v.pos} {v.call}</span>{k < voi.length - 1 ? ' · ' : ''}{/each}</p>
                {/if}
                {#if g.gene === 'CYP2D6' && det}
                  <p class="note">Called from your reads by <Term t="Cyrius" /> (it counts gene copies and resolves the look-alike CYP2D7 pseudogene). Copies of CYP2D6 found: <strong>{det.total_cn != null ? det.total_cn - 2 : '?'}</strong> (2 is typical). {#if det.warnings?.length}Caution: {det.warnings.join('; ')}.{' '}{/if}Cyrius was validated on Illumina data; yours is DNBSEQ, so call confidence is capped at Medium.</p>
                {/if}
                {#if g.gene.startsWith('HLA') && det}
                  <p class="note">Typed from your reads by <Term t="T1K" /> against IPD-IMGT/HLA {det.imgt_hla}. PharmCAT only checks for the specific risk versions named in guidelines. Your result: {g.phenotype}.</p>
                {/if}
                {#if g.gene === 'MT-RNR1' && det}
                  <p class="note">Read from your mitochondrial variants ({(det.callable_fraction * 100).toFixed(0)}% of the gene readable). Result: {det.call}{det.function ? ` — ${det.function}` : ''}.{#if det.variants_in_gene?.length}{' '}Other variants in the gene (not linked to drug response): {det.variants_in_gene.map((v: any) => v.name).join(', ')}.{/if}</p>
                {/if}
                {#if g.messages?.length}
                  <details><summary class="faint">PharmCAT notes ({g.messages.length})</summary>{#each g.messages as msg}<p class="faint">{msg}</p>{/each}</details>
                {/if}
                {#if g.related_drugs?.length}
                  <p>Related drugs: {#each g.related_drugs.slice(0, 30) as d, k}<button class="link" onclick={() => { search = d; openDrug = d; }}>{d}</button>{k < Math.min(g.related_drugs.length, 30) - 1 ? ', ' : ''}{/each}{#if g.related_drugs.length > 30}<span class="faint">{' '}and {g.related_drugs.length - 30} more — search below.</span>{/if}</p>
                {/if}
                <button class="link" onclick={() => go('variants', g.gene)}>Browse your variants in {g.gene} →</button>
              </div>
            {/if}
          </div>
        {/each}
      </div>
      {#if untyped.length}
        <p class="small faint">Not typed: {untyped.map((g) => g.gene).join(', ')} — too many defining positions were unreadable.</p>
      {/if}
      <Explain title="Star alleles, diplotypes and why some genes need special tools" learn="Pharmacogenomics">
        <p>Pharmacogenes are described by named versions called <Term t="Star allele" label="star alleles" />: <span class="mono">*1</span> is usually the typical one, other numbers are specific combinations of variants. You have two copies, so your result is a pair — a <Term t="Diplotype" label="diplotype" /> such as <span class="mono">*1/*4</span>.</p>
        <p>Most genes can be read straight from your variant list. Three can’t: <strong>CYP2D6</strong> has a near-identical neighbour (CYP2D7) and is often duplicated or deleted, so it needs copy-number counting from the raw reads (<Term t="Cyrius" />); <strong>HLA</strong> genes have thousands of versions and need a dedicated typer (<Term t="T1K" />); <strong>MT-RNR1</strong> sits in mitochondrial DNA. These are passed to PharmCAT as <Term t="Outside call" label="outside calls" />.</p>
      </Explain>
    </Card>

    <Card title="Medicines" subtitle="Every published guideline that matches your genotypes">
      <div class="row search">
        <input type="search" placeholder="Search a drug or gene (e.g. codeine, CYP2C19)" bind:value={search} aria-label="Search drugs" />
        <span class="small faint">{actionDrugs.length} different · {noteDrugs.length} worth knowing · {standardDrugs.length} standard</span>
      </div>

      <h4>⚠️ Guidance differs for you <span class="faint small">({actionDrugs.length})</span></h4>
      <p class="small muted">A guideline recommends a different dose or a different drug for someone with your genes.</p>
      {@render drugList(actionDrugs, 'No matching drugs.')}

      <h4 class="mt">📝 Worth knowing <span class="faint small">({noteDrugs.length})</span></h4>
      <p class="small muted">No change by default, but the guideline adds advice — monitoring, what to do if the drug doesn’t work, or a dose calculator.</p>
      {@render drugList(noteDrugs, 'No matching drugs.')}

      <h4 class="mt"><button class="link big" onclick={() => (showStandard = !showStandard)}>{showStandard ? '▾' : '▸'} ✓ Standard guidance <span class="faint small">({standardDrugs.length})</span></button></h4>
      {#if showStandard || q}{@render drugList(standardDrugs, 'No matching drugs.')}{:else}<p class="small muted">Guidelines exist and say “use as normal” for your genotype. Reassuring, and still worth knowing.</p>{/if}

      {#if fdaOnly.length}
        <h4 class="mt"><button class="link big" onclick={() => (showFda = !showFda)}>{showFda ? '▾' : '▸'} 🏷️ Mentioned only on FDA labels <span class="faint small">({fdaOnly.length})</span></button></h4>
        {#if showFda || q}{@render drugList(fdaOnly, 'No matching drugs.')}{:else}<p class="small muted">The US drug label mentions one of these genes, but no CPIC or DPWG guideline exists. Weaker, and not graded.</p>{/if}
      {/if}

      {#if noneDrugs.length}
        <h4 class="mt"><button class="link big" onclick={() => (showNone = !showNone)}>{showNone ? '▾' : '▸'} ∅ Guideline exists, but not for your result <span class="faint small">({noneDrugs.length})</span></button></h4>
        {#if showNone || q}{@render drugList(noneDrugs, 'No matching drugs.')}{:else}<p class="small muted">A guideline or FDA label covers this drug, but only for genotypes you don’t have (for example a risk version you don’t carry, or G6PD deficiency).</p>{/if}
      {/if}

      <Explain title="Where the advice comes from, and how it is graded" learn="How this project grades">
        <p><strong>CPIC</strong> rates each recommendation Strong / Moderate / Optional; <strong>DPWG</strong> recommendations are graded Moderate here; FDA label text is shown but not graded. The pill on each drug is the <em>overall</em> grade — the weaker of the guideline’s strength and how sure we are of your genotype (e.g. CYP2D6 is capped at Moderate because its call is Medium confidence). Some drugs appear more than once because the guideline differs by situation (e.g. clopidogrel for heart stents vs. stroke).</p>
      </Explain>
    </Card>

    <div class="grid cols-2">
      <Card title="More detail" subtitle="The full PharmCAT report and the evidence behind it">
        <ul class="small">
          {#if reportUrl}<li><a href={reportUrl} target="_blank" rel="noreferrer">Open the full PharmCAT report ↗</a> — the same report a clinical lab would review, with every allele and citation.</li>{/if}
          <li>Knowledge used: PharmCAT {m.knowledge?.versions?.pharmcat ?? '—'} (CPIC/DPWG/FDA via ClinPGx){#if m.knowledge?.versions?.imgt_hla}, IPD-IMGT/HLA {m.knowledge.versions.imgt_hla}{/if}. Refreshed with <span class="mono">wgs knowledge refresh</span>; changes show on <button class="link inline" onclick={() => go('changes')}>What changed</button>.</li>
          <li>Your data is GRCh37; PharmCAT is built for GRCh38. Only the ~1,200 positions it checks are translated, each verified against the reference — see <span class="mono">pgx_positions</span> in the SQL tab.</li>
        </ul>
      </Card>
      <Card title="Compare with tellmeGen" subtitle="Coming at this checkpoint">
        <p class="small muted">Your tellmeGen download didn’t include their pharmacogenetics report. If you download it from your tellmeGen account and drop it in <span class="mono">WGS/05_Reports</span>, it can be compared gene by gene here.</p>
      </Card>
    </div>

    <Card title="Single-variant drug-response entries (ClinVar)" subtitle="Secondary: single variants ClinVar links to a drug response, many from older research studies. 0★ means a submitter gave no review criteria, so it is graded Limited — a lead, not a finding. The gene-level results above are what prescribers use.">
      <ClaimList {m} section="pgx" groups={['drug_response']} limit={10} empty="No ClinVar drug-response variants." />
    </Card>
  {/if}
</div>

{#snippet drugList(list: any[], empty: string)}
  {#if !list.length}
    <p class="small muted">{empty}</p>
  {:else}
    <div class="drugs">
      {#each list as d (d.drug)}
        {@const open = openDrug === d.drug}
        <div class="drug" class:open>
          <button class="dhd" onclick={() => (openDrug = open ? null : d.drug)}>
            <span class="dname">{d.drug}</span>
            <span class="mono faint small">{d.genes.join(' · ')}</span>
            <span class="spacer"></span>
            {#if d.differ}<span class="pill warn" title="The US (CPIC) and Dutch (DPWG) panels reach different conclusions for your genotype — open to compare. The bucket follows the more cautious one.">⚖️ guidelines differ</span>{/if}
            {#if d.notes.length}<span class="pill soon" title="Has a guideline note — open for details">📝 note</span>{/if}
            {#if d.overall}<span class="pill {OVERALL_PILL[d.overall]}">{d.overall}</span>{/if}
          </button>
          {#if open}
            <div class="dbody small">
              {#each d.rows as r}
                <div class="src">
                  <div class="row">
                    <span class="pill">{r.source}</span>
                    {#if r.classification && r.classification !== 'Unspecified'}<span class="faint">strength: {r.classification}</span>{/if}
                    {#if r.population && r.population !== 'general'}<span class="faint" title={r.population}>· {popLabel(r.population)}</span>{/if}
                    <span class="spacer"></span>
                    {#if r.url}<a href={r.url} target="_blank" rel="noreferrer">source ↗</a>{/if}
                  </div>
                  <p class="faint">For: {(r.phenotypes ?? []).join('; ')}</p>
                  <p class="rec">{plain(r.recommendation)}</p>
                  {#if r.implications?.length}<p class="faint">Why: {r.implications.map(plain).join(' ')}</p>{/if}
                </div>
              {/each}
              {#each d.notes as r}
                <div class="src">
                  <div class="row"><span class="pill">{r.source}</span><span class="faint">note — no single recommendation</span><span class="spacer"></span>{#if r.url}<a href={r.url} target="_blank" rel="noreferrer">source ↗</a>{/if}</div>
                  {#each r.messages as msg}<p class="rec">{plain(msg)}</p>{/each}
                </div>
              {/each}
              {#each d.unmatched as r}
                <div class="src">
                  <div class="row"><span class="pill">{r.source}</span><span class="faint">no recommendation for your genotype</span><span class="spacer"></span>{#if r.url}<a href={r.url} target="_blank" rel="noreferrer">source ↗</a>{/if}</div>
                </div>
              {/each}
              {#each d.genes as gname}
                {@const gr = geneRow(gname)}
                {#if gr}<button class="link" onclick={() => { openGene = gname; document.getElementById(`gene-${gname}`)?.scrollIntoView({ behavior: 'smooth', block: 'center' }); }}>Your {gname}: {gr.diplotype ?? '—'} ({gr.phenotype}) →</button><br />{/if}
              {/each}
            </div>
          {/if}
        </div>
      {/each}
    </div>
  {/if}
{/snippet}

<style>
  .banner { display: flex; gap: 14px; align-items: center; background: var(--panel); border: 1px solid var(--border); border-radius: 14px; padding: 12px 16px; }
  .banner > div { flex: 1; }
  .bi { font-size: 1.8rem; }
  .gauge { display: grid; gap: 2px; margin-bottom: 8px; }
  .ghead, .grow { display: grid; grid-template-columns: 90px repeat(5, minmax(0, 1fr)); align-items: center; }
  .ghead span { font-size: 0.75rem; font-weight: 600; text-align: center; }
  .grow { background: none; border: 0; padding: 0; color: var(--text); cursor: pointer; border-radius: 8px; text-align: left; }
  .grow:hover, .grow.sel { background: var(--panel-2); }
  .gname { font-weight: 600; font-size: 0.85rem; padding-left: 6px; }
  .cell { position: relative; height: 30px; display: flex; align-items: center; justify-content: center; gap: 6px; border-left: 1px dashed var(--border); }
  .cell.normal { background: color-mix(in srgb, var(--good) 7%, transparent); }
  .dot { width: 14px; height: 14px; border-radius: 50%; box-shadow: 0 0 0 3px color-mix(in srgb, var(--text) 10%, transparent); flex: none; }
  .dip { font-size: 0.72rem; color: var(--muted); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 70%; }
  .genes { display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 8px; align-items: start; }
  .gene { border: 1px solid var(--border); border-radius: 10px; background: var(--panel-2); }
  .gene.open { grid-column: 1 / -1; }
  .gene.untyped { opacity: 0.7; }
  .ghd { display: grid; grid-template-columns: 72px minmax(70px, auto) 1fr auto; gap: 8px; align-items: center; width: 100%; background: none; border: 0; padding: 8px 10px; color: var(--text); text-align: left; cursor: pointer; }
  .gsym { font-weight: 700; }
  .dip2 { font-size: 0.8rem; color: var(--muted); overflow: hidden; text-overflow: ellipsis; white-space: nowrap; max-width: 180px; }
  .ph { font-size: 0.85rem; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  .gbody { padding: 0 12px 10px; border-top: 1px solid var(--border); }
  .gbody p { margin: 8px 0; }
  .reasons { margin: 6px 0; padding-left: 18px; color: var(--muted); }
  .note { padding: 6px 10px; border-left: 3px solid var(--accent); background: color-mix(in srgb, var(--accent) 6%, transparent); border-radius: 4px; }
  .search input { flex: 1; min-width: 240px; }
  h4 { margin: 14px 0 4px; }
  h4.mt { margin-top: 18px; }
  .drugs { display: grid; gap: 4px; }
  .drug { border: 1px solid var(--border); border-radius: 8px; }
  .drug.open { background: var(--panel-2); }
  .dhd { display: flex; gap: 10px; align-items: center; width: 100%; background: none; border: 0; padding: 7px 10px; color: var(--text); cursor: pointer; text-align: left; }
  .dname { font-weight: 600; text-transform: capitalize; min-width: 150px; }
  .dbody { padding: 0 12px 10px; }
  .src { border-top: 1px dashed var(--border); padding-top: 8px; margin-top: 6px; }
  .src p { margin: 4px 0; }
  .rec { white-space: pre-line; }
  .link { background: none; border: 0; padding: 0; color: var(--accent); cursor: pointer; }
  .link.big { font: inherit; font-weight: 700; color: var(--text); }
  .link.inline { font-size: inherit; }
  @media (max-width: 700px) { .ghead, .grow { grid-template-columns: 70px repeat(5, minmax(0, 1fr)); } .dip { display: none; } }
</style>
