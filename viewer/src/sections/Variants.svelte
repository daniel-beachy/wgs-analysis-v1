<script lang="ts">
  import type { Manifest } from '../lib/data';
  import { query, useRelease, sqlString, type Row } from '../lib/db';
  import { int, compact } from '../lib/format';
  import { route, go } from '../lib/router.svelte';
  import Card from '../components/Card.svelte';
  import Term from '../components/Term.svelte';
  import Explain from '../components/Explain.svelte';
  import { IMPACT_PILL, IMPACT_WORD, CLASS_LABEL, CLASS_PILL, OVERALL_PILL, CALL_PILL, STAR_WORD, CONTINENTS,
    consequenceWord, consequenceLabel, aaShort, stars, pct, rarity, carriersOf, parseSources } from '../lib/evidence';

  let { m, sex = '' }: { m: Manifest; sex?: string } = $props();

  const CHROMS = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y', 'MT'];
  const PAGE = 50;

  let search = $state(route.rest[0] ?? '');
  let chrom = $state('');
  let vtype = $state('');
  let zyg = $state('');
  let filter = $state('PASS');
  let known = $state('');
  let minGq = $state(20);
  let impact = $state('');
  let clin = $state('');
  let freq = $state('');
  let geneHit = $state<Row | null | undefined>(undefined);
  const annotated = $derived(!!m.tables.annotations);
  const annFilter = $derived(!!(impact || clin || freq));
  const tbl = $derived(annotated && (filter === 'PASS' || annFilter) ? 'annotations' : 'variants');
  let detail = $state<{ v?: Row; a?: Row; claim?: Row }>({});
  let page = $state(0);

  let rows = $state<Row[]>([]);
  let total = $state<number | null>(null);
  let ms = $state(0);
  let error = $state('');
  let busy = $state(false);
  let selected = $state<Row | null>(null);
  let callable = $state<string>('');
  let searchNote = $state('');

  // Parse the free-text box: rsID, chr:pos, chr:start-end.
  function parseSearch(s: string): { where: string[]; note: string } {
    const t = s.trim().replace(/,/g, '');
    if (!t) return { where: [], note: '' };
    let mm = /^rs\d+$/i.exec(t);
    if (mm) return { where: [`rsid = ${sqlString(t.toLowerCase())}`], note: '' };
    mm = /^(?:chr)?([0-9]{1,2}|x|y|mt?)[:\s]+(\d+)(?:\s*-\s*(\d+))?$/i.exec(t);
    if (mm) {
      let c = mm[1].toUpperCase();
      if (c === 'M') c = 'MT';
      const a = +mm[2], b = mm[3] ? +mm[3] : +mm[2];
      return { where: [`chrom = ${sqlString(c)}`, `pos BETWEEN ${Math.min(a, b)} AND ${Math.max(a, b)}`], note: '' };
    }
    if (GENE_RE.test(t) && m.tables.genes) {
      if (geneHit === undefined) return { where: ['false'], note: '' };
      if (geneHit) return { where: [`chrom = ${sqlString(geneHit.chrom)}`, `pos BETWEEN ${geneHit.start} AND ${geneHit.end}`], note: '' };
      return { where: ['false'], note: `No gene called “${t}” in Ensembl GRCh37. Gene symbols look like LCT, MTHFR, BRCA2.` };
    }
    return { where: ['false'], note: `Couldn’t read “${s}”. Try an rsID (rs4988235), a gene (LCT), a position (2:136608646) or a range (chr2:136,500,000-136,700,000).` };
  }
  const GENE_RE = /^[A-Za-z][A-Za-z0-9.-]{1,19}$/;

  // Resolve a gene symbol to its coordinates (small `genes` table), so the variant query can use a range.
  $effect(() => {
    const t = search.trim();
    if (!GENE_RE.test(t) || /^rs\d+$/i.test(t) || !m.tables.genes) { geneHit = undefined; return; }
    geneHit = undefined;
    (async () => {
      await useRelease(m);
      const r = await query(`SELECT * FROM genes WHERE upper(symbol) = ${sqlString(t.toUpperCase())} ORDER BY biotype = 'protein_coding' DESC LIMIT 1`);
      if (search.trim() === t) geneHit = r.rows[0] ?? null;
    })();
  });

  const where = $derived.by(() => {
    const p = parseSearch(search);
    const w = [...p.where];
    if (chrom) w.push(`chrom = ${sqlString(chrom)}`);
    if (vtype) w.push(`vtype = ${sqlString(vtype)}`);
    if (zyg) w.push(`zygosity = ${sqlString(zyg)}`);
    if (filter) w.push(`filter = ${sqlString(filter)}`);
    if (known === 'known') w.push('rsid IS NOT NULL');
    if (known === 'novel') w.push('rsid IS NULL');
    if (minGq > 0) w.push(`gq >= ${Math.floor(minGq)}`);
    if (tbl === 'annotations') {
      if (impact) w.push(`impact = ${sqlString(impact)}`);
      if (clin === 'any') w.push('clinvar_class IS NOT NULL');
      else if (clin === 'plp') w.push("clinvar_class IN ('pathogenic', 'likely_pathogenic')");
      else if (clin === 'risk') w.push("clinvar_class IN ('risk_factor', 'association', 'protective')");
      else if (clin === 'vus') w.push("clinvar_class IN ('uncertain', 'conflicting')");
      else if (clin) w.push(`clinvar_class = ${sqlString(clin)}`);
      if (freq === 'rare') w.push('coalesce(af_1kg, 0) < 0.01');
      if (freq === 'absent') w.push('af_1kg IS NULL');
      if (freq === 'common') w.push('af_1kg >= 0.05');
    }
    return { sql: w.length ? 'WHERE ' + w.join(' AND ') : '', note: p.note };
  });

  let seq = 0;
  $effect(() => {
    const w = where.sql;
    const pg = page;
    const from = tbl;
    searchNote = where.note;
    const my = ++seq;
    busy = true;
    error = '';
    (async () => {
      try {
        await useRelease(m);
        const r = await query(`SELECT * FROM ${from} ${w} ORDER BY chrom_order, pos LIMIT ${PAGE} OFFSET ${pg * PAGE}`);
        if (my !== seq) return;
        rows = r.rows;
        ms = r.ms;
        busy = false;
        const c = await query(`SELECT count(*) AS n FROM ${from} ${w}`);
        if (my === seq) total = c.rows[0].n;
      } catch (e) {
        if (my === seq) { error = String(e); busy = false; }
      }
    })();
  });

  function reset() { page = 0; selected = null; total = null; }

  async function pick(r: Row) {
    selected = r;
    callable = '';
    detail = {};
    const key = `chrom = ${sqlString(r.chrom)} AND pos = ${r.pos} AND alt = ${sqlString(r.alt)}`;
    (async () => {
      const v = (await query(`SELECT * FROM variants WHERE ${key} LIMIT 1`)).rows[0];
      const a = annotated ? (await query(`SELECT * FROM annotations WHERE ${key} LIMIT 1`)).rows[0] : undefined;
      const claim = m.tables.claims && a?.clinvar_class ? (await query(`SELECT * FROM claims WHERE claim_id = ${sqlString(`variant:${r.chrom}-${r.pos}-${r.ref}-${r.alt}:clinvar`)}`)).rows[0] : undefined;
      if (selected === r) detail = { v, a, claim };
    })();
    if (!m.tables.callable) return;
    const c = await query(`SELECT state FROM callable WHERE chrom = ${sqlString(r.chrom)} AND start < ${r.pos} AND "end" >= ${r.pos} LIMIT 1`);
    callable = c.rows[0]?.state ?? 'unknown';
  }

  const zygWord = (r: Row) => {
    if (r.zygosity === 'het') return 'one of your two copies';
    if (r.zygosity === 'hom') return (r.chrom === 'X' || r.chrom === 'Y') && sex === 'XY' && r.chrom !== 'MT' ? 'your single copy (hemizygous)' : r.chrom === 'MT' ? 'your mitochondrial DNA' : 'both of your copies';
    return 'neither copy (reference)';
  };
  const typeWord: Record<string, string> = { SNV: 'a single-letter change', INS: 'an insertion', DEL: 'a deletion', MNV: 'a multi-letter change' };
  const fmtAllele = (s: string) => (s.length > 14 ? s.slice(0, 12) + '…' : s);
  const gnomad = (r: Row) => `https://gnomad.broadinstitute.org/variant/${r.chrom}-${r.pos}-${r.ref}-${r.alt}?dataset=gnomad_r2_1`;
  const ucsc = (r: Row) => `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr${r.chrom === 'MT' ? 'M' : r.chrom}:${r.pos - 25}-${r.pos + 25}`;
  const callableWord: Record<string, string> = {
    CALLABLE: 'Well covered (≥ 10 reads): a confident call.', LOW_5_9: 'Thin coverage (5–9 reads): treat with some caution.',
    LOW_1_4: 'Very thin coverage (1–4 reads): unreliable.', NO_COVERAGE: 'No reads here.', HIGH: 'Unusually deep coverage — often a repetitive region where reads pile up wrongly.',
  };
</script>

<div class="grid">
  <Card title="Variant explorer" subtitle="Every position where your genome differs from the reference — about 4.6 million confident calls.">
    <Explain learn="Variants: how your genome differs">
      <p>Each row is one <Term t="Variant" label="variant" />: a position (<em>chromosome : position</em>), the <Term t="REF" label="reference letter(s)" /> and <Term t="ALT" label="your letter(s)" />, and your <Term t="Genotype" label="genotype" /> — whether you carry it on one copy (<Term t="Heterozygous" label="het" />) or both (<Term t="Homozygous" label="hom" />). Most variants are common and harmless; this table shows the raw evidence, and later sections decide which ones mean something. Click a row for a plain-language reading.</p>
    </Explain>
    <div class="filters">
      <input class="search" placeholder="LCT · rs4988235 · 2:136608646 · chr2:136,500,000-136,700,000" bind:value={search} oninput={reset} aria-label="Search" />
      <select bind:value={chrom} onchange={reset} aria-label="Chromosome"><option value="">All chromosomes</option>{#each CHROMS as c}<option value={c}>chr{c}</option>{/each}</select>
      <select bind:value={vtype} onchange={reset} aria-label="Type"><option value="">All types</option><option value="SNV">SNV</option><option value="INS">Insertion</option><option value="DEL">Deletion</option><option value="MNV">MNV</option></select>
      <select bind:value={zyg} onchange={reset} aria-label="Zygosity"><option value="">Het + hom</option><option value="het">Heterozygous (one copy)</option><option value="hom">Homozygous (both copies)</option></select>
      <select bind:value={filter} onchange={reset} aria-label="Call"><option value="PASS">Confident (PASS)</option><option value="RefCall">RefCall (matches reference)</option><option value="NoCall">NoCall (undecided)</option><option value="">All records</option></select>
      <select bind:value={known} onchange={reset} aria-label="Known"><option value="">Known + novel</option><option value="known">With rsID (catalogued)</option><option value="novel">No rsID</option></select>
      {#if annotated}
        <select bind:value={impact} onchange={reset} aria-label="Impact"><option value="">Any effect</option><option value="HIGH">High impact (likely breaks protein)</option><option value="MODERATE">Moderate (changes protein)</option><option value="LOW">Low (silent / splice region)</option><option value="MODIFIER">Modifier (non-coding)</option></select>
        <select bind:value={clin} onchange={reset} aria-label="ClinVar"><option value="">ClinVar: any or none</option><option value="any">In ClinVar</option><option value="plp">Pathogenic / likely pathogenic</option><option value="vus">Uncertain / conflicting</option><option value="drug_response">Drug response</option><option value="risk">Risk factor / association</option><option value="benign">Benign</option></select>
        <select bind:value={freq} onchange={reset} aria-label="Frequency"><option value="">Any frequency</option><option value="rare">Rare (&lt; 1%)</option><option value="absent">Not in 1000 Genomes</option><option value="common">Common (≥ 5%)</option></select>
      {/if}
      <label class="small muted gq">min <Term t="GQ" /> <input type="number" min="0" max="99" bind:value={minGq} oninput={reset} /></label>
    </div>
    {#if searchNote}<p class="small warn">{searchNote}</p>{/if}
    <div class="row small muted status">
      <span>{total == null ? 'Counting…' : `${int(total)} matching records`}</span>
      <span class="faint">· page {page + 1}{total != null ? ` of ${Math.max(1, Math.ceil(total / PAGE))}` : ''} · {ms.toFixed(0)} ms</span>
      <span class="spacer"></span>
      <button disabled={page === 0} onclick={() => { page--; selected = null; }}>← Prev</button>
      <button disabled={total != null && (page + 1) * PAGE >= total} onclick={() => { page++; selected = null; }}>Next →</button>
    </div>
    {#if error}<p class="warn small">{error}</p>{/if}
    <div class="tablewrap" class:busy>
      <table class="small">
        {#if tbl === 'annotations'}
          <thead><tr><th>Position</th><th>rsID</th><th>Gene</th><th>Effect</th><th>Ref → Alt</th><th>Genotype</th><th class="num"><Term t="Allele frequency" label="Frequency" /></th><th>ClinVar</th><th class="num"><Term t="GQ" /></th></tr></thead>
          <tbody>
            {#each rows as r (r.chrom + r.pos + r.alt)}
              <tr class:sel={selected === r} onclick={() => pick(r)}>
                <td class="mono">{r.chrom}:{int(r.pos)}</td>
                <td class="mono">{r.rsid ?? ''}</td>
                <td><strong>{r.gene ?? ''}</strong></td>
                <td>{#if r.impact && r.impact !== 'MODIFIER'}<span class="pill {IMPACT_PILL[r.impact]}">{consequenceLabel(r.consequence)}</span>{:else}<span class="faint">{consequenceLabel(r.consequence)}</span>{/if}
                  {#if r.aa_change}<span class="mono faint"> {aaShort(r.aa_change)}</span>{/if}</td>
                <td class="mono">{fmtAllele(r.ref)} → {fmtAllele(r.alt)}</td>
                <td>{r.zygosity}</td>
                <td class="num" title={rarity(r.af_1kg)}>{pct(r.af_1kg)}</td>
                <td>{#if r.clinvar_class}<span class="pill {CLASS_PILL[r.clinvar_class] ?? ''}">{CLASS_LABEL[r.clinvar_class] ?? r.clinvar_class}</span>{/if}</td>
                <td class="num">{r.gq ?? ''}</td>
              </tr>
            {:else}
              <tr><td colspan="9" class="muted">{busy ? 'Searching…' : 'No variants match these filters.'}</td></tr>
            {/each}
          </tbody>
        {:else}
          <thead><tr><th>Position</th><th>rsID</th><th>Ref → Alt</th><th>Type</th><th>Genotype</th><th class="num"><Term t="DP" label="Depth" /></th><th class="num"><Term t="AD" label="Reads ref / alt" /></th><th class="num"><Term t="GQ" /></th><th>Call</th></tr></thead>
          <tbody>
            {#each rows as r (r.chrom + r.pos + r.alt)}
              <tr class:sel={selected === r} onclick={() => pick(r)}>
                <td class="mono">{r.chrom}:{int(r.pos)}</td>
                <td class="mono">{r.rsid ?? ''}</td>
                <td class="mono">{fmtAllele(r.ref)} → {fmtAllele(r.alt)}</td>
                <td>{r.vtype}</td>
                <td><span class="mono faint">{r.gt}</span> {r.zygosity}</td>
                <td class="num">{r.dp ?? ''}</td>
                <td class="num">{r.ad_ref ?? ''} / {r.ad_alt ?? ''}</td>
                <td class="num">{r.gq ?? ''}</td>
                <td><span class="pill {r.filter === 'PASS' ? 'pass' : ''}">{r.filter}</span></td>
              </tr>
            {:else}
              <tr><td colspan="9" class="muted">{busy ? 'Searching…' : 'No records match these filters.'}</td></tr>
            {/each}
          </tbody>
        {/if}
      </table>
    </div>
    {#if geneHit}<p class="small faint">Showing {geneHit.symbol} ({geneHit.description?.replace(/ \[Source.*$/, '') ?? geneHit.biotype}) · chr{geneHit.chrom}:{int(geneHit.start)}–{int(geneHit.end)}</p>{/if}
  </Card>

  {#if selected}
    {@const r = { ...selected, ...(detail.v ?? {}), ...(detail.a ?? {}) }}
    {@const a = detail.a}
    <Card title={`chr${r.chrom}:${int(r.pos)} ${fmtAllele(r.ref)}→${fmtAllele(r.alt)}`} subtitle={r.rsid ? `Catalogued in dbSNP as ${r.rsid}` : 'Not catalogued with an rsID in the provider file'}>
      {#snippet aside()}<button onclick={() => (selected = null)} aria-label="Close">✕</button>{/snippet}
      <p class="reading">
        At position <strong>{int(r.pos)}</strong> on chromosome <strong>{r.chrom}</strong>, the reference genome has <code>{fmtAllele(r.ref)}</code>.
        You carry <code>{fmtAllele(r.alt)}</code> — {typeWord[r.vtype] ?? 'a change'} — on <strong>{zygWord(r)}</strong>.
        {#if r.dp}{r.ad_alt} of {r.dp} reads showed your version{r.vaf != null ? ` (${Math.round(r.vaf * 100)}%)` : ''}{#if r.zygosity === 'het'}, close to the ~50% expected for one copy{/if}.{/if}
        {#if r.filter !== 'PASS'} DeepVariant labelled this <strong>{r.filter}</strong>, so it is not treated as a confident variant.{/if}
      </p>
      {#if r.dp}
        <div class="reads" title="Reads supporting reference vs your allele">
          <div class="ref" style:flex={r.ad_ref || 0.0001}>{r.ad_ref} ref</div>
          <div class="alt" style:flex={r.ad_alt || 0.0001}>{r.ad_alt} alt</div>
        </div>
      {/if}
      <div class="grid cols-3 small" style="margin-top: 10px">
        <div><div class="faint">Genotype quality</div><strong>{r.gq ?? '—'}</strong> <span class="faint">{r.gq != null ? `(≈${(100 - 100 * Math.pow(10, -r.gq / 10)).toFixed(r.gq >= 30 ? 2 : 1)}% confident)` : ''}</span></div>
        <div><div class="faint">Coverage at this spot</div>{callable ? (callableWord[callable] ?? callable) : '…'}</div>
        <div><div class="faint">Multi-allelic site</div>{r.multiallelic ? 'Yes — more than one alternative letter here' : 'No'}</div>
      </div>
      <div class="row" style="margin-top: 12px">
        <span class="faint small">Look it up:</span>
        {#if r.rsid}<a href={`https://www.ncbi.nlm.nih.gov/snp/${r.rsid}`} target="_blank" rel="noreferrer">dbSNP</a>{/if}
        <a href={gnomad(r)} target="_blank" rel="noreferrer">gnomAD (population frequency)</a>
        <a href={`https://www.ncbi.nlm.nih.gov/clinvar/?term=${r.rsid ?? `${r.chrom}[chr] AND ${r.pos}[chrpos37]`}`} target="_blank" rel="noreferrer">ClinVar</a>
        <a href={ucsc(r)} target="_blank" rel="noreferrer">UCSC browser (hg19)</a>
      </div>
      {#if a}
        <div class="ann">
          <h4>What it does</h4>
          {#if a.gene}
            <p>It is in the gene <strong>{a.gene}</strong>{a.genes && a.genes.length > 1 ? ` (also overlaps ${a.genes.filter((g: string) => g !== a.gene).join(', ')})` : ''}, where it {consequenceWord(a.consequence)}{a.aa_change ? ` — amino acid ${aaShort(a.aa_change)}` : ''}.
              {#if a.impact}<span class="pill {IMPACT_PILL[a.impact]}"><Term t="Impact" label={`${a.impact.toLowerCase()} impact`} /></span> <span class="faint small">({IMPACT_WORD[a.impact]})</span>{/if}</p>
          {:else}<p>It {consequenceWord(a.consequence)} — most such variants have no known effect.</p>{/if}

          <h4>How common it is</h4>
          {#if a.af_1kg != null}
            <p>{pct(a.af_1kg)} of chromosomes in 1000 Genomes carry it — <strong>{rarity(a.af_1kg)}</strong>. About {pct(carriersOf(a.af_1kg))} of people have at least one copy{a.zygosity === 'hom' && r.chrom !== 'X' && r.chrom !== 'Y' && r.chrom !== 'MT' ? `, and about ${pct(a.af_1kg ** 2)} have two like you` : ''}.</p>
            <div class="afbars">
              {#each CONTINENTS as [k, label]}
                {@const f = a[`af_1kg_${k}`]}
                <div class="afrow"><span class="faint">{label}</span><div class="track"><div class="fill" style:width={`${Math.max(f ?? 0, 0.004) * 100}%`}></div></div><span class="num">{pct(f)}</span></div>
              {/each}
            </div>
          {:else}
            <p>Not seen in 1000 Genomes (2,504 people){a.gnomad_popmax_af == null ? ' — either very rare, or in a region that panel couldn’t measure.' : '.'}</p>
          {/if}
          {#if a.gnomad_popmax_af != null}<p class="small faint">gnomAD (≈140,000 people): highest frequency in any one population is {pct(a.gnomad_popmax_af)}.</p>{/if}

          <h4>What ClinVar says</h4>
          {#if a.clinvar_class}
            <p><span class="pill {CLASS_PILL[a.clinvar_class] ?? ''}">{CLASS_LABEL[a.clinvar_class] ?? a.clinvar_class}</span>
              <span class="stars" title={STAR_WORD[a.clinvar_stars]}>{stars(a.clinvar_stars)}</span> <span class="faint small">{STAR_WORD[a.clinvar_stars] ?? ''}</span>
              {#if a.clinvar_conditions}<br /><span class="small">{a.clinvar_conditions.split('|').filter((c: string) => !/^not (provided|specified)$/i.test(c)).slice(0, 6).join(' · ')}</span>{/if}
              <br /><a class="small" href={`https://www.ncbi.nlm.nih.gov/clinvar/variation/${a.clinvar_id}/`} target="_blank" rel="noreferrer">ClinVar record {a.clinvar_id} ↗</a></p>
          {:else}<p class="faint">No ClinVar entry for this exact change — the usual situation; ClinVar covers about 4 million of the hundreds of millions of known variants.</p>{/if}

          {#if detail.claim}
            {@const c = detail.claim}
            <div class="claim">
              <div class="row"><strong>Graded finding</strong>
                <span class="pill {OVERALL_PILL[c.evidence_level]}">evidence: {c.evidence_level}</span>
                <span class="pill {CALL_PILL[c.call_level]}">your call: {c.call_level}</span>
                <span class="pill {OVERALL_PILL[c.overall]}">overall: {c.overall}</span></div>
              <p>{c.statement}</p>
              <ul class="small">{#each [...(c.evidence_reasons ?? []), ...(c.call_reasons ?? [])] as x}<li>{x}</li>{/each}</ul>
              <p class="small faint">Sources: {parseSources(c).map((s) => `${s.source} ${s.version ?? ''}`).join(' · ')}</p>
            </div>
          {/if}
        </div>
      {:else if annotated}
        <p class="small faint" style="margin-top: 10px">{detail.v ? 'This record carries no alternative allele (or was not annotated), so there is nothing to interpret.' : 'Loading annotation…'}</p>
      {:else}
        <p class="small faint" style="margin-top: 10px">Run <code>wgs knowledge refresh</code> and re-run the pipeline to add gene, frequency and ClinVar evidence to every variant.</p>
      {/if}
    </Card>
  {/if}

  <p class="small faint">Tip: power users can query these tables directly in the <a href="#/sql" onclick={(e) => { e.preventDefault(); go('sql'); }}>SQL tab</a>.</p>
</div>

<style>
  .filters { display: flex; flex-wrap: wrap; gap: 8px; margin: 10px 0; }
  .search { flex: 1 1 320px; }
  .gq input { width: 64px; margin-left: 4px; }
  .status { margin: 4px 0 8px; }
  .tablewrap { overflow-x: auto; transition: opacity 0.15s; }
  .tablewrap.busy { opacity: 0.55; }
  tbody tr { cursor: pointer; }
  tbody tr:hover { background: var(--panel-2); }
  tr.sel { background: color-mix(in srgb, var(--accent) 14%, transparent) !important; }
  .warn { color: var(--warn); }
  .reading { font-size: 1rem; max-width: 80ch; }
  .reads { display: flex; height: 26px; border-radius: 6px; overflow: hidden; font-size: 0.78rem; font-weight: 600; max-width: 520px; }
  .reads > div { display: flex; align-items: center; justify-content: center; min-width: 46px; color: #04121f; }
  .ref { background: #94a3b8; }
  .alt { background: var(--accent-2); }
  button:disabled { opacity: 0.4; cursor: default; }
  .ann h4 { margin: 16px 0 4px; font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.05em; color: var(--muted); }
  .ann p { max-width: 85ch; }
  .afbars { display: grid; gap: 4px; max-width: 520px; font-size: 0.82rem; }
  .afrow { display: grid; grid-template-columns: 130px 1fr 56px; gap: 8px; align-items: center; }
  .track { height: 8px; background: var(--chip); border-radius: 4px; overflow: hidden; }
  .fill { height: 100%; background: var(--accent); }
  .stars { color: var(--warn); letter-spacing: 1px; }
  .claim { margin-top: 12px; padding: 10px 14px; border: 1px solid var(--border); border-radius: 10px; background: var(--panel-2); }
  .claim ul { margin: 4px 0; padding-left: 18px; }
</style>
