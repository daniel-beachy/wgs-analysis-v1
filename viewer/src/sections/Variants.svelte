<script lang="ts">
  import type { Manifest } from '../lib/data';
  import { query, useRelease, sqlString, type Row } from '../lib/db';
  import { int, compact } from '../lib/format';
  import { route, go } from '../lib/router.svelte';
  import Card from '../components/Card.svelte';
  import Term from '../components/Term.svelte';
  import Explain from '../components/Explain.svelte';

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
    return { where: ['false'], note: `Couldn’t read “${s}”. Try an rsID (rs4988235), a position (2:136608646) or a range (chr2:136,500,000-136,700,000). Gene names arrive with the knowledge layer.` };
  }

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
    return { sql: w.length ? 'WHERE ' + w.join(' AND ') : '', note: p.note };
  });

  let seq = 0;
  $effect(() => {
    const w = where.sql;
    const pg = page;
    searchNote = where.note;
    const my = ++seq;
    busy = true;
    error = '';
    (async () => {
      try {
        await useRelease(m);
        const r = await query(`SELECT * FROM variants ${w} ORDER BY chrom_order, pos LIMIT ${PAGE} OFFSET ${pg * PAGE}`);
        if (my !== seq) return;
        rows = r.rows;
        ms = r.ms;
        busy = false;
        const c = await query(`SELECT count(*) AS n FROM variants ${w}`);
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
      <input class="search" placeholder="rs4988235 · 2:136608646 · chr2:136,500,000-136,700,000" bind:value={search} oninput={reset} aria-label="Search" />
      <select bind:value={chrom} onchange={reset} aria-label="Chromosome"><option value="">All chromosomes</option>{#each CHROMS as c}<option value={c}>chr{c}</option>{/each}</select>
      <select bind:value={vtype} onchange={reset} aria-label="Type"><option value="">All types</option><option value="SNV">SNV</option><option value="INS">Insertion</option><option value="DEL">Deletion</option><option value="MNV">MNV</option></select>
      <select bind:value={zyg} onchange={reset} aria-label="Zygosity"><option value="">Het + hom</option><option value="het">Heterozygous (one copy)</option><option value="hom">Homozygous (both copies)</option></select>
      <select bind:value={filter} onchange={reset} aria-label="Call"><option value="PASS">Confident (PASS)</option><option value="RefCall">RefCall (matches reference)</option><option value="NoCall">NoCall (undecided)</option><option value="">All records</option></select>
      <select bind:value={known} onchange={reset} aria-label="Known"><option value="">Known + novel</option><option value="known">With rsID (catalogued)</option><option value="novel">No rsID</option></select>
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
      </table>
    </div>
  </Card>

  {#if selected}
    {@const r = selected}
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
      <p class="small faint" style="margin-top: 10px">Gene, consequence, population frequency and clinical evidence for every variant are added by the knowledge layer (Step 3); until then this view shows only what your own data says.</p>
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
</style>
