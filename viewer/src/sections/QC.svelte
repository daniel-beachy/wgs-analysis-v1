<script lang="ts">
  import type { Manifest } from '../lib/data';
  import { loadDoc, loadAudit } from '../lib/data';
  import { query, useRelease } from '../lib/db';
  import { int, pct, compact } from '../lib/format';
  import Card from '../components/Card.svelte';
  import Term from '../components/Term.svelte';
  import Explain from '../components/Explain.svelte';
  import RangeBar from '../components/RangeBar.svelte';
  import Bars from '../components/Bars.svelte';
  import Line from '../components/Line.svelte';
  import GenomeDepth from '../components/GenomeDepth.svelte';
  import { RECOUNT, DIFFICULT, difficultText } from '../lib/evidence';

  let { m }: { m: Manifest } = $props();

  type Check = { metric: string; value: number; low: number; high: number; status: string; explanation: string };

  // Plain-language framing for each automated check.
  const META: Record<string, { name: string; term: string; fmt: (v: number) => string; why: string }> = {
    snv_count: { name: 'Single-letter variants', term: 'SNV', fmt: compact, why: 'Everyone carries 3.5–4.5 million single-letter differences from the reference. Far fewer suggests missing data; far more suggests noise.' },
    indel_count: { name: 'Small insertions / deletions', term: 'Indel', fmt: compact, why: 'Short added or missing stretches of DNA. Their count depends on the caller, but should be in this range at 30×.' },
    titv: { name: 'Ti/Tv ratio', term: 'Ti/Tv', fmt: (v) => v.toFixed(2), why: 'Real mutations favour "transitions" about 2:1. Random sequencing errors don’t, so a ratio near 2.0–2.1 means the calls are real biology, not noise.' },
    het_hom_ratio: { name: 'Het / hom ratio', term: 'Heterozygous', fmt: (v) => v.toFixed(2), why: 'How many variants are on one copy versus both copies. It reflects ancestry and would shift if the sample were contaminated.' },
    het_vaf_median: { name: 'Allele balance', term: 'VAF', fmt: (v) => v.toFixed(2), why: 'At a heterozygous site about half the reads should show each version. A median near 0.50 means one clean sample.' },
    autosomal_mean_depth: { name: 'Average depth', term: 'Coverage', fmt: (v) => `${v.toFixed(1)}×`, why: 'How many reads cover a typical position. You ordered 30×; more reads means more confident calls.' },
    fraction_20x: { name: 'Genome at ≥ 20×', term: 'Coverage', fmt: (v) => pct(v), why: 'Share of the sequenceable genome with enough reads to call a genotype confidently.' },
    duplicate_pct: { name: 'Duplicate reads', term: 'Duplicate reads', fmt: (v) => `${v.toFixed(2)}%`, why: 'Copies of the same DNA fragment add no new evidence. Low is good.' },
    mapped_pct: { name: 'Reads placed on the map', term: 'Alignment', fmt: (v) => `${v.toFixed(2)}%`, why: 'Share of reads the aligner could place on the reference genome. Near 100% means a clean human sample.' },
    error_rate: { name: 'Mismatch rate', term: 'Phred score', fmt: (v) => `${(v * 100).toFixed(2)}%`, why: 'How often a read letter disagrees with the reference (includes your real variants). Under 1% is healthy.' },
    concordance_pct: { name: 'Agreement with tellmeGen file', term: 'Concordance', fmt: (v) => `${v.toFixed(2)}%`, why: 'The variant file and the provider’s genotype file were made separately; they should agree almost perfectly.' },
  };

  const SRC: Record<string, string> = { clinvar: 'ClinVar records (classification and review stars)', gwas: 'GWAS Catalog associations (study, effect allele)',
    pgs: 'PGS Catalog scores (trait, publication)', pubmed: 'PubMed citations (title)', links: 'Source links (sampled)', curated: 'Cited baseline risks (quote still on the page)' };
  const audit = $derived(loadAudit(m));
  const load = $derived(Promise.all([
    loadDoc(m, 'qc'), loadDoc(m, 'coverage'), loadDoc(m, 'alignment_stats'), loadDoc(m, 'fastq_stats'),
    loadDoc(m, 'reference_full'), loadDoc(m, 'inventory'), loadDoc(m, 'sections'), loadDoc(m, 'checks'),
    loadDoc(m, 'benchmark'), loadDoc(m, 'secondcall'),
  ]));

  let bins = $state<{ chrom: string; start: number; depth: number }[]>([]);
  $effect(() => {
    if (!m.tables.coverage_bins) return;
    useRelease(m).then(() => query('SELECT chrom, start, depth FROM coverage_bins WHERE chrom_order < 22 ORDER BY chrom_order, start'))
      .then((r) => (bins = r.rows as any));
  });

  type Calib = { kind: string; band: number; panel_af_band: string; sites: number; hom_ref_pct: number; het_pct: number; hom_alt_pct: number; missing_pct: number; your_af: number; panel_af: number };
  let calib = $state<Calib[]>([]);
  $effect(() => {
    if (!m.tables.pgs_genotype_qc) return;
    useRelease(m).then(() => query("SELECT * FROM pgs_genotype_qc ORDER BY kind = 'indel', band"))
      .then((r) => (calib = r.rows as any));
  });

  // Read re-count (ADR-022): every position a finding or a medicine result rests on, counted again from the CRAM.
  let rc = $state<any[]>([]);
  $effect(() => {
    if (!m.tables.recount) return;
    useRelease(m).then(() => query('SELECT * FROM recount')).then((r) => (rc = r.rows as any));
  });
  const RC_ORDER = ['confirmed', 'weak', 'unclear', 'disagrees', 'not_called', 'not_checkable'];
  const RC_SOURCE: Record<string, string> = { claim: 'Behind your findings', pgx: 'Defining your medicine results' };
  const RC_FLAG: Record<string, string> = { look_alike: 'many reads could belong to look-alike DNA elsewhere', one_strand: 'variant reads come from one reading direction only' };
  const copiesWord = (n: number | null | undefined) => (n == null ? '—' : ['no copies', 'one copy', 'two copies'][n] ?? `${n} copies`);
  const shortAllele = (a: string) => (a.length > 6 ? `${a.slice(0, 4)}…(${a.length})` : a);
  const RC_PRIORITY = ['disagrees', 'unclear', 'weak', 'not_called', 'not_checkable', 'confirmed'];
  const rcExceptions = $derived(rc.filter((r) => r.status !== 'confirmed' || (r.flags ?? []).length)
    .sort((a, b) => RC_PRIORITY.indexOf(a.status) - RC_PRIORITY.indexOf(b.status)));

  // Accuracy benchmark on Genome in a Bottle's HG002, and a second variant caller on your own reads (ADR-023).
  let sc = $state<any[]>([]);
  $effect(() => {
    if (!m.tables.second_caller) return;
    useRelease(m).then(() => query('SELECT * FROM second_caller ORDER BY TRY_CAST(chrom AS INTEGER), pos')).then((r) => (sc = r.rows as any));
  });
  const KIND: Record<string, string> = { all: 'All changes', snv: 'Single letters', indel: 'Insertions / deletions' };
  const REGION: Record<string, string> = { all: 'everywhere', easy: 'ordinary', difficult: 'hard to sequence' };
  const SC_VERDICT: Record<string, { word: string; pill: string }> = {
    reads_support: { word: 'reads show it', pill: 'pass' },
    reads_flagged: { word: 'reads show it, placed uncertainly', pill: 'partial' },
    unclear: { word: 'unclear', pill: 'warn' },
    reads_disagree: { word: 'reads don’t show it', pill: 'fail' },
  };
  const POP: Record<string, string> = { gnomad_filtered: 'gnomAD flags the site', common: 'common (1% or more)', rare: 'rare (under 1%)', not_seen: 'not seen' };
  const VCF_SAID: Record<string, string> = { NoCall: 'looked, couldn’t decide', RefCall: 'looked, saw no variant', PASS: 'called a different change here' };
  const smallPct = (f: number | null | undefined) => (f == null ? '—' : f >= 0.1 ? pct(f, 0) : f >= 0.01 ? pct(f, 1) : pct(f, 2));
  const seenText = (r: any) => (r.population === 'gnomad_filtered' ? 'gnomAD flags this site as unreliable'
    : [r.gnomad_af != null ? `${smallPct(r.gnomad_af)} in gnomAD` : '', r.af_1kg != null ? `${smallPct(r.af_1kg)} in 1000 Genomes` : ''].filter(Boolean).join(', ') || 'not seen in gnomAD or 1000 Genomes');
  const scShown = $derived(sc.filter((r) => ['reads_support', 'reads_flagged'].includes(r.verdict) && r.population !== 'gnomad_filtered'));
  const scOther = $derived(sc.filter((r) => !scShown.includes(r)));

  const fmtBp = (v: number) => `${(v / 1e6).toFixed(0)} Mb`;
  const subColor = (s: string) => (['A>G', 'G>A', 'C>T', 'T>C'].includes(s) ? 'var(--accent)' : 'var(--accent-2)');
  const chromLabel = (c: string) => (c === 'MT' ? 'M' : c);
</script>

{#await load}
  <p class="muted">Loading QC…</p>
{:then [qc, cov, aln, fq, full, inv, sections, selfc, bench, second]}
  {#if !qc}
    <Card title="No QC results yet"><p>Run <code>pixi run wgs run</code> to compute data-quality checks.</p></Card>
  {:else}
    {@const checks = qc.assessment as Check[]}
    {@const passed = checks.filter((c) => c.status === 'pass').length}
    {@const sex = qc.inferred_sex}
    {@const v = qc.sections.vcf}
    <div class="grid">
      <section class="hero" class:allpass={passed === checks.length}>
        <div class="big">{passed === checks.length ? '✅' : '⚠️'}</div>
        <div>
          <h2>{passed} of {checks.length} quality checks pass</h2>
          <p class="lead">
            {#if passed === checks.length}
              Your sequencing data looks like a clean, complete, single-person human genome. That matters because every later section trusts these calls — if the data were noisy, contaminated or patchy, the interpretations would be too.
            {:else}
              Some checks fall outside the usual range. Results in other sections that depend on the affected data are flagged accordingly.
            {/if}
          </p>
          <Explain learn="Reading the QC report">
            <p>Quality control (QC) asks “can we trust this data?” before asking “what does it mean?”. Each check below compares one number from your data with the range seen in good-quality genomes, and explains why it matters. A check that passes doesn’t say anything about your health — it says the measurement worked.</p>
          </Explain>
        </div>
      </section>

      <Card title="The checks" subtitle="Each bar shows the healthy range in green and your value as the marker.">
        <div class="grid cols-3 checks">
          {#each checks as c}
            {@const meta = META[c.metric] ?? { name: c.metric, term: c.metric, fmt: String, why: '' }}
            <div class="check">
              <div class="row"><span class="pill {c.status}">{c.status === 'pass' ? '✓ pass' : c.status}</span><span class="spacer"></span><Term t={meta.term} label="?" learn="Reading the QC report" /></div>
              <div class="name">{meta.name}</div>
              <div class="kpi">{meta.fmt(c.value)}</div>
              <RangeBar value={c.value} low={c.low} high={c.high} fmt={meta.fmt} />
              <p class="small muted why">{meta.why}</p>
            </div>
          {/each}
        </div>
      </Card>

      {#if cov}
        <Card title="Coverage: how many reads cover each position" subtitle="Measured from every read in your CRAM file.">
          <div class="grid cols-4">
            <div><div class="kpi">{cov.autosomal_mean_depth.toFixed(1)}×</div><div class="muted small">average <Term t="Coverage" label="depth" learn="How whole-genome sequencing works" /> (ordered 30×)</div></div>
            <div><div class="kpi">{pct(cov.fraction_at_least['20'])}</div><div class="muted small">of the genome at ≥ 20 reads</div></div>
            <div><div class="kpi">{pct(cov.callable_fraction)}</div><div class="muted small"><Term t="Callable" label="callable" /> (enough good reads to trust a genotype)</div></div>
            <div><div class="kpi">{int(cov.mt_mean_depth)}×</div><div class="muted small"><Term t="chrM" label="mitochondrial DNA" /> depth — each cell has hundreds of copies</div></div>
          </div>
          <div class="grid cols-2" style="margin-top: 14px">
            <div>
              <h4>Share of genome with at least N reads</h4>
              <Bars data={Object.entries(cov.fraction_at_least).map(([k, f]) => ({ label: `≥${k}×`, value: (f as number) * 100, color: +k <= 20 ? 'var(--good)' : 'var(--accent)' }))} fmt={(x) => `${x.toFixed(0)}%`} height={170} />
            </div>
            <div>
              <h4>Average depth per chromosome</h4>
              <Bars data={[...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'].filter((c) => cov.per_chrom[c]).map((c) => [c, cov.per_chrom[c]] as [string, any]).map(([c, s]: [string, any]) => ({ label: chromLabel(c), value: s.mean, color: c === 'X' || c === 'Y' ? 'var(--accent-2)' : 'var(--accent)' }))} fmt={(x) => `${x.toFixed(0)}×`} height={170} />
              <p class="small faint">X and Y (pink) sit at about half depth because you have one copy of each instead of two.</p>
            </div>
          </div>
          <h4 style="margin-top: 12px">Depth along every chromosome</h4>
          {#if bins.length}<GenomeDepth {bins} mean={cov.autosomal_mean_depth} />{:else}<p class="muted small">Loading genome-wide depth…</p>{/if}
          <Explain learn="How whole-genome sequencing works">
            <p>Sequencing reads random short pieces of your DNA many times over. <strong>Depth</strong> is how many of those pieces overlap a given letter. At 30–40× each position is seen dozens of times, so a single misread letter is easily outvoted. Regions with very low depth can’t be called confidently either way — the dashboard marks results there as “not covered” rather than “normal”.</p>
          </Explain>
        </Card>
      {/if}

      <Card title="Your variant calls" subtitle="From the DeepVariant VCF: positions where your DNA differs from the reference.">
        <div class="grid cols-4">
          <div><div class="kpi">{compact(v.snv_count)}</div><div class="muted small"><Term t="SNV" label="SNVs" /> (single letters)</div></div>
          <div><div class="kpi">{compact(v.indel_count)}</div><div class="muted small"><Term t="Indel" label="indels" /> ({compact(v.insertion_count)} ins, {compact(v.deletion_count)} del)</div></div>
          <div><div class="kpi">{v.titv.toFixed(2)}</div><div class="muted small"><Term t="Ti/Tv" /> ratio</div></div>
          <div><div class="kpi">{v.het_hom_ratio.toFixed(2)}</div><div class="muted small"><Term t="Heterozygous" label="het" /> : <Term t="Homozygous" label="hom" /> ratio</div></div>
        </div>
        <div class="grid cols-2" style="margin-top: 14px">
          <div>
            <h4>Which letter changed to which</h4>
            <Bars data={v.hist.substitutions.map(([s, n]: [string, number]) => ({ label: s, value: n, color: subColor(s) }))} fmt={compact} height={180} rotate />
            <p class="small faint"><span style="color: var(--accent)">Blue</span> = <Term t="Transition" label="transitions" /> (A↔G, C↔T), <span style="color: var(--accent-2)">pink</span> = transversions. Biology makes transitions about twice as common.</p>
          </div>
          <div>
            <h4>Allele balance at heterozygous sites</h4>
            {#if v.hist.het_vaf?.length}
              <Line points={v.hist.het_vaf} xFmt={(x) => x.toFixed(1)} yFmt={compact} xLabel="share of reads showing the variant (VAF)" marker={{ x: 0.5, label: '50%' }} height={180} />
              <p class="small faint">A clean single sample peaks at 50% — one copy from each parent.</p>
            {:else}<p class="muted small">This VCF has no read counts.</p>{/if}
          </div>
          <div>
            <h4>Read depth at variant sites</h4>
            {#if v.hist.depth?.length}
              <Line points={v.hist.depth} xFmt={(x) => `${x.toFixed(0)}×`} yFmt={compact} xLabel="reads covering the variant" height={170} />
            {:else}<p class="muted small">This VCF has no read counts.</p>{/if}
          </div>
          <div>
            <h4>Genotype quality (GQ)</h4>
            <Line points={v.hist.gq} xFmt={(x) => x.toFixed(0)} yFmt={compact} xLabel="GQ (Phred: 20 = 99%, 30 = 99.9% confident)" height={170} />
          </div>
          <div>
            <h4>Indel sizes</h4>
            <Bars data={v.hist.indel_length.filter(([l]: [number]) => Math.abs(l) <= 12 && l !== 0).map(([l, n]: [number, number]) => ({ label: String(l), value: n, color: l < 0 ? 'var(--accent-2)' : 'var(--accent)' }))} fmt={compact} height={170} />
            <p class="small faint">Negative = deletion (pink), positive = insertion (blue), in letters. Single-letter indels dominate, as expected.</p>
          </div>
          <div>
            <h4>Record types in the VCF</h4>
            <Bars data={Object.entries(v.records_by_filter).map(([k, n]) => ({ label: k, value: n as number, color: k === 'PASS' ? 'var(--good)' : 'var(--faint)' }))} fmt={compact} height={170} />
            <p class="small faint"><Term t="PASS" label="PASS" /> = confident variant. RefCall = checked, matches the reference. NoCall = too little evidence to decide.</p>
          </div>
        </div>
      </Card>

      <div class="grid cols-2">
        <Card title="Genetic sex check" subtitle="Inferred from the data, independent of what the provider was told.">
          <div class="row"><span class="kpi">{sex.call ?? '—'}</span><span class="muted">chromosomal sex</span></div>
          {#if !sex.call}<p class="muted small">Not checked: this release has no X or Y chromosome data to compare.</p>{/if}
          <table class="small">
            <tbody>
              <tr><td>X depth relative to autosomes</td><td class="num">{sex.x_ratio?.toFixed(2)}</td><td class="faint">≈0.5 for one X, ≈1.0 for two</td></tr>
              <tr><td>Y depth relative to autosomes</td><td class="num">{sex.y_ratio?.toFixed(2)}</td><td class="faint">≈0 without a Y</td></tr>
              <tr><td>Heterozygous fraction on X (outside <Term t="PAR" />)</td><td class="num">{pct(sex.chrX_nonpar_het_fraction)}</td><td class="faint">near 0 with one X</td></tr>
              <tr><td>Confident SNVs on Y</td><td class="num">{int(sex.chrY_pass_snvs)}</td><td class="faint">thousands with a Y</td></tr>
            </tbody>
          </table>
          <Explain learn="Two copies">
            <p>With one X and one Y, almost the whole X is present in a single copy (<Term t="Hemizygous" label="hemizygous" />), so it shows half the depth and almost no heterozygous sites. This check catches sample swaps and matters later: X-linked conditions are read differently with one X.</p>
          </Explain>
        </Card>

        {#if qc.sections.concordance}
        <Card title="Agreement with the provider’s genotype file" subtitle="Two independently produced files compared position by position.">
          {@const c = qc.sections.concordance}
          <div class="row"><span class="kpi">{c.concordance_pct.toFixed(2)}%</span><span class="muted">of {int(c.compared)} positions agree</span></div>
          <table class="small">
            <tbody>
              <tr><td>Both say “variant” with the same genotype</td><td class="num">{int(c.concordant_variant)}</td></tr>
              <tr><td>Both say “matches the reference”</td><td class="num">{int(c.concordant_reference)}</td></tr>
              <tr><td>Disagree</td><td class="num">{int(c.discordant)}</td></tr>
              <tr><td>Couldn’t be checked</td><td class="num">{int(c.unverifiable)}</td></tr>
            </tbody>
          </table>
          <Explain>
            <p>tellmeGen’s 23andMe-style text file and the VCF were produced by separate steps. Full agreement means the files belong together and were converted without errors, so results from either file can be trusted interchangeably.</p>
          </Explain>
        </Card>
        {/if}
      </div>

      {#if rc.length}
        {@const judged = rc.filter((r) => r.expected != null)}
        {@const ok = judged.filter((r) => r.status === 'confirmed').length}
        <Card title="Re-counted from the raw reads" subtitle="At every position a finding rests on, each read in your CRAM was counted again, separately from the variant caller.">
          <div class="row"><span class="kpi">{int(ok)} of {int(judged.length)}</span><span class="muted">re-counted genotypes match the call</span></div>
          {@const sources = Object.entries(RC_SOURCE).filter(([src]) => rc.some((r) => r.source === src))}
          <div class="tablewrap"><table class="small rc-table">
            <thead><tr><th>Verdict</th>{#each sources as [, label]}<th class="num">{label}</th>{/each}</tr></thead>
            <tbody>
              {#each RC_ORDER as st}
                {#if rc.some((r) => r.status === st)}
                  <tr><td><span class="pill {RECOUNT[st].pill}" title={RECOUNT[st].title}>{RECOUNT[st].word}</span></td>
                    {#each sources as [src]}{@const n = rc.filter((r) => r.source === src && r.status === st).length}<td class="num" class:faint={!n}>{int(n)}</td>{/each}</tr>
                {/if}
              {/each}
            </tbody>
          </table></div>
          {#if rcExceptions.length}
            <h4>Where the re-count adds something</h4>
            <div class="tablewrap"><table class="small rc-table rc-stack">
              <thead><tr><th>Gene / position</th><th>Change</th><th class="num">Reads: usual / variant</th><th>Call says</th><th>Reads say</th><th>Verdict</th></tr></thead>
              <tbody>
                {#each rcExceptions as r (r.key)}
                  <tr>
                    <td><strong>{r.gene ?? '—'}</strong><div class="mono faint">{r.rsid ?? `${r.chrom}:${r.pos}`}</div></td>
                    <td class="mono" data-label="Change">{shortAllele(r.ref)}→{shortAllele(r.alt)}</td>
                    <td class="num" data-label="Reads: usual / variant">{r.ref_reads} / {r.alt_reads}{#if r.low_mapq_reads}<div class="faint">+{r.low_mapq_reads} placed uncertainly</div>{/if}</td>
                    <td data-label="Call says">{copiesWord(r.expected)}</td>
                    <td data-label="Reads say">{copiesWord(r.observed)}{#if r.observed != null && r.clear === false}<div class="faint">too few reads to be sure</div>{/if}</td>
                    <td data-label="Verdict"><span class="pill {RECOUNT[r.status]?.pill}" title={RECOUNT[r.status]?.title}>{RECOUNT[r.status]?.word}</span>
                      {#each r.flags ?? [] as f}<div class="faint">{RC_FLAG[f] ?? f}</div>{/each}
                      {#if r.note}<div class="faint">{r.note}</div>{/if}
                      {#if (r.difficult ?? []).length}<div class="faint">hard to sequence: {difficultText(r.difficult)}</div>{/if}</td>
                  </tr>
                {/each}
              </tbody>
            </table></div>
          {/if}
          <Explain learn="How this project grades its confidence">
            <p>Your variant file was made by tellmeGen’s variant caller (DeepVariant) from these same reads. Counting them again with a simple, open method checks that the caller read them correctly. Single-letter changes are read straight off each read. For insertions and deletions, each whole read is laid against both versions of the DNA there and counted for the one it fits better, which also catches reads the aligner partly set aside (<Term t="Soft clip" label="soft clips" />). A binomial model (<a href="https://pubmed.ncbi.nlm.nih.gov/21903627/" target="_blank" rel="noreferrer">Li 2011</a>) then asks whether the counts fit 0, 1 or 2 copies, and needs that answer to be at least 100 times likelier than the next, with at least 8 usable reads.</p>
            <p>If the reads fit a different genotype, the finding’s call drops to Low. If they are too few or unclear, or many could belong to a look-alike stretch of DNA elsewhere (low <Term t="MAPQ" />), the call is capped at Medium. A re-count never raises a grade. Positions the variant file left uncalled are shown here for information only. What this can’t catch is an error in the reads themselves; the Genome in a Bottle benchmark measures that.</p>
          </Explain>
        </Card>
      {/if}

      {#if bench?.recount}
        {@const grp = (k: string, r: string) => bench.recount.groups.find((x: any) => x.kind === k && x.region === r)}
        {@const easy = grp('all', 'easy')}
        {@const hard = grp('all', 'difficult')}
        {@const d = bench.data}
        {@const cal = bench.caller}
        <Card title="Checked on a genome whose answer is known" subtitle={`Genome in a Bottle’s reference person HG002, whose DNA letters are known almost exactly, went through the same checks: ${int(d.truth_alleles)} true variants and ${int(d.decoys)} variants HG002 doesn’t have, in the same genes and medicine positions as yours.`}>
          <div class="grid cols-4">
            <div><div class="kpi">{pct(easy.right_confirmed)}</div><div class="muted small">of correct calls confirmed by the re-count in ordinary stretches</div></div>
            <div><div class="kpi">{pct(easy.wrong_caught)}</div><div class="muted small">of wrong claims caught in ordinary stretches</div></div>
            <div><div class="kpi">{pct(hard.right_confirmed)}</div><div class="muted small">of correct calls confirmed in <Term t="Hard-to-sequence region" label="hard-to-sequence stretches" /></div></div>
            <div><div class="kpi">{pct(hard.wrong_caught)}</div><div class="muted small">of wrong claims caught there</div></div>
          </div>
          <h4 class="gap">The re-count, by kind of change</h4>
          <div class="tablewrap"><table class="small rc-table">
            <thead><tr><th>Change</th><th>Stretch</th><th class="num">Positions</th><th class="num">Correct calls confirmed</th><th class="num">Correct calls rejected</th><th class="num">Wrong claims caught</th><th class="num">Wrong claims let through</th></tr></thead>
            <tbody>
              {#each ['snv', 'indel'] as k}
                {#each ['easy', 'difficult'] as r}
                  {@const x = grp(k, r)}
                  {#if x}<tr><td>{KIND[k]}</td><td>{REGION[r]}</td><td class="num">{int(x.sites)}</td><td class="num">{pct(x.right_confirmed)}</td><td class="num">{pct(x.right_false_alarm, 2)}</td><td class="num">{pct(x.wrong_caught)}</td><td class="num">{pct(x.wrong_confirmed, 2)}</td></tr>{/if}
                {/each}
              {/each}
            </tbody>
          </table></div>
          <h4 class="gap">Hard-to-sequence stretches, by reason</h4>
          <div class="tablewrap"><table class="small rc-table">
            <thead><tr><th>Why the stretch is hard</th><th class="num">Positions</th><th class="num">Correct calls confirmed</th><th class="num">Wrong claims let through</th></tr></thead>
            <tbody>
              {#each bench.recount.by_reason as x}
                <tr><td>{DIFFICULT[x.reason] ?? x.reason}</td><td class="num">{int(x.sites)}</td><td class="num">{pct(x.right_confirmed)}</td><td class="num">{pct(x.wrong_confirmed, 2)}</td></tr>
              {/each}
            </tbody>
          </table></div>
          {#if cal?.groups}
            {@const cg = (k: string, r: string) => cal.groups.find((x: any) => x.kind === k && x.region === r)}
            {@const all = cg('all', 'all')}
            <h4 class="gap">A second variant caller (bcftools)</h4>
            <p class="small">Of the variants bcftools called, <strong>{pct(all.precision)}</strong> were real (its <Term t="Precision / recall / F1" label="precision" />), and it found <strong>{pct(all.recall)}</strong> of HG002’s true variants (recall). Where the re-count also confirmed a call, {pct(all.precision_if_confirmed)} were real.</p>
            <div class="tablewrap"><table class="small rc-table">
              <thead><tr><th>Change</th><th>Stretch</th><th class="num">Calls that were real</th><th class="num">True variants found</th></tr></thead>
              <tbody>
                {#each ['all', 'snv', 'indel'] as k}
                  {#each ['easy', 'difficult'] as r}
                    {@const x = cg(k, r)}
                    {#if x}<tr><td>{KIND[k]}</td><td>{REGION[r]}</td><td class="num">{pct(x.precision)}</td><td class="num">{pct(x.recall)}</td></tr>{/if}
                  {/each}
                {/each}
              </tbody>
            </table></div>
            {#if cal.by_population?.length}
              <h4 class="gap">Calls by how often other people carry them</h4>
              <div class="tablewrap"><table class="small rc-table">
                <thead><tr><th>In gnomAD and 1000 Genomes</th><th class="num">bcftools calls</th><th class="num">Real</th></tr></thead>
                <tbody>
                  {#each cal.by_population as x}
                    <tr><td>{POP[x.population] ?? x.population}</td><td class="num">{int(x.calls)}</td><td class="num">{int(x.tp)} ({pct(x.precision, 0)})</td></tr>
                  {/each}
                </tbody>
              </table></div>
            {/if}
          {/if}
          <Explain learn="Accuracy benchmarking (Genome in a Bottle)">
            <p>Genome in a Bottle, run by the US National Institute of Standards and Technology, combined many sequencing technologies to work out HG002’s true variants (<a href="https://pubmed.ncbi.nlm.nih.gov/36452119/" target="_blank" rel="noreferrer">Wagner 2022</a>). HG002’s reads from a DNBSEQ sequencer, the same family as yours, were downloaded for these regions only and thinned from {d.full_depth}× to {d.depth}× to match your genome. The re-count was given every true variant with its true number of copies, and also wrong claims: one copy too many or too few, and ClinVar variants HG002 doesn’t carry. bcftools’ calls were matched to the truth with RTG vcfeval, which recognises the same change written in different ways (<a href="https://pubmed.ncbi.nlm.nih.gov/30858580/" target="_blank" rel="noreferrer">Krusche 2019</a>).</p>
            <p>Hard-to-sequence stretches come from Genome in a Bottle’s map of the genome (<a href="https://pubmed.ncbi.nlm.nih.gov/39424793/" target="_blank" rel="noreferrer">Dwarshuis 2024</a>). Errors gather there, mostly insertions and deletions in runs of one repeated letter, whose length short reads miscount. Your findings in such stretches say so in their call reasons. Their grade stays the same, because the re-count still caught {pct(hard.wrong_caught, 0)} of wrong claims there and an unclear re-count already caps a call at Medium.</p>
            <p>gnomAD marks sites whose calls fail its own quality checks (<a href="https://pubmed.ncbi.nlm.nih.gov/32461654/" target="_blank" rel="noreferrer">Karczewski 2020</a>). On HG002 a call at such a site was usually wrong, which is why the second opinion below sorts by it. Not measured here: DeepVariant itself (it doesn’t run on this computer), the X chromosome, and stretches where even Genome in a Bottle’s truth is uncertain.</p>
          </Explain>
        </Card>
      {/if}

      {#if second}
        {@const flagged = bench?.caller?.by_population?.find((x: any) => x.population === 'gnomad_filtered')}
        {@const common = bench?.caller?.by_population?.find((x: any) => x.population === 'common')}
        {@const indelEasy = bench?.caller?.groups?.find((x: any) => x.kind === 'indel' && x.region === 'easy')}
        <Card title="A second opinion on your genome" subtitle="bcftools, a simpler variant caller, read your own reads in the ACMG genes and at the medicine positions; RTG vcfeval then compared its calls with your variant file.">
          <div class="grid cols-3">
            <div><div class="kpi">{int(second.matched)}</div><div class="muted small">variants both callers found</div></div>
            <div><div class="kpi">{int(second.only_bcftools)}</div><div class="muted small">found only by bcftools</div></div>
            <div><div class="kpi">{int(second.only_vcf)}</div><div class="muted small">only in your variant file ({int(second.only_vcf_by_kind?.indel ?? 0)} of them insertions or deletions)</div></div>
          </div>
          {#if second.only_bcftools}
            <h4 class="gap">The {int(second.only_bcftools)} variants only bcftools found</h4>
            <div class="tablewrap"><table class="small rc-table">
              <thead><tr><th>Re-counted reads</th>{#each Object.keys(POP) as p}<th class="num">{POP[p]}</th>{/each}</tr></thead>
              <tbody>
                {#each Object.keys(SC_VERDICT) as v}
                  {#if second.verdicts?.[v]}
                    <tr><td><span class="pill {SC_VERDICT[v].pill}">{SC_VERDICT[v].word}</span></td>
                      {#each Object.keys(POP) as p}{@const n = second.population?.[v]?.[p] ?? 0}<td class="num" class:faint={!n}>{int(n)}</td>{/each}</tr>
                  {/if}
                {/each}
              </tbody>
            </table></div>
            <p class="small">{second.clinvar_pathogenic ? `${second.clinvar_pathogenic} of them are listed in ClinVar as pathogenic or likely pathogenic.` : 'None of them is listed in ClinVar as pathogenic or likely pathogenic.'}
              {#if flagged && common}On HG002, bcftools calls at sites gnomAD flags were real {flagged.tp} times out of {flagged.calls}; at common sites, {pct(common.precision)} were.{/if}</p>
          {/if}
          {#if scShown.length}
            <h4 class="gap">Shown by your reads and not flagged by gnomAD</h4>
            <div class="tablewrap"><table class="small rc-table rc-stack">
              <thead><tr><th>Gene / position</th><th>Change</th><th class="num">Reads: usual / variant</th><th>Your variant file</th><th>Seen in other people</th><th>Notes</th></tr></thead>
              <tbody>
                {#each scShown as r (`${r.chrom}:${r.pos}:${r.alt}`)}
                  <tr>
                    <td><strong>{r.gene ?? '—'}</strong><div class="mono faint">{r.chrom}:{r.pos}</div></td>
                    <td class="mono" data-label="Change">{shortAllele(r.ref)}→{shortAllele(r.alt)}</td>
                    <td class="num" data-label="Reads: usual / variant">{r.ref_reads} / {r.alt_reads}{#if r.low_mapq_reads}<div class="faint">+{r.low_mapq_reads} placed uncertainly</div>{/if}</td>
                    <td data-label="Your variant file">{VCF_SAID[r.vcf_filter] ?? 'no record'}{#if r.vcf_pass_nearby != null}<div class="faint">a different change {r.vcf_pass_nearby} letter{r.vcf_pass_nearby === 1 ? '' : 's'} away</div>{/if}</td>
                    <td data-label="Seen in other people">{seenText(r)}</td>
                    <td class="notes" data-label="Notes">{#if r.clinvar_significance}<div>ClinVar: {r.clinvar_significance.replaceAll('_', ' ')}</div>{/if}
                      {#if (r.difficult ?? []).length}<div class="faint">hard to sequence: {difficultText(r.difficult)}</div>{/if}
                      {#if r.verdict === 'reads_flagged'}<div class="faint">{RC_FLAG.look_alike}</div>{/if}</td>
                  </tr>
                {/each}
              </tbody>
            </table></div>
          {/if}
          {#if scOther.length}
            <details class="small"><summary>The other {scOther.length}: unclear, or at sites gnomAD flags</summary>
              <div class="tablewrap"><table class="small rc-table">
                <thead><tr><th>Gene / position</th><th>Change</th><th class="num">Reads</th><th>Re-count</th><th>Seen in other people</th></tr></thead>
                <tbody>
                  {#each scOther as r (`${r.chrom}:${r.pos}:${r.alt}`)}
                    <tr><td><strong>{r.gene ?? '—'}</strong> <span class="mono faint">{r.chrom}:{r.pos}</span></td><td class="mono">{shortAllele(r.ref)}→{shortAllele(r.alt)}</td>
                      <td class="num">{r.ref_reads} / {r.alt_reads}</td><td>{SC_VERDICT[r.verdict]?.word ?? r.verdict}</td><td>{seenText(r)}</td></tr>
                  {/each}
                </tbody>
              </table></div>
            </details>
          {/if}
          <Explain learn="Accuracy benchmarking (Genome in a Bottle)">
            <p>Two <Term t="Variant caller" label="variant callers" /> reading the same reads mostly disagree where the reads are hard to interpret: runs of a repeated letter, and stretches with a look-alike copy elsewhere. A variant only bcftools found is a lead, not a finding. It has no annotation or grade, and nothing in the other tabs changes because of it. “Looked, couldn’t decide” means DeepVariant saw the evidence but gave no genotype (NoCall); only a clinical test could settle such a position if it ever mattered.</p>
            <p>Several differences close together that gnomAD flags, with almost nobody carrying two copies, are the usual sign of reads from a look-alike copy placed in the wrong spot. The re-count can’t tell, because it counts the same placed reads. Most variants only in your file are insertions or deletions, which bcftools handles less well{#if indelEasy} (on HG002 it found {pct(indelEasy.recall, 0)} of them even in ordinary stretches){/if}.</p>
          </Explain>
        </Card>
      {/if}

      <div class="grid cols-2">
        {#if fq}
          <Card title="Raw read quality (FASTQ)" subtitle={`Sampled ${int(fq.sampled_reads_per_file)} reads from each file.`}>
            {#each Object.entries(fq.files as Record<string, any>) as [mate, f]}
              <h4>{mate === 'R1' ? 'Read 1' : 'Read 2'} <span class="faint small">— <Term t="Phred" label={`${f.q30_pct}% of letters ≥ Q30`} learn="How whole-genome sequencing works" />, GC {f.gc_pct}%</span></h4>
              <Line points={f.mean_quality_by_cycle.map((q: number, i: number) => [i + 1, q])} yMin={20} yMax={42} xFmt={(x) => x.toFixed(0)} yFmt={(y) => `Q${y.toFixed(0)}`} xLabel="position in read (cycle)" height={140} color={mate === 'R1' ? 'var(--accent)' : 'var(--accent-2)'} />
            {/each}
            <p class="small faint">Instrument: {fq.files.R1?.platform_guess} (flowcell {fq.files.R1?.flowcell}). Q30 means a 1-in-1,000 chance a letter is wrong; Q35 ≈ 1 in 3,000.</p>
          </Card>
        {/if}
        {#if aln}
          <Card title="How reads aligned" subtitle={`samtools stats on ${aln.region} (${compact(aln.reads)} reads) as a representative sample.`}>
            <table class="small">
              <tbody>
                <tr><td>Mapped to the reference</td><td class="num">{aln.mapped_pct.toFixed(2)}%</td></tr>
                <tr><td>Both ends mapped as a proper pair</td><td class="num">{aln.properly_paired_pct.toFixed(1)}%</td></tr>
                <tr><td><Term t="Duplicate reads" label="Duplicates" /></td><td class="num">{aln.duplicate_pct.toFixed(2)}%</td></tr>
                <tr><td>Ambiguous placement (<Term t="MAPQ" /> 0)</td><td class="num">{aln.mapq0_pct.toFixed(2)}%</td></tr>
                <tr><td>Fragment (insert) size</td><td class="num">{aln.insert_size_mean.toFixed(0)} ± {aln.insert_size_sd.toFixed(0)} bp</td></tr>
              </tbody>
            </table>
            <h4 style="margin-top: 10px">Fragment sizes</h4>
            <Line points={aln.insert_size_hist.filter(([s]: [number]) => s >= 50 && s <= 700)} xFmt={(x) => `${x.toFixed(0)}`} yFmt={compact} xLabel="DNA fragment length (bp) — both ends were read" height={150} />
          </Card>
        {/if}
      </div>

      {#if calib.length}
        <Card title="Genotype calibration" subtitle="Your genotypes at the ~millions of positions used for polygenic scores, against how common each variant is in 1000 Genomes.">
          <div class="tablewrap"><table class="small">
            <thead><tr><th>type</th><th>variant in panel</th><th class="num">sites</th><th class="num">you: none</th><th class="num">one copy</th><th class="num">two copies</th><th class="num">not callable</th><th class="num">your frequency</th><th class="num">panel</th></tr></thead>
            <tbody>
              {#each calib as r}
                <tr class:alarm={r.band === 6 && r.hom_ref_pct > 2}>
                  <td>{r.kind}</td><td>{r.panel_af_band}</td><td class="num">{compact(r.sites)}</td>
                  <td class="num">{r.hom_ref_pct.toFixed(1)}%</td><td class="num">{r.het_pct.toFixed(1)}%</td><td class="num">{r.hom_alt_pct.toFixed(1)}%</td>
                  <td class="num faint">{r.missing_pct.toFixed(1)}%</td><td class="num">{r.your_af != null ? pct(r.your_af) : '—'}</td><td class="num faint">{pct(r.panel_af)}</td>
                </tr>
              {/each}
            </tbody>
          </table></div>
          <Explain>
            <p>A single person’s genotypes should roughly follow population frequencies: where more than 99% of people carry a variant, you should almost never be called as not carrying it. A high “you: none” share in that last band would mean genotypes are being read wrongly — this once caught an indel bug, and the self-check <em>Genotype calibration</em> now guards it on every run.</p>
            <p>“Not callable” positions (common for insertions/deletions described differently by different tools) are left out of your score and filled with the population average, as the PGS Catalog’s own pipeline does.</p>
          </Explain>
        </Card>
      {/if}

      <Card title="Reference genome proof" subtitle={qc.sections.reference?.reads_decoded == null ? 'Every reference letter the VCF names must match the reference genome used here.' : 'The CRAM can only be decoded with the exact reference it was made with.'}>
        {@const r = qc.sections.reference}
        {#if r?.reads_decoded == null}
          <p>The reference was <strong>{r?.origin ?? 'located'}</strong>.</p>
        {:else}
          <p>The reference was <strong>{r?.origin ?? 'located'}</strong>, then checked against fingerprints stored inside your CRAM file.</p>
        {/if}
        <table class="small">
          <tbody>
            {#if r?.reads_decoded != null}<tr><td>Quick check: reads decoded across all {r?.contigs} <Term t="Contig" label="contigs" /></td><td class="num">{int(r?.reads_decoded)}</td></tr>{/if}
            <tr><td>VCF reference letters that disagree with the FASTA</td><td class="num">{int(r?.vcf_ref_mismatches)}</td></tr>
            {#if full}
              <tr><td><strong>Full proof: every read decoded</strong></td><td class="num">{int(full.reads_decoded)}</td></tr>
              <tr><td>Fingerprint (<Term t="MD5" />) mismatches</td><td class="num">{int(full.md5_mismatches)}</td></tr>
            {/if}
          </tbody>
        </table>
        <Explain learn="The reference genome">
          <p>A CRAM file saves space by storing only how each read differs from the reference. Each block also stores an MD5 fingerprint of the reference stretch it used, so decoding with even one wrong letter fails loudly. Zero mismatches across every read proves the rebuilt reference is the one tellmeGen used.</p>
        </Explain>
      </Card>

      {#if selfc}
        {@const ok = selfc.counts.fail + selfc.counts.error === 0}
        <Card title={`Self-checks: ${selfc.counts.pass} of ${selfc.checks.length} passed`} subtitle="Automatic consistency checks run on this release’s own results before it is published — they catch bugs, not biology.">
          <ul class="selfc small">
            {#each selfc.checks as c}
              <li class={c.status}>
                <span class="pill {c.status === 'pass' ? 'pass' : c.status === 'skipped' ? 'soon' : 'fail'}">{c.status === 'pass' ? '✓' : c.status === 'skipped' ? '–' : '✗'} {c.status}</span>
                <span><strong>{c.title}</strong> <span class="muted">— {c.why}</span>
                  {#if c.status === 'fail'}<br /><span class="warn">{c.count} problem{c.count > 1 ? 's' : ''}: {c.examples.join(', ')}{c.count > c.examples.length ? '…' : ''}</span>{/if}
                  {#if c.detail}<br /><span class="faint">{c.detail}</span>{/if}</span>
              </li>
            {/each}
          </ul>
          {#if !ok}<p class="small warn">Some checks failed: results they mention may be shown wrongly. Re-run the pipeline after updating, or report it.</p>{/if}
        </Card>
      {/if}

      {#await audit then au}
        <Card title="Source re-verification" subtitle={au ? `This release re-checked against live ClinVar, GWAS Catalog, PGS Catalog and PubMed on ${String(au.audited_at).slice(0, 10)}.` : 'Every finding names its exact source record; wgs audit re-fetches those records and confirms they still say what this release says.'}>
          {#if au}
            <ul class="selfc small">
              {#each Object.entries(au.counts as Record<string, Record<string, number>>) as [k, c]}
                {@const bad = (c.mismatch ?? 0) + (c.missing ?? 0) + (c.broken ?? 0) + (c.error ?? 0)}
                {@const tot = Object.values(c).reduce((a: number, b: number) => a + b, 0)}
                <li class={bad ? 'fail' : 'pass'}><span class="pill {bad ? 'warn' : 'pass'}" style="white-space: nowrap">{bad ? '!' : '✓'} {tot - bad}/{tot}</span>
                  <span><strong>{SRC[k] ?? k}</strong> <span class="muted">— {Object.entries(c).map(([s, n]) => `${n} ${s}`).join(', ')}</span></span></li>
              {/each}
            </ul>
            {#if au.problems?.length}
              <details class="small"><summary>{au.problems.length} item{au.problems.length > 1 ? 's' : ''} differ from the live source</summary>
                <ul>{#each au.problems.slice(0, 50) as p}<li><span class="mono">{p.record}</span> — {p.status}{p.live ? `: live says ${p.live}` : ''}{p.stored ? `; this release says ${p.stored}` : ''}{p.detail ? ` (${p.detail})` : ''}{#if p.url} · <a href={p.url} target="_blank" rel="noreferrer">source</a>{/if}</li>{/each}</ul>
                <p class="faint">“Changed upstream” means the source was updated after this release’s snapshot; refresh knowledge and re-run to pick it up.</p>
              </details>
            {/if}
          {:else}
            <p class="small muted">Not audited yet. Run <code>pixi run wgs audit</code> (needs internet).</p>
          {/if}
        </Card>
      {/await}

      {#if inv}
        <Card title="Input files found" subtitle="The pipeline finds files by content, not by name or folder; sections adapt to whatever is present.">
          <div class="tablewrap"><table class="small files">
            <thead><tr><th>kind</th><th>file</th><th class="num">size</th></tr></thead>
            <tbody>
              {#each inv.files.filter((f: any) => f.kind !== 'report') as f}
                <tr><td><span class="pill">{f.kind}</span></td><td class="mono">{f.folder}/{f.name}</td><td class="num">{compact(f.size)}B</td></tr>
              {/each}
              <tr><td><span class="pill">report</span></td><td class="muted">{inv.files.filter((f: any) => f.kind === 'report').length} provider PDF/zip reports</td><td></td></tr>
            </tbody>
          </table></div>
        </Card>
      {/if}
    </div>
  {/if}
{:catch e}
  <Card title="Could not load QC"><p class="muted">{String(e)}</p></Card>
{/await}

<style>
  .selfc { list-style: none; padding: 0; margin: 0; display: grid; gap: 6px; }
  .selfc li { display: grid; grid-template-columns: 80px 1fr; gap: 10px; align-items: start; }
  .hero { display: flex; gap: 18px; align-items: flex-start; padding: 22px; border-radius: 16px; border: 1px solid var(--border);
    background: linear-gradient(135deg, color-mix(in srgb, var(--warn) 12%, var(--panel)), var(--panel)); }
  .hero.allpass { background: linear-gradient(135deg, color-mix(in srgb, var(--good) 14%, var(--panel)), var(--panel)); }
  .hero h2 { font-size: 1.5rem; }
  .big { font-size: 2.6rem; line-height: 1; }
  .check { border: 1px solid var(--border); border-radius: 12px; padding: 12px 14px; background: var(--panel-2); }
  .name { font-weight: 600; margin-top: 6px; }
  .why { margin: 10px 0 0; }
  h4 { font-size: 0.9rem; margin: 6px 0; color: var(--muted); font-weight: 600; }
  h4.gap { margin-top: 18px; }
  details > summary { cursor: pointer; margin-top: 12px; }
  td.num { white-space: nowrap; }
  tr.alarm td { color: var(--warn); font-weight: 600; }
  .tablewrap { overflow-x: auto; }
  .files td.mono { overflow-wrap: anywhere; }
  .rc-table th.num { white-space: normal; }
  .rc-table .pill { white-space: nowrap; }
  .rc-table td:last-child { white-space: normal; min-width: 14ch; }
  @media (max-width: 700px) {
    .rc-stack thead { display: none; }
    .rc-stack, .rc-stack tbody, .rc-stack tr, .rc-stack td { display: block; }
    .rc-stack td { border: 0; padding: 2px 6px; white-space: normal; text-align: left; }
    .rc-stack td[data-label]::before { content: attr(data-label) ': '; color: var(--faint); }
    .rc-stack td .faint, .rc-stack td .pill { display: inline; margin-left: 4px; }
    .rc-stack td.notes .faint { display: block; margin-left: 0; }
    .rc-stack tr { border-bottom: 1px solid var(--border); padding: 6px 0; }
  }
</style>
