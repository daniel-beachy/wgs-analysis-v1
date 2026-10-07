<script lang="ts">
  import type { Manifest } from '../lib/data';
  import { loadDoc } from '../lib/data';
  import { query, useRelease } from '../lib/db';
  import { int, pct, compact } from '../lib/format';
  import Card from '../components/Card.svelte';
  import Term from '../components/Term.svelte';
  import Explain from '../components/Explain.svelte';
  import RangeBar from '../components/RangeBar.svelte';
  import Bars from '../components/Bars.svelte';
  import Line from '../components/Line.svelte';
  import GenomeDepth from '../components/GenomeDepth.svelte';

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

  const load = $derived(Promise.all([
    loadDoc(m, 'qc'), loadDoc(m, 'coverage'), loadDoc(m, 'alignment_stats'), loadDoc(m, 'fastq_stats'),
    loadDoc(m, 'reference_full'), loadDoc(m, 'inventory'), loadDoc(m, 'sections'), loadDoc(m, 'checks'),
  ]));

  let bins = $state<{ chrom: string; start: number; depth: number }[]>([]);
  $effect(() => {
    if (!m.tables.coverage_bins) return;
    useRelease(m).then(() => query('SELECT chrom, start, depth FROM coverage_bins WHERE chrom_order < 22 ORDER BY chrom_order, start'))
      .then((r) => (bins = r.rows as any));
  });

  const fmtBp = (v: number) => `${(v / 1e6).toFixed(0)} Mb`;
  const subColor = (s: string) => (['A>G', 'G>A', 'C>T', 'T>C'].includes(s) ? 'var(--accent)' : 'var(--accent-2)');
  const chromLabel = (c: string) => (c === 'MT' ? 'M' : c);
</script>

{#await load}
  <p class="muted">Loading QC…</p>
{:then [qc, cov, aln, fq, full, inv, sections, selfc]}
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
            <Line points={v.hist.het_vaf} xFmt={(x) => x.toFixed(1)} yFmt={compact} xLabel="share of reads showing the variant (VAF)" marker={{ x: 0.5, label: '50%' }} height={180} />
            <p class="small faint">A clean single sample peaks at 50% — one copy from each parent.</p>
          </div>
          <div>
            <h4>Read depth at variant sites</h4>
            <Line points={v.hist.depth} xFmt={(x) => `${x.toFixed(0)}×`} yFmt={compact} xLabel="reads covering the variant" height={170} />
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
          <div class="row"><span class="kpi">{sex.call}</span><span class="muted">chromosomal sex</span></div>
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
      </div>

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

      <Card title="Reference genome proof" subtitle="The CRAM can only be decoded with the exact reference it was made with.">
        {@const r = qc.sections.reference}
        <p>The reference was <strong>{r?.origin ?? 'located'}</strong>, then checked against fingerprints stored inside your CRAM file.</p>
        <table class="small">
          <tbody>
            <tr><td>Quick check: reads decoded across all {r?.contigs} <Term t="Contig" label="contigs" /></td><td class="num">{int(r?.reads_decoded)}</td></tr>
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

      {#if inv}
        <Card title="Input files found" subtitle="The pipeline finds files by content, not by name or folder; sections adapt to whatever is present.">
          <table class="small">
            <thead><tr><th>kind</th><th>file</th><th class="num">size</th></tr></thead>
            <tbody>
              {#each inv.files.filter((f: any) => f.kind !== 'report') as f}
                <tr><td><span class="pill">{f.kind}</span></td><td class="mono">{f.folder}/{f.name}</td><td class="num">{compact(f.size)}B</td></tr>
              {/each}
              <tr><td><span class="pill">report</span></td><td class="muted">{inv.files.filter((f: any) => f.kind === 'report').length} provider PDF/zip reports</td><td></td></tr>
            </tbody>
          </table>
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
  td.num { white-space: nowrap; }
</style>
