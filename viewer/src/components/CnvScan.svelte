<script lang="ts">
  // Deletion/duplication scan of the ACMG genes (ADR-021): two independent signals (breakpoint reads from Delly
  // and GC-corrected exon depth) must agree before anything is reported; what couldn't be judged is listed.
  import type { Manifest } from '../lib/data';
  import { query, useRelease, type Row } from '../lib/db';
  import { OVERALL_PILL } from '../lib/evidence';

  let { m }: { m: Manifest } = $props();
  let genes = $state<Row[]>([]);
  let events = $state<Row[]>([]);
  let loaded = $state(false);
  $effect(() => {
    (async () => {
      await useRelease(m);
      if (m.tables.cnv_genes) genes = (await query('SELECT * FROM cnv_genes ORDER BY gene')).rows;
      if (m.tables.cnv_events) events = (await query(`SELECT * FROM cnv_events ORDER BY CASE status WHEN 'possible_exon_deletion' THEN 0
        WHEN 'possible_exon_duplication' THEN 1 WHEN 'long_depth_signal' THEN 2 WHEN 'unconfirmed_signal' THEN 3 ELSE 4 END`)).rows;
      loaded = true;
    })();
  });
  const STATUS: Record<string, [string, string, string]> = {
    possible_exon_deletion: ['Possible deletion', 'fail', 'Both methods agree; needs a clinical test (MLPA or microarray)'],
    possible_exon_duplication: ['Possible duplication', 'partial', 'Both methods agree; meaning often uncertain'],
    long_depth_signal: ['Depth drop over many exons', 'warn', 'Seen in depth only, over enough exons to deserve a clinical look'],
    unconfirmed_signal: ['Depth signal only', 'soon', 'Seen by one method only; usually noise'],
    common: ['Common', 'pass', 'Common in gnomAD-SV: not a cause of rare disease'],
    no_exon: ['Between exons', 'pass', 'Removes no exon of the gene'],
  };
  const exons = $derived(genes.reduce((a, g) => a + (g.exons ?? 0), 0));
  const skipped = $derived(genes.reduce((a, g) => a + (g.extreme_gc_exons ?? 0) + (g.low_mappability_exons ?? 0), 0));
  const found = $derived(events.filter((e) => String(e.status).startsWith('possible') || e.status === 'long_depth_signal'));
  const list = (x: any) => (Array.isArray(x) ? x : x ? Array.from(x) : []);
  const gn = (e: Row) => { try { return e.gnomad ? JSON.parse(e.gnomad) : null; } catch { return null; } };
</script>

{#if loaded && genes.length}
  <div class="summary">
    <div class="verdict {found.length ? 'warn' : 'good'}">{found.length ? '!' : '✓'}</div>
    <div>
      <p><strong>Deletions and duplications:</strong> {genes.length} genes, {exons.toLocaleString()} exons scanned. {#if found.length}<strong>{found.length}</strong> need a closer look (below).{:else}No deletion or duplication of any exon was found by both methods.{/if}</p>
      <p class="small faint">{skipped} exons couldn’t be judged by read depth (very GC-rich or GC-poor, or too similar to another part of the genome), so a change there would be missed. Two methods must agree: reads that span the change’s edges (Delly) and a drop in read depth after correcting for GC content. Single-genome depth is noisy, so depth alone is never reported as a finding.</p>
    </div>
  </div>
  {#if events.length}
    <details class="small"><summary>{events.length} candidate{events.length === 1 ? '' : 's'} looked at and set aside or flagged</summary>
      <table class="t"><thead><tr><th>Where</th><th>Type</th><th>Result</th><th>Evidence</th></tr></thead><tbody>
        {#each events as e (e.event_id)}
          {@const st = STATUS[e.status] ?? [e.status, '', '']}
          {@const g = gn(e)}
          <tr><td><strong>{list(e.genes).join(', ')}</strong><div class="faint mono">{e.chrom}:{Number(e.start).toLocaleString()}–{Number(e.end).toLocaleString()} ({(e.end - e.start + 1).toLocaleString()} bp)</div></td>
            <td>{e.svtype === 'DEL' ? 'Deletion' : 'Duplication'}{list(e.exons).length ? ` · ${list(e.exons).length} exon${list(e.exons).length > 1 ? 's' : ''}` : ''}</td>
            <td><span class="pill {st[1]}">{st[0]}</span><div class="faint">{st[2]}</div></td>
            <td class="faint">{e.delly ? `Delly: ${e.pe} spanning pairs, ${e.sr} split reads` : 'No breakpoint reads'}{e.mean_ratio != null ? ` · depth ${Number(e.mean_ratio).toFixed(2)}×` : ''}{g ? ` · gnomAD-SV ${(100 * (g.popmax_af ?? g.af)).toFixed(1)}%` : ''} <span class="pill {OVERALL_PILL[e.overall]}">{e.overall}</span></td></tr>
        {/each}
      </tbody></table>
    </details>
  {/if}
{/if}

<style>
  .summary { display: flex; gap: 14px; align-items: center; margin: 12px 0 6px; }
  .summary p { margin: 2px 0; }
  .verdict { flex: none; width: 40px; height: 40px; border-radius: 50%; display: grid; place-items: center; font-size: 1.3rem; font-weight: 700; }
  .verdict.good { color: var(--good); background: color-mix(in srgb, var(--good) 16%, transparent); }
  .verdict.warn { color: var(--warn); background: color-mix(in srgb, var(--warn) 16%, transparent); }
  .t { width: 100%; border-collapse: collapse; margin-top: 6px; }
  .t th { text-align: left; color: var(--faint); font-weight: 600; padding: 4px 6px; border-bottom: 1px solid var(--border); }
  .t td { padding: 6px; border-bottom: 1px solid var(--border); vertical-align: top; }
</style>
