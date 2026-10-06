<script lang="ts">
  // One strip per chromosome showing read depth along it, scaled to chromosome length.
  import { int } from '../lib/format';
  export interface Bin { chrom: string; start: number; depth: number }
  let { bins, mean }: { bins: Bin[]; mean: number } = $props();
  const byChrom = $derived.by(() => {
    const m = new Map<string, Bin[]>();
    for (const b of bins) (m.get(b.chrom) ?? m.set(b.chrom, []).get(b.chrom)!).push(b);
    return [...m.entries()];
  });
  const longest = $derived(Math.max(...bins.map((b) => b.start)) + 100_000);
  const cap = $derived(mean * 2);
  const W = 1000, H = 22;
  function path(list: Bin[]) {
    return list.map((b, i) => `${i ? 'L' : 'M'}${((b.start / longest) * W).toFixed(1)},${(H - (Math.min(b.depth, cap) / cap) * H).toFixed(1)}`).join('');
  }
</script>

<div class="gd">
  {#each byChrom as [chrom, list]}
    {@const len = list[list.length - 1].start + 100_000}
    {@const avg = list.reduce((s, b) => s + b.depth, 0) / list.length}
    <div class="row">
      <span class="lab mono">{chrom}</span>
      <svg viewBox="0 0 {W} {H}" preserveAspectRatio="none" class="strip">
        <title>chr{chrom}: {int(len)} bp, average depth {avg.toFixed(1)}×</title>
        <rect x="0" y="0" width={(len / longest) * W} height={H} class="bg" />
        <line x1="0" x2={(len / longest) * W} y1={H / 2} y2={H / 2} class="mean" />
        <path d={path(list)} class="ln" />
      </svg>
      <span class="avg mono">{avg.toFixed(0)}×</span>
    </div>
  {/each}
  <p class="faint small legend">Each strip is one chromosome drawn to scale; the dashed line is your genome-wide mean ({mean.toFixed(1)}×) and the top of a strip is 2× the mean. Flat gaps at zero are regions nobody can sequence well (centromeres, the short arms of 13, 14, 15, 21, 22, and repetitive blocks), not problems with your sample.</p>
</div>

<style>
  .row { display: grid; grid-template-columns: 26px 1fr 40px; align-items: center; gap: 8px; height: 24px; }
  .lab { color: var(--muted); font-size: 0.75rem; text-align: right; }
  .avg { color: var(--faint); font-size: 0.72rem; }
  .strip { width: 100%; height: 22px; display: block; }
  .bg { fill: var(--panel-2); }
  .mean { stroke: var(--faint); stroke-dasharray: 3 4; stroke-width: 1; vector-effect: non-scaling-stroke; }
  .ln { fill: none; stroke: var(--accent); stroke-width: 1; vector-effect: non-scaling-stroke; }
  .legend { margin-top: 8px; }
</style>
