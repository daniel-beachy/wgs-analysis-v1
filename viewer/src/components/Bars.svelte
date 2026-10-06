<script lang="ts">
  export interface Bar { label: string; value: number; color?: string; note?: string }
  let { data, height = 180, fmt = (v: number) => String(v), yLabel = '', rotate = false }: { data: Bar[]; height?: number; fmt?: (v: number) => string; yLabel?: string; rotate?: boolean } = $props();
  const W = 600;
  const pad = $derived({ l: 46, r: 8, t: 10, b: rotate ? 46 : 26 });
  const max = $derived(Math.max(...data.map((d) => d.value), 0) || 1);
  const bw = $derived((W - pad.l - pad.r) / Math.max(data.length, 1));
  const y = (v: number) => pad.t + (height - pad.t - pad.b) * (1 - v / max);
  const ticks = $derived([0, 0.5, 1].map((f) => f * max));
</script>

<svg viewBox="0 0 {W} {height}" class="chart" role="img" aria-label={yLabel}>
  {#each ticks as t}
    <line x1={pad.l} x2={W - pad.r} y1={y(t)} y2={y(t)} class="grid" />
    <text x={pad.l - 6} y={y(t) + 4} text-anchor="end" class="tick">{fmt(t)}</text>
  {/each}
  {#each data as d, i}
    <g>
      <title>{d.label}: {fmt(d.value)}{d.note ? ` — ${d.note}` : ''}</title>
      <rect x={pad.l + i * bw + bw * 0.12} y={y(d.value)} width={bw * 0.76} height={Math.max(0, height - pad.b - y(d.value))} rx="2" fill={d.color ?? 'var(--accent)'} />
      {#if rotate}
        <text transform="translate({pad.l + i * bw + bw / 2},{height - pad.b + 10}) rotate(-45)" text-anchor="end" class="tick">{d.label}</text>
      {:else}
        <text x={pad.l + i * bw + bw / 2} y={height - pad.b + 16} text-anchor="middle" class="tick">{d.label}</text>
      {/if}
    </g>
  {/each}
  {#if yLabel}<text x="10" y={pad.t + (height - pad.t - pad.b) / 2} transform="rotate(-90 10 {pad.t + (height - pad.t - pad.b) / 2})" text-anchor="middle" class="tick">{yLabel}</text>{/if}
</svg>

<style>
  .chart { width: 100%; height: auto; display: block; }
  .grid { stroke: var(--border); stroke-width: 1; }
  .tick { fill: var(--faint); font-size: 11px; }
</style>
