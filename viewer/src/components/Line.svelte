<script lang="ts">
  // Simple line/area chart. `band` shades an x-range (e.g. the expected window).
  let { points, height = 180, xFmt = (v: number) => String(v), yFmt = (v: number) => String(v), xLabel = '', yLabel = '', yMin = undefined, yMax = undefined, color = 'var(--accent)', marker = undefined }:
    { points: [number, number][]; height?: number; xFmt?: (v: number) => string; yFmt?: (v: number) => string; xLabel?: string; yLabel?: string; yMin?: number; yMax?: number; color?: string; marker?: { x: number; label: string } } = $props();
  const W = 600;
  const pad = { l: 50, r: 10, t: 10, b: 34 };
  const xs = $derived(points.map((p) => p[0]));
  const ys = $derived(points.map((p) => p[1]));
  const x0 = $derived(Math.min(...xs)), x1 = $derived(Math.max(...xs));
  const y0 = $derived(yMin ?? Math.min(0, ...ys)), y1 = $derived(yMax ?? (Math.max(...ys) || 1));
  const sx = (v: number) => pad.l + ((v - x0) / (x1 - x0 || 1)) * (W - pad.l - pad.r);
  const sy = (v: number) => pad.t + (1 - (Math.min(Math.max(v, y0), y1) - y0) / (y1 - y0 || 1)) * (height - pad.t - pad.b);
  const path = $derived(points.map((p, i) => `${i ? 'L' : 'M'}${sx(p[0]).toFixed(1)},${sy(p[1]).toFixed(1)}`).join(''));
  const area = $derived(points.length ? `${path}L${sx(x1)},${sy(y0)}L${sx(x0)},${sy(y0)}Z` : '');
  const yt = $derived([0, 0.25, 0.5, 0.75, 1].map((f) => y0 + f * (y1 - y0)));
  const xt = $derived([0, 0.25, 0.5, 0.75, 1].map((f) => x0 + f * (x1 - x0)));
</script>

<svg viewBox="0 0 {W} {height}" class="chart" role="img" aria-label={yLabel}>
  {#each yt as t}
    <line x1={pad.l} x2={W - pad.r} y1={sy(t)} y2={sy(t)} class="grid" />
    <text x={pad.l - 6} y={sy(t) + 4} text-anchor="end" class="tick">{yFmt(t)}</text>
  {/each}
  {#each xt as t}<text x={sx(t)} y={height - pad.b + 15} text-anchor="middle" class="tick">{xFmt(t)}</text>{/each}
  <path d={area} fill={color} opacity="0.15" />
  <path d={path} fill="none" stroke={color} stroke-width="1.8" />
  {#if marker}
    <line x1={sx(marker.x)} x2={sx(marker.x)} y1={pad.t} y2={height - pad.b} stroke="var(--accent-2)" stroke-dasharray="4 3" />
    <text x={sx(marker.x) + 5} y={pad.t + 12} class="tick mk">{marker.label}</text>
  {/if}
  {#if xLabel}<text x={(W + pad.l) / 2} y={height - 3} text-anchor="middle" class="tick">{xLabel}</text>{/if}
  {#if yLabel}<text x="11" y={(height - pad.b) / 2} transform="rotate(-90 11 {(height - pad.b) / 2})" text-anchor="middle" class="tick">{yLabel}</text>{/if}
</svg>

<style>
  .chart { width: 100%; height: auto; display: block; }
  .grid { stroke: var(--border); stroke-width: 1; }
  .tick { fill: var(--faint); font-size: 11px; }
  .mk { fill: var(--accent-2); }
</style>
