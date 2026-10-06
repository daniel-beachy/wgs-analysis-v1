<script lang="ts">
  // Shows where a value sits relative to the expected range (green band).
  let { value, low, high, fmt = (v: number) => String(v) }: { value: number; low: number; high: number; fmt?: (v: number) => string } = $props();
  const span = $derived(Math.max(high - low, Math.abs(value) * 0.1, 1e-9));
  const min = $derived(Math.min(low - span * 0.35, value - span * 0.1));
  const max = $derived(Math.max(high + span * 0.35, value + span * 0.1));
  const x = (v: number) => ((v - min) / (max - min)) * 100;
  const inside = $derived(value >= low && value <= high);
</script>

<div class="rb" title={`expected ${fmt(low)} – ${fmt(high)}`}>
  <div class="track"></div>
  <div class="band" style:left="{x(low)}%" style:width="{x(high) - x(low)}%"></div>
  <div class="mark" class:out={!inside} style:left="{x(value)}%"></div>
  <div class="labels"><span style:left="{x(low)}%">{fmt(low)}</span><span style:left="{x(high)}%">{fmt(high)}</span></div>
</div>

<style>
  .rb { position: relative; height: 30px; margin: 6px 4px 0; }
  .track { position: absolute; top: 6px; left: 0; right: 0; height: 6px; border-radius: 3px; background: var(--chip); }
  .band { position: absolute; top: 6px; height: 6px; border-radius: 3px; background: color-mix(in srgb, var(--good) 55%, transparent); }
  .mark { position: absolute; top: 1px; width: 4px; height: 16px; margin-left: -2px; border-radius: 2px; background: var(--text); box-shadow: 0 0 0 2px var(--panel); }
  .mark.out { background: var(--bad); }
  .labels span { position: absolute; top: 15px; transform: translateX(-50%); font-size: 0.7rem; color: var(--faint); white-space: nowrap; }
</style>
