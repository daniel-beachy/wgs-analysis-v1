<script lang="ts">
  import type { Manifest } from '../lib/data';
  import { query, useRelease, type Row } from '../lib/db';
  import Card from '../components/Card.svelte';
  import Explain from '../components/Explain.svelte';

  let { m }: { m: Manifest } = $props();
  const EXAMPLES: [string, string][] = [
    ['Variants per chromosome', "SELECT chrom, count(*) AS variants, count(*) FILTER (WHERE zygosity = 'hom') AS hom\nFROM variants WHERE filter = 'PASS'\nGROUP BY chrom, chrom_order ORDER BY chrom_order"],
    ['Lactase persistence (rs4988235)', "SELECT * FROM variants WHERE rsid = 'rs4988235'"],
    ['Largest deletions', "SELECT chrom, pos, length(ref) - length(alt) AS deleted_letters, zygosity, gq\nFROM variants WHERE filter = 'PASS' AND vtype = 'DEL'\nORDER BY deleted_letters DESC LIMIT 20"],
    ['Low-coverage stretches on chr1', "SELECT chrom, start, \"end\", \"end\" - start AS bp, state\nFROM callable WHERE chrom = '1' AND state <> 'CALLABLE'\nORDER BY bp DESC LIMIT 20"],
    ['Tables and columns', "SELECT table_name, column_name, data_type FROM information_schema.columns ORDER BY table_name, ordinal_position"],
  ];
  let sql = $state(EXAMPLES[0][1]);
  let rows = $state<Row[]>([]);
  let cols = $state<string[]>([]);
  let ms = $state(0);
  let error = $state('');
  let busy = $state(false);

  async function run() {
    busy = true; error = '';
    try {
      await useRelease(m);
      const r = await query(sql.trim().replace(/;$/, ''));
      rows = r.rows.slice(0, 1000); cols = r.columns; ms = r.ms;
    } catch (e) { error = String(e).replace(/^Error: /, ''); rows = []; cols = []; }
    busy = false;
  }
  const show = (v: any) => (v == null ? '' : typeof v === 'number' && !Number.isInteger(v) ? v.toFixed(4) : String(v));
</script>

<div class="grid">
  <Card title="Ask your own questions (SQL)" subtitle="The whole release is queryable here with DuckDB, running entirely in your browser — nothing leaves this computer.">
    <Explain title="New to SQL?">
      <p>SQL is a language for asking questions of tables: <code>SELECT</code> the columns you want <code>FROM</code> a table, keep rows <code>WHERE</code> a condition holds, and <code>GROUP BY</code> / <code>ORDER BY</code> to summarise and sort. Start from an example below and change one thing at a time. Tables in this release: {Object.keys(m.tables).map((t) => t).join(', ')}.</p>
    </Explain>
    <div class="row examples">{#each EXAMPLES as [label, q]}<button class="small" onclick={() => { sql = q; run(); }}>{label}</button>{/each}</div>
    <textarea bind:value={sql} rows="6" spellcheck="false" onkeydown={(e) => { if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') run(); }}></textarea>
    <div class="row"><button class="primary" onclick={run} disabled={busy}>{busy ? 'Running…' : 'Run'}</button><span class="faint small">⌘/Ctrl + Enter · {rows.length ? `${rows.length} rows · ${ms.toFixed(0)} ms` : ''}</span></div>
    {#if error}<pre class="err">{error}</pre>{/if}
    {#if cols.length}
      <div class="out">
        <table class="small mono">
          <thead><tr>{#each cols as c}<th>{c}</th>{/each}</tr></thead>
          <tbody>{#each rows as r}<tr>{#each cols as c}<td>{show(r[c])}</td>{/each}</tr>{/each}</tbody>
        </table>
      </div>
    {/if}
  </Card>
</div>

<style>
  textarea { width: 100%; font: 0.88rem/1.45 ui-monospace, Menlo, Consolas, monospace; background: var(--code); color: var(--text); border: 1px solid var(--border); border-radius: 8px; padding: 10px; margin: 10px 0; resize: vertical; }
  .examples { margin-top: 6px; }
  .out { overflow: auto; max-height: 560px; margin-top: 10px; }
  .err { color: var(--bad); white-space: pre-wrap; background: var(--code); padding: 10px; border-radius: 8px; }
</style>
