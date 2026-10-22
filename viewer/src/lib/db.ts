// DuckDB-WASM: SQL over the release's Parquet files, read with HTTP range requests (no full download).
import * as duckdb from '@duckdb/duckdb-wasm';
import eh_wasm from '@duckdb/duckdb-wasm/dist/duckdb-eh.wasm?url';
import eh_worker from '@duckdb/duckdb-wasm/dist/duckdb-browser-eh.worker.js?url';
import type { Manifest } from './data';
import { tableUrls } from './data';

let dbp: Promise<duckdb.AsyncDuckDB> | null = null;
let conn: duckdb.AsyncDuckDBConnection | null = null;
let release: { m: Manifest; views: Map<string, Promise<unknown>> } | null = null;

async function init(): Promise<duckdb.AsyncDuckDB> {
  const worker = new Worker(eh_worker, { type: 'classic' });
  const db = new duckdb.AsyncDuckDB(new duckdb.VoidLogger(), worker);
  await db.instantiate(eh_wasm);
  // forceFullHTTPReads defaults to on in this build; without these flags every Parquet file is downloaded whole.
  await db.open({
    query: { castBigIntToDouble: true },
    filesystem: { forceFullHTTPReads: false, reliableHeadRequests: true, allowFullHTTPReads: true },
  });
  // Load the Parquet extension from the bundled copy (public/duckdb-ext) so the dashboard works offline.
  const c = await db.connect();
  await c.query(`SET custom_extension_repository = '${new URL('./duckdb-ext', document.baseURI).href}'`);
  // Release files never change (content-addressed), so each Parquet footer is read once, not once per query.
  await c.query('SET GLOBAL parquet_metadata_cache = true');
  await c.close();
  return db;
}

export async function db(): Promise<duckdb.AsyncDuckDB> {
  dbp ??= init();
  return dbp;
}

/** Query this release's tables (variants, coverage_bins, ...) by name. */
export async function useRelease(m: Manifest): Promise<void> {
  if (release?.m.id !== m.id) release = { m, views: new Map() };
}

// A table's view is created the first time a query names it: creating one reads its Parquet footer, a few round
// trips per file over the web, so making all of them up front held up the first page by seconds (ADR-028).
async function viewsFor(c: duckdb.AsyncDuckDBConnection, sql: string): Promise<void> {
  if (!release) return;
  const { m, views } = release;
  const catalog = /information_schema|duckdb_\w+\s*\(|\bshow\s+(all\s+)?tables\b/i.test(sql);
  for (const name of Object.keys(m.tables)) {
    if (!catalog && !new RegExp(`\\b${name}\\b`, 'i').test(sql)) continue;
    if (!views.has(name)) {
      const files = tableUrls(m, name).map((u) => `'${u}'`).join(', ');
      views.set(name, c.query(`CREATE OR REPLACE VIEW ${name} AS SELECT * FROM read_parquet([${files}])`)
        .catch((e) => { views.delete(name); throw e; }));
    }
    await views.get(name);
  }
}

export type Row = Record<string, any>;

function toPlain(v: any): any {
  if (typeof v === 'bigint') return Number(v);
  if (v && typeof v === 'object' && 'toArray' in v) return Array.from(v.toArray(), toPlain);
  return v;
}

export async function query(sql: string): Promise<{ rows: Row[]; columns: string[]; ms: number }> {
  const d = await db();
  conn ??= await d.connect();
  const t0 = performance.now();
  await viewsFor(conn, sql);
  const res = await conn.query(sql);
  const columns = res.schema.fields.map((f) => f.name);
  const rows = res.toArray().map((r: any) => {
    const o: Row = {};
    for (const c of columns) o[c] = toPlain(r[c]);
    return o;
  });
  return { rows, columns, ms: performance.now() - t0 };
}

export const sqlString = (s: string) => `'${s.replaceAll("'", "''")}'`;
