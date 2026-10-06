// DuckDB-WASM: SQL over the release's Parquet files, read with HTTP range requests (no full download).
import * as duckdb from '@duckdb/duckdb-wasm';
import eh_wasm from '@duckdb/duckdb-wasm/dist/duckdb-eh.wasm?url';
import eh_worker from '@duckdb/duckdb-wasm/dist/duckdb-browser-eh.worker.js?url';
import type { Manifest } from './data';
import { tableUrl } from './data';

let dbp: Promise<duckdb.AsyncDuckDB> | null = null;
let conn: duckdb.AsyncDuckDBConnection | null = null;
let registered = '';

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
  await c.close();
  return db;
}

export async function db(): Promise<duckdb.AsyncDuckDB> {
  dbp ??= init();
  return dbp;
}

/** Create one SQL view per table in the release (variants, coverage_bins, ...). */
export async function useRelease(m: Manifest): Promise<void> {
  if (registered === m.id) return;
  const d = await db();
  conn ??= await d.connect();
  for (const name of Object.keys(m.tables)) {
    await conn.query(`CREATE OR REPLACE VIEW ${name} AS SELECT * FROM read_parquet('${tableUrl(m, name)!}')`);
  }
  registered = m.id;
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
