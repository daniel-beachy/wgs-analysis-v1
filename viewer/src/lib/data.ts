// Release discovery: the launcher serves the workspace `releases/` folder at ./data/.

export interface ObjectRef {
  path: string; bytes: number; sha256: string; rows?: number; description?: string;
  /** A table split by chromosome for a static web host (`wgs export-site`, ADR-028): read every part. */
  parts?: ObjectRef[];
  /** The rows (and parts) are stored in this order, so a query in this order needs no ORDER BY (ADR-028). */
  sorted_by?: string[];
}
export interface ReleaseEntry { id: string; sample: string; created: string; digest: string; manifest: string; previous: string | null; summary?: string; knowledge?: Record<string, string> }
export interface Index { schema: number; releases: ReleaseEntry[]; latest: string }
export interface Manifest {
  schema: number; id: string; sample: string; created: string; digest: string; previous: string | null;
  pipeline_version: string;
  tables: Record<string, ObjectRef>;
  documents: Record<string, ObjectRef>;
  knowledge: { versions?: Record<string, string>; evidence_model?: string | null };
  /** Set when the release was built from a public truth set (the demo, Step 9c), not a personal genome. */
  truth_set?: { label: string; version?: string; url?: string };
  /** Set on the public web copy of the demo: what was left out or split (ADR-028). */
  export?: { created: string; from_release: string; dropped_columns: Record<string, Record<string, string>>; split: Record<string, number> };
}

export const DATA = new URL('./data/', document.baseURI).href;

async function getJSON<T>(rel: string): Promise<T> {
  const r = await fetch(DATA + rel, { cache: 'no-cache' });
  if (!r.ok) throw new Error(`${r.status} ${r.statusText} for data/${rel}`);
  return r.json() as Promise<T>;
}

export const loadIndex = () => getJSON<Index>('index.json');
export const loadManifest = (entry: ReleaseEntry) => getJSON<Manifest>(entry.manifest);

const docCache = new Map<string, Promise<any>>();
export function loadDoc<T = any>(m: Manifest, name: string): Promise<T | null> {
  const ref = m.documents[name];
  if (!ref) return Promise.resolve(null);
  if (!docCache.has(ref.path)) docCache.set(ref.path, getJSON(ref.path));
  return docCache.get(ref.path)!;
}

/** The latest `wgs audit` of this release against live sources, if one has been run (not part of the release). */
export async function loadAudit(m: Manifest): Promise<any | null> {
  try { return await getJSON(`audits/${m.id}.json`); } catch { return null; }
}

export function tableUrls(m: Manifest, name: string): string[] {
  const ref = m.tables[name];
  return ref ? (ref.parts ?? [ref]).map((p) => DATA + p.path) : [];
}

// Tell the launcher the tab is still open; it exits a while after the last tab closes. A static web host has no
// launcher, so the first refused beat stops the timer.
export function startHeartbeat() {
  let timer = 0;
  const beat = () => fetch('./api/heartbeat', { method: 'POST' })
    .then((r) => { if (!r.ok) { clearInterval(timer); window.removeEventListener('pagehide', bye); } })
    .catch(() => {});
  const bye = () => navigator.sendBeacon?.('./api/bye');
  timer = window.setInterval(beat, 10_000);
  window.addEventListener('pagehide', bye);
  beat();
}
