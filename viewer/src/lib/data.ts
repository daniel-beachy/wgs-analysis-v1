// Release discovery: the launcher serves the workspace `releases/` folder at ./data/.

export interface ObjectRef { path: string; bytes: number; sha256: string; rows?: number; description?: string }
export interface ReleaseEntry { id: string; sample: string; created: string; digest: string; manifest: string; previous: string | null; summary?: string; knowledge?: Record<string, string> }
export interface Index { schema: number; releases: ReleaseEntry[]; latest: string }
export interface Manifest {
  schema: number; id: string; sample: string; created: string; digest: string; previous: string | null;
  pipeline_version: string;
  tables: Record<string, ObjectRef>;
  documents: Record<string, ObjectRef>;
  knowledge: { versions?: Record<string, string>; evidence_model?: string | null };
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

export function tableUrl(m: Manifest, name: string): string | null {
  const ref = m.tables[name];
  return ref ? DATA + ref.path : null;
}

// Tell the launcher the tab is still open; it exits a while after the last tab closes.
export function startHeartbeat() {
  const beat = () => fetch('./api/heartbeat', { method: 'POST' }).catch(() => {});
  beat();
  setInterval(beat, 10_000);
  window.addEventListener('pagehide', () => navigator.sendBeacon?.('./api/bye'));
}
