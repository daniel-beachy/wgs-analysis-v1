// The primer (docs/GENETICS_PRIMER.md) is bundled so the dashboard can explain itself offline.
// Its glossary table is the single source for every tooltip in the app.
import primer from '../../../docs/GENETICS_PRIMER.md?raw';
import zoomSvg from '../../../docs/img/genome-zoom.svg?url';

export const PRIMER: string = primer.replaceAll('img/genome-zoom.svg', zoomSvg);

export const slug = (s: string) =>
  s.toLowerCase().replace(/<[^>]+>/g, '').replace(/[`*_]/g, '').replace(/[^a-z0-9 -]/g, '').trim().replace(/\s+/g, '-');

export interface Gloss { term: string; names: string[]; meaning: string }

export const GLOSSARY: Gloss[] = (() => {
  const start = PRIMER.indexOf('## Glossary');
  const out: Gloss[] = [];
  if (start < 0) return out;
  for (const line of PRIMER.slice(start).split('\n')) {
    const m = /^\|\s*\*\*(.+?)\*\*\s*\|\s*(.+?)\s*\|\s*$/.exec(line);
    if (!m) continue;
    const term = m[1];
    const names = term.split(/\s*\/\s*|\s*\(|\)/).map((s) => s.trim().toLowerCase()).filter(Boolean);
    out.push({ term, names: [term.toLowerCase(), ...names], meaning: m[2].replace(/`/g, '') });
  }
  return out;
})();

export function lookup(t: string): Gloss | undefined {
  const k = t.toLowerCase();
  return GLOSSARY.find((g) => g.names.includes(k)) ?? GLOSSARY.find((g) => g.names.some((n) => n.startsWith(k)));
}

export const SECTIONS = Array.from(PRIMER.matchAll(/^(#{2,3}) (.+)$/gm), (m) => ({
  level: m[1].length,
  title: m[2].replace(/[`*]/g, ''),
  slug: slug(m[2]),
}));

/** Find a primer section whose title contains `needle` (case-insensitive). */
export function sectionFor(needle: string): string | undefined {
  const n = needle.toLowerCase();
  return SECTIONS.find((s) => s.title.toLowerCase().includes(n))?.slug;
}
