// "Reveal sensitive findings" is remembered on this device (the owner asked to see everything), and falls back to
// hidden whenever browser storage is unavailable.
const KEY = 'reveal-sensitive';

export function getReveal(): boolean {
  try { return localStorage.getItem(KEY) === '1' || sessionStorage.getItem(KEY) === '1'; } catch { return false; }
}

export function setReveal(on: boolean): void {
  try { localStorage.setItem(KEY, on ? '1' : '0'); sessionStorage.setItem(KEY, on ? '1' : '0'); } catch { /* ignore */ }
}
