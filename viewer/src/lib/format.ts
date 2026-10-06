export const nf = new Intl.NumberFormat('en-US');
export const int = (n: number | null | undefined) => (n == null ? '—' : nf.format(Math.round(n)));
export const pct = (f: number | null | undefined, digits = 1) => (f == null ? '—' : `${(f * 100).toFixed(digits)}%`);
export function bytes(n: number): string {
  const u = ['B', 'KB', 'MB', 'GB', 'TB'];
  let i = 0;
  while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
  return `${n.toFixed(i ? 1 : 0)} ${u[i]}`;
}
export function compact(n: number): string {
  if (n >= 1e9) return (n / 1e9).toFixed(2).replace(/\.?0+$/, '') + ' B';
  if (n >= 1e6) return (n / 1e6).toFixed(2).replace(/\.?0+$/, '') + ' M';
  if (n >= 1e3) return (n / 1e3).toFixed(1).replace(/\.0$/, '') + ' k';
  return String(n);
}
export const dateTime = (iso: string) => new Date(iso).toLocaleString(undefined, { dateStyle: 'medium', timeStyle: 'short' });
