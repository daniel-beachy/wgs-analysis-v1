import { defineConfig, type Plugin } from 'vite';
import { svelte } from '@sveltejs/vite-plugin-svelte';
import { createReadStream, existsSync, statSync } from 'node:fs';
import { resolve, join, normalize } from 'node:path';

// In `npm run dev`, serve the workspace releases at /data like the launcher does.
function devData(): Plugin {
  const root = resolve(process.env.WGS_RELEASES ?? resolve(import.meta.dirname, '../../wgs-data/releases'));
  return {
    name: 'wgs-dev-data',
    configureServer(server) {
      // The launcher implements /api/* (heartbeat for auto-exit); in dev just acknowledge.
      server.middlewares.use('/api', (_req, res) => { res.statusCode = 204; res.end(); });
      server.middlewares.use('/data', (req, res, next) => {
        const rel = normalize(decodeURIComponent((req.url ?? '/').split('?')[0])).replace(/^([/\\])+/, '');
        const file = join(root, rel);
        if (!file.startsWith(root) || !existsSync(file) || statSync(file).isDirectory()) return next();
        const size = statSync(file).size;
        const range = /bytes=(\d+)-(\d*)/.exec(req.headers.range ?? '');
        res.setHeader('Accept-Ranges', 'bytes');
        if (range) {
          const start = +range[1];
          const end = range[2] ? +range[2] : size - 1;
          res.statusCode = 206;
          res.setHeader('Content-Range', `bytes ${start}-${end}/${size}`);
          res.setHeader('Content-Length', end - start + 1);
          createReadStream(file, { start, end }).pipe(res);
        } else {
          res.setHeader('Content-Length', size);
          createReadStream(file).pipe(res);
        }
      });
    },
  };
}

export default defineConfig({
  base: './',
  plugins: [svelte(), devData()],
  server: { fs: { allow: ['..'] } },
  build: { outDir: process.env.WGS_VIEWER_OUT ?? 'dist', emptyOutDir: true, chunkSizeWarningLimit: 4000 },
  worker: { format: 'es' },
});
