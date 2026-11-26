// The static site on Cloudflare: dist/ is the Worker's static assets, except
// its tiles. 114,218 tiles are over the 100,000 files a Worker version may
// carry, so they live in R2 -- the same bucket tiles.joelint.com serves -- and
// only /tiles/* runs this (wrangler.jsonc, run_worker_first). Everything else,
// index.html and data/ with its _headers, is served as assets without it.
const TILE_CACHE = 'public, max-age=31536000, immutable';   // a tile's path names its content

export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    if (!url.pathname.startsWith('/tiles/')) return env.ASSETS.fetch(request);
    if (request.method !== 'GET' && request.method !== 'HEAD') {
      return new Response('method not allowed\n', { status: 405, headers: { allow: 'GET, HEAD' } });
    }
    const key = url.pathname.slice('/tiles/'.length);      // b1949/16/17600/24210.webp
    const obj = request.method === 'HEAD'
      ? await env.TILES.head(key)
      : await env.TILES.get(key, { onlyIf: request.headers });
    if (obj === null) {
      // outside a layer's coverage; short, so a rebuild that adds tiles shows up
      return new Response('not found\n', { status: 404, headers: { 'cache-control': 'public, max-age=300' } });
    }
    const headers = new Headers();
    obj.writeHttpMetadata(headers);
    headers.set('etag', obj.httpEtag);
    headers.set('cache-control', TILE_CACHE);
    if (request.method === 'HEAD') {
      headers.set('content-length', String(obj.size));
      return new Response(null, { headers });
    }
    if (!('body' in obj)) return new Response(null, { status: 304, headers });   // If-None-Match hit
    return new Response(obj.body, { headers });
  },
};
