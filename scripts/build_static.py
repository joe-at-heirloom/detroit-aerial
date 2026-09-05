#!/usr/bin/env python
"""Render the survey into a folder that any static host can serve.

The local viewer reads windows out of 54 GB of GeoTIFFs on demand. GitHub Pages
and Vercel serve files, so the tiles have to exist as files first: this walks
every layer's bbox through the same renderer `serve.py` uses (hand alignments in
`data/adjust.json` included) and writes a WebP pyramid, then lays out `dist/`
with the viewer, the manifest and the places next to it.

Size is the whole question. Measured on real tiles, WebP q80 is 4-12 KB where
PNG is 30-140, so a native-resolution (z18) pyramid of the four film blocks and
downtown is about 0.8 GB -- inside what a Pages repo will hold -- and the ten
modern layers to z17 add ~2.7 GB, which wants object storage (`--tiles-base`).
Beyond native zoom Leaflet upsamples client-side, so `maxNativeZoom` is written
into the manifest and nothing is lost by capping there.

    ./.venv/bin/python scripts/build_static.py                     # film + downtown, z18
    ./.venv/bin/python scripts/build_static.py --layers all --zmax modern:17
    ./.venv/bin/python scripts/build_static.py --today remote      # Today = Esri tiles, 0 bytes
    ./.venv/bin/python scripts/build_static.py --tiles-base https://tiles.example.com

Resumable: tiles already on disk are skipped, so a killed run picks up where it
stopped and a second run only renders what changed.
"""
import argparse, io, json, math, os, shutil, sys, time
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))

FILM = {'b1949', 'b1956', 'b1961', 'b1967', 'dt1949', 'dt1956', 'dt1961'}
TODAY = {'dtnow', 'modwest'}
ESRI_IMAGERY = ('https://server.arcgisonline.com/ArcGIS/rest/services/'
                'World_Imagery/MapServer/tile/{z}/{y}/{x}')


def kind(lid):
    if lid in FILM:
        return 'film'
    if lid in TODAY:
        return 'today'
    return 'modern'


def tile_range(bbox, z):
    S, W, N, E = bbox
    n = 2 ** z
    tx = lambda lon: int((lon + 180.0) / 360.0 * n)
    ty = lambda lat: int((1.0 - math.log(math.tan(math.radians(lat)) +
                          1.0 / math.cos(math.radians(lat))) / math.pi) / 2.0 * n)
    return range(tx(W), tx(E) + 1), range(ty(N), ty(S) + 1)


# ------------------------------------------------------------------ workers
_serve = None


def _init():
    global _serve
    import serve as s
    _serve = s


def _render(job):
    """One tile. Returns bytes written (0 = empty or already there)."""
    lid, z, x, y, path, quality = job
    if os.path.exists(path):
        return 0
    from PIL import Image
    try:
        if lid == '__roads__':
            png = _serve.render_roads(z, x, y)
        else:
            adj = _serve.adjust_for(lid)
            png = (_serve.render_tile_adj(lid, z, x, y, adj) if adj
                   else _serve.render_tile(lid, z, x, y))
    except Exception as exc:
        sys.stderr.write('tile %s/%d/%d/%d: %s\n' % (lid, z, x, y, exc))
        return 0
    if not png:
        return 0
    im = Image.open(io.BytesIO(png))
    if im.mode == 'RGBA' and not im.getchannel('A').getbbox():
        return 0                                   # fully transparent
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if path.endswith('.png'):
        with open(path, 'wb') as fh:
            fh.write(png)
        return len(png)
    buf = io.BytesIO()
    im.save(buf, format='WEBP', quality=quality, method=4)
    with open(path, 'wb') as fh:
        fh.write(buf.getvalue())
    return buf.tell()


# ------------------------------------------------------------------ driver
def parse_zmax(spec, default):
    """--zmax 18 | --zmax film:18,modern:17,dt1949:19"""
    out = dict(default)
    if not spec:
        return out
    for part in spec.split(','):
        if ':' in part:
            k, v = part.split(':')
            out[k] = int(v)
        else:
            for k in out:
                out[k] = int(part)
    return out


def zmax_for(lid, zmax):
    return zmax.get(lid, zmax.get(kind(lid), 18))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(ROOT, 'dist'))
    ap.add_argument('--layers', default='film',
                    help="'film' (blocks + downtown, default), 'all', or comma list of layer ids")
    ap.add_argument('--zmin', type=int, default=8,
                    help='low zooms are a handful of tiles; the viewer fits a whole block at ~z11')
    ap.add_argument('--zmax', default='', help='e.g. 18, or film:18,modern:17,dt1949:19')
    ap.add_argument('--quality', type=int, default=80)
    ap.add_argument('--roads-max', type=int, default=17)
    ap.add_argument('--today', choices=['tiles', 'remote'], default='remote',
                    help="'remote' points Today at Esri World Imagery (no tiles written)")
    ap.add_argument('--tiles-base', default='',
                    help='URL prefix if the tiles are hosted elsewhere (still written under --out)')
    ap.add_argument('--jobs', type=int, default=max(1, os.cpu_count() - 1))
    ap.add_argument('--no-tiles', action='store_true', help='lay out dist/ without rendering')
    a = ap.parse_args()

    man = json.load(open(os.path.join(ROOT, 'data', 'manifest.json')))
    zmax = parse_zmax(a.zmax, {'film': 18, 'modern': 17, 'today': 17})
    ids = [l['id'] for l in man['layers']]
    if a.layers == 'film':
        want = [i for i in ids if i in FILM or (i in TODAY and a.today == 'tiles')]
    elif a.layers == 'all':
        want = [i for i in ids if not (i in TODAY and a.today == 'remote')]
    else:
        want = [i for i in a.layers.split(',') if i in ids]

    out = a.out
    tiles_dir = os.path.join(out, 'tiles')
    os.makedirs(tiles_dir, exist_ok=True)

    # ---- the manifest the static viewer will read
    layers = []
    for l in man['layers']:
        lid = l['id']
        m = {k: v for k, v in l.items() if k not in ('file',)}
        if lid in TODAY and a.today == 'remote':
            m['url'] = ESRI_IMAGERY
            m['maxNativeZoom'] = 19
        elif lid in want:
            m['maxNativeZoom'] = zmax_for(lid, zmax)
            m['minNativeZoom'] = a.zmin
        else:
            continue                                     # not built: not offered
        layers.append(m)
    kept = {m['id'] for m in layers}
    groups = {}
    for gid, g in man['groups'].items():
        gl = [i for i in g['layers'] if i in kept]
        if gl:
            groups[gid] = dict(g, layers=gl)
    static_man = dict(layers=layers, groups=groups, static=True,
                      roads_max=a.roads_max, serve=man.get('serve', {}),
                      built=time.strftime('%Y-%m-%d %H:%M'))
    os.makedirs(os.path.join(out, 'data'), exist_ok=True)
    json.dump(static_man, open(os.path.join(out, 'data', 'manifest.json'), 'w'), indent=1)
    shutil.copy(os.path.join(ROOT, 'data', 'places.json'), os.path.join(out, 'data', 'places.json'))

    # ---- the viewer, told where it is
    html = open(os.path.join(ROOT, 'viewer', 'index.html'), encoding='utf-8').read()
    cfg = dict(tiles=(a.tiles_base.rstrip('/') or 'tiles'), ext='webp',
               data='data', api=None, static=True)
    inject = '<script>window.DAS = %s;</script>\n' % json.dumps(cfg)
    marker = '<script src="https://unpkg.com/leaflet'
    if marker not in html:
        raise SystemExit('viewer/index.html: cannot find where to inject config')
    html = html.replace(marker, inject + marker, 1)
    open(os.path.join(out, 'index.html'), 'w', encoding='utf-8').write(html)
    open(os.path.join(out, '.nojekyll'), 'w').close()      # Pages: serve dotfiles/underscores as-is

    if a.no_tiles:
        print('laid out %s without tiles' % out)
        return

    # ---- the jobs
    jobs = []
    for m in layers:
        lid = m['id']
        if 'url' in m:
            continue
        for z in range(a.zmin, m['maxNativeZoom'] + 1):
            xs, ys = tile_range(m['bbox'], z)
            for x in xs:
                for y in ys:
                    jobs.append((lid, z, x, y,
                                 os.path.join(tiles_dir, lid, str(z), str(x), '%d.webp' % y),
                                 a.quality))
    # roads over everything the viewer can show, capped: beyond it Leaflet scales the lines
    union = None
    for g in groups.values():
        b = g['bbox']
        union = b if union is None else [min(union[0], b[0]), min(union[1], b[1]),
                                          max(union[2], b[2]), max(union[3], b[3])]
    if union:
        for z in range(14, a.roads_max + 1):
            xs, ys = tile_range(union, z)
            for x in xs:
                for y in ys:
                    jobs.append(('__roads__', z, x, y,
                                 os.path.join(tiles_dir, '__roads__', str(z), str(x), '%d.png' % y),
                                 a.quality))
    print('%d tiles to consider across %d layers, %d workers' % (len(jobs), len(layers), a.jobs), flush=True)

    t0 = time.time()
    done = written = nbytes = 0
    per = {}
    with Pool(a.jobs, initializer=_init) as pool:
        for job, n in zip(jobs, pool.imap(_render, jobs, chunksize=64)):
            done += 1
            if n:
                written += 1; nbytes += n
                per[job[0]] = per.get(job[0], 0) + n
            if done % 5000 == 0 or done == len(jobs):
                el = time.time() - t0
                print('  %6d / %d   %5.1f min   %.2f GB' % (done, len(jobs), el / 60, nbytes / 1e9), flush=True)

    # ---- what is on disk now (including tiles from earlier runs)
    print('\nlayer sizes on disk:')
    total = 0
    for lid in sorted(os.listdir(tiles_dir)):
        d = os.path.join(tiles_dir, lid)
        sz = sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(d) for f in fs)
        total += sz
        print('  %-10s %7.1f MB' % (lid, sz / 1e6))
    print('  %-10s %7.1f MB   (%d tiles written this run in %.1f min)'
          % ('TOTAL', total / 1e6, written, (time.time() - t0) / 60))
    print('\nserve it:  python3 -m http.server -d %s 8780' % out)


if __name__ == '__main__':
    main()
