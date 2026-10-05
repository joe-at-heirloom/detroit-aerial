#!/usr/bin/env python
"""One spot, every year: a contact sheet of what the year view shows there.

For each year in data/timeline.json, the picture on top at the point -- the
newest layer at or before that year whose outline holds it, exactly the one the
year view draws -- cut to the same ground window and captioned with its year and
publisher. The fastest way to see whether the years line up, and what changed.

    ./.venv/bin/python scripts/yearsheet.py 42.452 -83.206 [--km 2.4] [--out sheet.jpg]
"""
import sys, os, json, io, math, urllib.request
import numpy as np
import rasterio, rasterio.transform
from rasterio.windows import from_bounds
from PIL import Image, ImageDraw, ImageFont
from shapely.geometry import shape, Point, box
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import dtmap
def P(*a): return os.path.join(ROOT, *a)
def arg(n, d=None): return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d

CELL = 300          # px per year


def remote_window(url, S, W, N, E, z=16):
    """A remote layer (Esri Wayback) at its own true zoom: a handful of tiles."""
    n = 2 ** z
    tx = lambda lon: (lon + 180) / 360 * n
    ty = lambda lat: (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n
    x0, x1 = int(tx(W)), int(tx(E)); y0, y1 = int(ty(N)), int(ty(S))
    mos = Image.new('RGB', ((x1 - x0 + 1) * 256, (y1 - y0 + 1) * 256))
    for x in range(x0, x1 + 1):
        for y in range(y0, y1 + 1):
            u = url.replace('{z}', str(z)).replace('{x}', str(x)).replace('{y}', str(y))
            try:
                rq = urllib.request.Request(u, headers={'User-Agent': 'detroit-air-survey/1.0'})
                mos.paste(Image.open(io.BytesIO(urllib.request.urlopen(rq, timeout=60).read())).convert('RGB'),
                          ((x - x0) * 256, (y - y0) * 256))
            except Exception:
                pass
    px = lambda v0, v: (v - v0) * 256
    return mos.crop((int(px(x0, tx(W))), int(px(y0, ty(N))), int(px(x0, tx(E))), int(px(y0, ty(S)))))


def main():
    lat, lon = float(sys.argv[1]), float(sys.argv[2])
    km = float(arg('--km', 2.4))
    dlat = km * 500 / dtmap.MLAT; dlon = km * 500 / dtmap.MLON
    S, W, N, E = lat - dlat, lon - dlon, lat + dlat, lon + dlon
    tl = json.load(open(P('data', 'timeline.json')))
    man = {l['id']: l for l in json.load(open(P('data', 'manifest.json')))['layers']}
    pt = Point(lon, lat)
    win = box(W, S, E, N)
    cache = {}

    def picture(lid):
        if lid not in cache:
            l = man[lid]
            if l.get('url'):
                im = np.asarray(remote_window(l['url'], S, W, N, E).resize((CELL, CELL)))
            else:
                with rasterio.open(P(l['file'])) as ds:
                    bands = [1, 2, 3] if ds.count >= 3 else [1, 1, 1]
                    # a plain image (downtown's Today JPEG) has no georeferencing of
                    # its own: it is placed by its manifest bbox, as serve.py does
                    tr = ds.transform
                    if ds.crs is None or tr.is_identity:
                        b = l['bbox']
                        tr = rasterio.transform.from_bounds(b[1], b[0], b[3], b[2], ds.width, ds.height)
                    a = ds.read(bands, window=from_bounds(W, S, E, N, tr),
                                out_shape=(3, CELL, CELL), boundless=True, fill_value=0)
                im = np.transpose(a, (1, 2, 0))
            cache[lid] = im
        return cache[lid]

    cells = []
    for step in tl['years']:
        top = next((s['id'] for s in step['stack']
                    if shape(tl['layers'][s['id']]['footprint']).contains(pt)), None)
        if top is None:
            continue
        # what the year view draws: newest first, each layer filling only what is
        # still empty, until nothing is
        canvas = np.zeros((CELL, CELL, 3), np.uint8)
        for s in step['stack']:
            if not shape(tl['layers'][s['id']]['footprint']).intersects(win):
                continue
            hole = canvas.max(axis=2) == 0
            if not hole.any():
                break
            im = picture(s['id'])
            fill = hole & (im.max(axis=2) > 0)
            canvas[fill] = im[fill]
        meta = tl['layers'][top]
        cells.append((step['label'] if meta['year'] == step['year'] else f"{step['label']}: {meta['label']}",
                      meta.get('credit', ''), Image.fromarray(canvas)))
        print(f"  {step['label']:>6}  {top:<10} {meta.get('credit', '')}", flush=True)
    cols = int(arg('--cols', 7)); rows = math.ceil(len(cells) / cols)
    sheet = Image.new('RGB', (cols * CELL, rows * (CELL + 34)), (11, 15, 18))
    d = ImageDraw.Draw(sheet)
    try:
        f1 = ImageFont.truetype('/System/Library/Fonts/Menlo.ttc', 15)
        f2 = ImageFont.truetype('/System/Library/Fonts/Menlo.ttc', 10)
    except OSError:
        f1 = f2 = ImageFont.load_default()
    for i, (yr, cr, im) in enumerate(cells):
        x, y = (i % cols) * CELL, (i // cols) * (CELL + 34)
        sheet.paste(im, (x, y + 34))
        d.text((x + 6, y + 3), yr, fill=(224, 169, 74), font=f1)
        d.text((x + 6, y + 20), cr[:46], fill=(124, 143, 155), font=f2)
    out = arg('--out', P('logs', f'yearsheet_{lat:.4f}_{lon:.4f}.jpg'))
    sheet.save(out, quality=88)
    print(f"  wrote {out}: {len(cells)} years, {km} km window")


if __name__ == '__main__':
    main()
