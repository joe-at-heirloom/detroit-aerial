#!/usr/bin/env python
"""The same ground in every year we have, side by side.

Given a point, renders a row of square crops -- one per block that covers it,
plus modern imagery -- at the same scale and extent, so a site can be followed
forward and backward in time. The build shown for each year is whatever the
manifest currently serves, so the picture never disagrees with the viewer.

Usage:  ./.venv/bin/python scripts/sitepic.py 42.3598 -83.2225 [--span 600] [--out /tmp] [--name joy-southfield]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import rasterio
from rasterio.windows import from_bounds
from PIL import Image, ImageDraw
from validate import P, modern_for


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def served():
    """(year, mosaic path) for each West Detroit block, in the manifest's own order."""
    m = json.load(open(P('data', 'manifest.json')))
    out = []
    for l in m['layers']:
        if l.get('group') == 'west' and l['id'].startswith('b'):
            out.append((l['id'][1:], os.path.join(ROOT, l['file'])))
    return sorted(out)


def crop(path, lat, lon, span, px):
    ds = rasterio.open(path)
    dlat = span / 111320.0; dlon = span / (111320.0 * math.cos(math.radians(lat)))
    s, w, n, e = lat - dlat / 2, lon - dlon / 2, lat + dlat / 2, lon + dlon / 2
    a = ds.read(1, window=from_bounds(w, s, e, n, ds.transform), out_shape=(px, px),
                boundless=True, fill_value=0)
    return a.astype(np.uint8), (s, w, n, e)


def main():
    lat, lon = float(sys.argv[1]), float(sys.argv[2])
    span = float(arg('--span', 600)); out = arg('--out', '/tmp'); name = arg('--name', f'{lat:.4f}_{lon:.4f}')
    px = 460; pad = 26
    panels = []
    for year, path in served():
        a, bb = crop(path, lat, lon, span, px)
        if (a > 0).mean() < 0.5:
            continue
        panels.append((year, a))
    m = modern_for((lat - span / 111320.0 / 2, lon - span / (111320.0 * math.cos(math.radians(lat))) / 2,
                    lat + span / 111320.0 / 2, lon + span / (111320.0 * math.cos(math.radians(lat))) / 2), px, px)
    if m.ndim == 3: m = m.mean(axis=2)
    panels.append(('today', np.clip(m, 0, 255).astype(np.uint8)))

    W = len(panels) * (px + 8) + 8
    canvas = Image.new('RGB', (W, px + pad + 48), (18, 18, 18))
    d = ImageDraw.Draw(canvas)
    d.text((10, 10), f'{lat:.5f}, {lon:.5f}   |   {span:.0f} m across each frame', fill=(240, 200, 80))
    for i, (year, a) in enumerate(panels):
        im = Image.fromarray(a).convert('RGB'); dd = ImageDraw.Draw(im)
        # centre cross, so the same point is identifiable in every frame
        c = px // 2
        for (x0, y0, x1, y1) in ((c - 14, c, c - 5, c), (c + 5, c, c + 14, c),
                                 (c, c - 14, c, c - 5), (c, c + 5, c, c + 14)):
            dd.line([(x0, y0), (x1, y1)], fill=(255, 170, 0), width=2)
        bar = int(px * 100 / span)
        dd.line([(12, px - 14), (12 + bar, px - 14)], fill=(255, 255, 255), width=3)
        dd.text((14, px - 30), '100 m', fill=(255, 255, 255))
        x = 8 + i * (px + 8)
        canvas.paste(im, (x, 40 + pad))
        d.text((x, 44), year.upper(), fill=(235, 235, 235))
    fn = os.path.join(out, f'site_{name}.png'); canvas.save(fn)
    print('wrote', fn, '-', ', '.join(y for y, _ in panels))


if __name__ == '__main__':
    main()
