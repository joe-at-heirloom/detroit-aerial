#!/usr/bin/env python
"""A build wiped against today at mile-road crossings.

Each tile is a square of ground centred on an arterial intersection: the left
half is the historical build, the right half is modern imagery, and the seam is
a vertical line through the crossing. A road that runs straight through the
seam is a road the build has put in the right place; a jog at the seam is the
placement error, readable directly in metres from the scale bar.

Usage:  ./.venv/bin/python scripts/wipepic.py 1949 placedC [--span 500] [--out /tmp]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import rasterio
from rasterio.windows import from_bounds
from PIL import Image, ImageDraw
from validate import P, modern_for

# Detroit's mile roads on the west side (WGS84), so the crossings are real ones.
EW = {'Ford Rd': 42.3325, 'Warren Ave': 42.3468, 'Joy Rd': 42.3598, 'Plymouth Rd': 42.3737,
      'Schoolcraft Rd': 42.3875, 'Fenkell Ave': 42.4008, 'McNichols Rd': 42.4155, 'Seven Mile Rd': 42.4295}
NS = {'Telegraph Rd': -83.2735, 'Evergreen Rd': -83.2445, 'Southfield Fwy': -83.2225,
      'Greenfield Rd': -83.2010, 'Outer Dr': -83.2300, 'Wyoming Ave': -83.1650}


def tile(tag, label, lat, lon, span, px):
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{label}.tif'))
    dlat = span / 111320.0; dlon = span / (111320.0 * math.cos(math.radians(lat)))
    s, w, n, e = lat - dlat / 2, lon - dlon / 2, lat + dlat / 2, lon + dlon / 2
    a = ds.read(1, window=from_bounds(w, s, e, n, ds.transform), out_shape=(px, px), boundless=True, fill_value=0)
    m = modern_for((s, w, n, e), px, px)
    if m.ndim == 3: m = m.mean(axis=2)
    return a.astype(np.uint8), np.clip(m, 0, 255).astype(np.uint8)


def main():
    tag, label = sys.argv[1], sys.argv[2]
    span = int(sys.argv[sys.argv.index('--span') + 1]) if '--span' in sys.argv else 500
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else '/tmp'
    g = json.load(open(P('data', f'{tag}_{label}_geo.json'))); S, W, N, E = g['bbox']
    px = 420
    cands = []
    for en, lat in EW.items():
        for nn, lon in NS.items():
            if S + 0.01 < lat < N - 0.01 and W + 0.01 < lon < E - 0.01:
                a, _ = tile(tag, label, lat, lon, span, 60)
                if (a > 0).mean() > 0.9:
                    cands.append((en, nn, lat, lon))
    # spread the picks over the block: sort by latitude then take evenly spaced
    cands.sort(key=lambda c: (c[2], c[3]))
    k = min(6, len(cands)); picks = [cands[int(i * (len(cands) - 1) / max(k - 1, 1))] for i in range(k)]
    cols = 3; rows = math.ceil(len(picks) / cols)
    pad = 26
    canvas = Image.new('RGB', (cols * (px + 8) + 8, rows * (px + pad + 8) + 40), (18, 18, 18))
    d = ImageDraw.Draw(canvas)
    d.text((10, 10), f'{tag} {label}  |  left half {tag}, right half today, seam through the crossing  |  tile {span} m', fill=(240, 200, 80))
    for i, (en, nn, lat, lon) in enumerate(picks):
        a, m = tile(tag, label, lat, lon, span, px)
        img = np.zeros((px, px), np.uint8); img[:, :px // 2] = a[:, :px // 2]; img[:, px // 2:] = m[:, px // 2:]
        im = Image.fromarray(img).convert('RGB'); dd = ImageDraw.Draw(im)
        dd.line([(px // 2, 0), (px // 2, px)], fill=(255, 170, 0), width=1)
        bar = int(px * 50 / span); dd.line([(10, px - 12), (10 + bar, px - 12)], fill=(255, 255, 255), width=3); dd.text((12, px - 26), '50 m', fill=(255, 255, 255))
        x = 8 + (i % cols) * (px + 8); y = 40 + (i // cols) * (px + pad + 8)
        canvas.paste(im, (x, y + pad)); d.text((x, y + 6), f'{en} x {nn}', fill=(230, 230, 230))
    fn = os.path.join(out, f'wipe_{tag}_{label}.png'); canvas.save(fn); print('wrote', fn, len(picks), 'crossings')


if __name__ == '__main__':
    main()
