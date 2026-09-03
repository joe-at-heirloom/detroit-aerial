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

# Real arterial crossings, computed from OpenStreetMap way geometry rather than
# guessed: an earlier hand-typed table was up to 900 m out, which put the label on
# a picture of the wrong junction. (scripts/roadcross.py regenerates this.)
CROSSINGS = [
    ('Ford Rd x Telegraph Rd', 42.32697, -83.27261),
    ('Ford Rd x Outer Dr', 42.32691, -83.26155),
    ('Ford Rd x Evergreen Rd', 42.32866, -83.23463),
    ('Ford Rd x Southfield Fwy', 42.32903, -83.21567),
    ('Ford Rd x Greenfield Rd', 42.32925, -83.19570),
    ('Warren Ave x Telegraph Rd', 42.34144, -83.27312),
    ('Warren Ave x Outer Dr', 42.34203, -83.25217),
    ('Warren Ave x Evergreen Rd', 42.34313, -83.23539),
    ('Warren Ave x Southfield Fwy', 42.34337, -83.21614),
    ('Warren Ave x Greenfield Rd', 42.34367, -83.19629),
    ('Joy Rd x Telegraph Rd', 42.35688, -83.27455),
    ('Joy Rd x Outer Dr', 42.35713, -83.26180),
    ('Joy Rd x Evergreen Rd', 42.35757, -83.23595),
    ('Joy Rd x Southfield Fwy', 42.35785, -83.21649),
    ('Joy Rd x Greenfield Rd', 42.35813, -83.19700),
    ('Plymouth Rd x Telegraph Rd', 42.37125, -83.27566),
    ('Plymouth Rd x Outer Dr', 42.37160, -83.26087),
    ('Plymouth Rd x Evergreen Rd', 42.37189, -83.23626),
    ('Plymouth Rd x Southfield Fwy', 42.37226, -83.21668),
    ('Plymouth Rd x Greenfield Rd', 42.37256, -83.19756),
    ('Schoolcraft Rd x Telegraph Rd', 42.38535, -83.27616),
    ('Schoolcraft Rd x Outer Dr', 42.38504, -83.25630),
    ('Fenkell Ave x Outer Dr', 42.40077, -83.23259),
    ('Fenkell Ave x Southfield Fwy', 42.40114, -83.21829),
    ('Fenkell Ave x Greenfield Rd', 42.40151, -83.19870),
    ('McNichols Rd x Telegraph Rd', 42.41454, -83.27717),
    ('McNichols Rd x Outer Dr', 42.41530, -83.23184),
    ('McNichols Rd x Southfield Fwy', 42.41555, -83.21881),
    ('McNichols Rd x Greenfield Rd', 42.41599, -83.19929),
    ('Grand River x Outer Dr', 42.40642, -83.23253),
    ('Grand River x Telegraph Rd', 42.42447, -83.27761),
]


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
    for name, lat, lon in CROSSINGS:
        if S + 0.005 < lat < N - 0.005 and W + 0.005 < lon < E - 0.005:
            a, _ = tile(tag, label, lat, lon, span, 60)
            if (a > 0).mean() > 0.9:
                cands.append((name, lat, lon))
    # spread the picks over the block: sort by latitude then take evenly spaced
    cands.sort(key=lambda c: (c[1], c[2]))
    k = min(6, len(cands)); picks = [cands[int(i * (len(cands) - 1) / max(k - 1, 1))] for i in range(k)]
    cols = 3; rows = math.ceil(len(picks) / cols)
    pad = 26
    canvas = Image.new('RGB', (cols * (px + 8) + 8, rows * (px + pad + 8) + 40), (18, 18, 18))
    d = ImageDraw.Draw(canvas)
    d.text((10, 10), f'{tag} {label}  |  left half {tag}, right half today, seam through the crossing  |  tile {span} m', fill=(240, 200, 80))
    for i, (name, lat, lon) in enumerate(picks):
        a, m = tile(tag, label, lat, lon, span, px)
        img = np.zeros((px, px), np.uint8); img[:, :px // 2] = a[:, :px // 2]; img[:, px // 2:] = m[:, px // 2:]
        im = Image.fromarray(img).convert('RGB'); dd = ImageDraw.Draw(im)
        dd.line([(px // 2, 0), (px // 2, px)], fill=(255, 170, 0), width=1)
        bar = int(px * 50 / span); dd.line([(10, px - 12), (10 + bar, px - 12)], fill=(255, 255, 255), width=3); dd.text((12, px - 26), '50 m', fill=(255, 255, 255))
        x = 8 + (i % cols) * (px + 8); y = 40 + (i // cols) * (px + pad + 8)
        canvas.paste(im, (x, y + pad)); d.text((x, y + 6), name, fill=(230, 230, 230))
    fn = os.path.join(out, f'wipe_{tag}_{label}.png'); canvas.save(fn); print('wrote', fn, len(picks), 'crossings')


if __name__ == '__main__':
    main()
