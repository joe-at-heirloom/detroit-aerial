#!/usr/bin/env python
"""Visual proof, spread over the whole block rather than cherry-picked.

Windows are chosen on a lattice over the covered area, rendered as
historical-in-red over modern-in-green, and tiled into one sheet. Where the two
agree the streets are yellow and single; where they do not, every street shows a
red ghost beside a green one. Areas that are red-only or green-only are land-use
change -- a freeway that did not exist in 1961 has nothing to align to.

Coverage is found on a coarse overview and only the chosen windows are read at
full resolution, so this stays cheap on an 800 megapixel mosaic.

Usage:  ./.venv/bin/python scripts/montage.py 1961 [--suffix final] [--span 700]
        [--cols 4] [--n 12] [--tilepx 560] [--out sheet.png] [--solo]

`--solo` renders the historical mosaic on its own, which is how you see stitching
seams: per-frame correction can open a step where two negatives meet, and the
red/green overlay hides it.
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image, ImageDraw
import rasterio
from rasterio.windows import Window
import dtmap
import validate
from validate import P
Image.MAX_IMAGE_PIXELS = None


def arg(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


def stretch(a):
    a = a.astype(np.float32); v = a > 0
    if v.sum() < 10:
        return np.zeros_like(a)
    lo, hi = np.percentile(a[v], [2, 98])
    o = np.clip((a - lo) / max(hi - lo, 1) * 255, 0, 255); o[~v] = 0
    return o


def main():
    tag = [a for a in sys.argv[1:] if not a.startswith('--')][0]
    suffix = arg('--suffix', 'final')
    span = float(arg('--span', 700))
    cols = int(arg('--cols', 4))
    n = int(arg('--n', 12))
    tilepx = int(arg('--tilepx', 560))
    out = arg('--out', f'/tmp/montage_{tag}_{suffix}.png')
    solo = '--solo' in sys.argv   # historical only, to look for seams

    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json')))
    S, W, N, E = geo['bbox']
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{suffix}.tif'))

    # coverage on a coarse overview
    SCAN = 25.0
    sw = max(1, int((E - W) * dtmap.MLON / SCAN)); sh = max(1, int((N - S) * dtmap.MLAT / SCAN))
    cov = ds.read(1, out_shape=(sh, sw)) > 0

    rows = int(math.ceil(n / cols))
    half_deg_lat = span / 2 / dtmap.MLAT
    half_deg_lon = span / 2 / dtmap.MLON
    picks = []
    ny, nx = rows * 3, cols * 2
    for j in range(ny):
        for i in range(nx):
            lat = N - (j + 0.5) / ny * (N - S)
            lon = W + (i + 0.5) / nx * (E - W)
            if not (S + half_deg_lat < lat < N - half_deg_lat
                    and W + half_deg_lon < lon < E - half_deg_lon):
                continue
            cy = int((N - lat) / (N - S) * sh); cx = int((lon - W) / (E - W) * sw)
            r = max(2, int(span / 2 / SCAN))
            if cov[max(0, cy - r):cy + r, max(0, cx - r):cx + r].mean() > 0.97:
                picks.append((lat, lon))
    if not picks:
        print("no fully covered windows"); return
    step = max(1, len(picks) // n)
    picks = picks[::step][:n]

    sheet = np.zeros((rows * tilepx, cols * tilepx, 3), np.uint8)
    for k, (lat, lon) in enumerate(picks):
        r, c = divmod(k, cols)
        if r >= rows:
            break
        bb = [lat - half_deg_lat, lon - half_deg_lon, lat + half_deg_lat, lon + half_deg_lon]
        col0 = (bb[1] - W) / (E - W) * ds.width
        row0 = (N - bb[2]) / (N - S) * ds.height
        wpx = (bb[3] - bb[1]) / (E - W) * ds.width
        hpx = (bb[2] - bb[0]) / (N - S) * ds.height
        a = ds.read(1, window=Window(col0, row0, wpx, hpx),
                    out_shape=(tilepx, tilepx), boundless=True, fill_value=0)
        rgb = np.zeros((tilepx, tilepx, 3), np.uint8)
        if solo:
            g = stretch(a).astype(np.uint8)
            rgb[..., 0] = g; rgb[..., 1] = g; rgb[..., 2] = g
        else:
            rgb[..., 0] = stretch(a).astype(np.uint8)
            rgb[..., 1] = stretch(validate.modern_for(bb, tilepx, tilepx)).astype(np.uint8)
        sheet[r * tilepx:(r + 1) * tilepx, c * tilepx:(c + 1) * tilepx] = rgb
    img = Image.fromarray(sheet)
    d = ImageDraw.Draw(img)
    for k, (lat, lon) in enumerate(picks):
        r, c = divmod(k, cols)
        if r >= rows:
            break
        d.rectangle([c * tilepx, r * tilepx, (c + 1) * tilepx - 1, (r + 1) * tilepx - 1],
                    outline=(90, 90, 90))
        d.text((c * tilepx + 6, r * tilepx + 6), f"{lat:.4f},{lon:.4f}", fill=(255, 255, 255))
    img.save(out)
    print(f"{out}  ({len(picks)} windows, {span:.0f} m each @ {span/tilepx:.2f} m/px, "
          f"{tag} {suffix})")


if __name__ == '__main__':
    main()
