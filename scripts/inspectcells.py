#!/usr/bin/env python
"""Look at the cells the numbers are arguing about.

Historical goes in the red channel, modern in the green. Streets that coincide
come out yellow; streets that do not split into a red ghost and a green one. Each
cell is rendered twice -- as built, and with the measured displacement applied --
so it is immediately obvious whether the measurement is describing the imagery or
inventing something.

Usage:  ./.venv/bin/python scripts/inspectcells.py 1961 [--n 4] [--suffix rbf]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image, ImageDraw
import rasterio, dtmap
from validate import load, P, MPP
Image.MAX_IMAGE_PIXELS = None

OUT = os.environ.get('INSPECT_OUT', '/tmp')
VIEW = 1.25   # m/px for the rendered crop
SPAN = 1200.0  # metres across each panel


def stretch(a):
    a = a.astype(np.float32); v = a > 0
    if v.sum() < 10:
        return np.zeros_like(a)
    lo, hi = np.percentile(a[v], [2, 98])
    out = np.clip((a - lo) / max(hi - lo, 1) * 255, 0, 255)
    out[~v] = 0
    return out


def main():
    tag = [a for a in sys.argv[1:] if not a.startswith('--')][0]
    n = int(sys.argv[sys.argv.index('--n') + 1]) if '--n' in sys.argv else 4
    suffix = sys.argv[sys.argv.index('--suffix') + 1] if '--suffix' in sys.argv else 'rbf'
    cells = json.load(open(P('data', f'gridval_{tag}.json')))
    worst = sorted([c for c in cells if 'skip' not in c], key=lambda c: -c['mag'])[:n]

    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json'))); bbox = geo['bbox']
    Wp = int((bbox[3] - bbox[1]) * dtmap.MLON / VIEW)
    Hp = int((bbox[2] - bbox[0]) * dtmap.MLAT / VIEW)
    hist = rasterio.open(P('mosaics', f'detroit_{tag}_{suffix}.tif')).read(1, out_shape=(Hp, Wp))
    import validate
    mod = validate.modern_for(bbox, Wp, Hp)
    NY, NX = 16, 3
    ch, cw = Hp // NY, Wp // NX
    half = int(SPAN / 2 / VIEW)
    for c in worst:
        cy = c['row'] * ch + ch // 2; cx = c['col'] * cw + cw // 2
        ys = slice(max(0, cy - half), cy + half); xs = slice(max(0, cx - half), cx + half)
        Hc = stretch(hist[ys, xs]); Mc = stretch(mod[ys, xs])
        panels = []
        for sE, sN in ((0.0, 0.0), (-c['dE'], -c['dN'])):
            px = int(round(sE / VIEW)); py = int(round(-sN / VIEW))
            Hs = np.roll(np.roll(Hc, py, 0), px, 1)
            rgb = np.zeros(Hc.shape + (3,), np.uint8)
            rgb[..., 0] = Hs.astype(np.uint8); rgb[..., 1] = Mc.astype(np.uint8)
            panels.append(rgb)
        img = Image.fromarray(np.concatenate(panels, 1))
        d = ImageDraw.Draw(img)
        d.text((8, 8), f"{tag} cell ({c['row']},{c['col']})  AS BUILT", fill=(255, 255, 255))
        d.text((panels[0].shape[1] + 8, 8),
               f"CORRECTED BY ({-c['dE']:+.1f},{-c['dN']:+.1f}) m  "
               f"[metric says {c['mag']:.0f} m, ratio {c['ratio']:.2f}]", fill=(255, 255, 255))
        p = os.path.join(OUT, f"insp_{tag}_{c['row']}_{c['col']}.png")
        img.save(p)
        print(p, flush=True)


if __name__ == '__main__':
    main()
