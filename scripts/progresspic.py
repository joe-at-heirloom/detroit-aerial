#!/usr/bin/env python
"""Old build versus closed block, same ground, side by side, frame outlines drawn.

The question this answers is the one that started all this: does a road stay
straight across the join between two negatives. Left is what the viewer serves
today; right is the block after the tie-point bundle. Frame outlines are drawn
faintly so a kink can be attributed to a seam.

Usage:  ./.venv/bin/python scripts/progresspic.py 1961 [--span 6000] [--out /tmp]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import rasterio, dtmap
from rasterio.windows import from_bounds
from PIL import Image, ImageDraw
from validate import P
from rebuild import placements


def read_window(tag, label, lat0, lon0, lat1, lon1, out_px):
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{label}.tif'))
    w = from_bounds(lon0, lat0, lon1, lat1, ds.transform)
    a = ds.read(1, window=w, out_shape=(out_px, out_px), boundless=True, fill_value=0)
    return a


def outlines(tag, stage, lat0, lon0, lat1, lon1, out_px, draw, col):
    sol, ppm = placements(tag)
    fp = P('data', f'{stage}_{tag}.json')
    if not os.path.exists(fp):
        return
    saved = json.load(open(fp))
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r].update(dE=v['dE'], dN=v['dN'], rot=v.get('rot', sol[r]['rot']),
                          gw=v.get('gw', sol[r]['gw']), gh=v.get('gh', sol[r]['gh']))
    E0 = (lon0 - dtmap.LON0) * dtmap.MLON; E1 = (lon1 - dtmap.LON0) * dtmap.MLON
    N0 = (lat0 - dtmap.LAT0) * dtmap.MLAT; N1 = (lat1 - dtmap.LAT0) * dtmap.MLAT
    sx = out_px / (E1 - E0); sy = out_px / (N1 - N0)
    for r, p in sol.items():
        if not p.get('ok'): continue
        e = p['e'] - p['dE']; n = p['n'] - p['dN']; th = math.radians(p['rot'])
        c, s = math.cos(th), math.sin(th); gw, gh = p['gw'] * 0.95, p['gh'] * 0.95
        pts = []
        for dx, dy in ((-gw/2, -gh/2), (gw/2, -gh/2), (gw/2, gh/2), (-gw/2, gh/2)):
            X = e + c*dx - s*dy; Y = n + s*dx + c*dy
            pts.append(((X - E0) * sx, (N1 - Y) * sy))
        draw.line(pts + [pts[0]], fill=col, width=1)


def stretch(a):
    a = a.astype(np.float32); v = a > 0
    if v.sum() < 10: return np.zeros_like(a, np.uint8)
    lo, hi = np.percentile(a[v], [1, 99]); o = np.clip((a - lo) / max(hi - lo, 1) * 255, 0, 255); o[~v] = 0
    return o.astype(np.uint8)


def main():
    tag = sys.argv[1]
    span = float(sys.argv[sys.argv.index('--span') + 1]) if '--span' in sys.argv else 6000.0
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else '/tmp'
    newlab, newstage = 'closedA3', 'stageA3'
    oldlab, oldstage = 'final', 'frameadj'
    g = json.load(open(P('data', f'{tag}_{newlab}_geo.json'))); S, W, N, E = g['bbox']
    dlat = span / dtmap.MLAT; dlon = span / dtmap.MLON
    # two windows: one third and two thirds down the block. Centred across the
    # block by default; with --sidelap, centred on the join between the two
    # flight lines that overlap most, where a cross-line seam would show.
    panels = []; px = 900
    lon_c = (W + E) / 2
    if '--sidelap' in sys.argv:
        import seamclass as SC
        sol, ppm = placements(tag)
        recs = sorted(r for r in sol if sol[r].get('ok')); lab = SC.lines_of(sol, recs)
        cE = {}
        for r in recs:
            cE.setdefault(lab[r], []).append(sol[r]['e'])
        ks = sorted(cE); mids = [(np.mean(cE[ks[i]]) + np.mean(cE[ks[i + 1]])) / 2 for i in range(len(ks) - 1)]
        lon_c = dtmap.LON0 + mids[len(mids) // 2] / dtmap.MLON
    for f in (0.33, 0.62):
        lat = N - f * (N - S); lon = lon_c
        lat0, lat1, lon0, lon1 = lat - dlat/2, lat + dlat/2, lon - dlon/2, lon + dlon/2
        row = []
        for lab, stage, title in ((oldlab, oldstage, 'CURRENT VIEWER BUILD'), (newlab, newstage, 'CLOSED BLOCK (tie-point bundle)')):
            try:
                a = read_window(tag, lab, lat0, lon0, lat1, lon1, px)
            except Exception as ex:
                a = np.zeros((px, px), np.uint8)
            im = Image.fromarray(stretch(a)).convert('RGB'); d = ImageDraw.Draw(im)
            outlines(tag, stage, lat0, lon0, lat1, lon1, px, d, (255, 120, 0))
            d.rectangle((0, 0, px, 20), fill=(0, 0, 0)); d.text((6, 4), f"{tag}  {title}   {span/1000:.0f} km across", fill=(255, 255, 255))
            row.append(im)
        panels.append(row)
    sheet = Image.new('RGB', (2 * px + 10, 2 * px + 10), (0, 0, 0))
    for j, row in enumerate(panels):
        for i, im in enumerate(row):
            sheet.paste(im, (i * (px + 10), j * (px + 10)))
    p = os.path.join(out, f"progress_{tag}{'_sidelap' if '--sidelap' in sys.argv else ''}_{int(span)}.png"); sheet.save(p); print(p)


if __name__ == '__main__':
    main()
