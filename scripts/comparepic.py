#!/usr/bin/env python
"""Same ground, three methods, side by side, frame outlines drawn.

  old      = per-frame correction + RBF warp (what the viewer served for a day)
  bundle   = hand-rolled per-frame tie-point bundle (close3), placed
  colmap   = COLMAP cameras projected onto the ground plane (pilot: 6 frames)

Usage:  ./.venv/bin/python scripts/comparepic.py 1961 [--span 2500] [--out /tmp]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import rasterio, dtmap
from rasterio.windows import from_bounds
from PIL import Image, ImageDraw
from validate import P
import progresspic as PP
import colmap_render as CR


def colmap_footprints(work, model_dir=None, crop=0.95):
    meta, cams, imgs, zg = CR.load_block(work, model_dir)
    E0, N0 = meta['E0'], meta['N0']; out = {}
    for n, img in imgs.items():
        cam = cams[img['cam']]; w, h = cam['w'], cam['h']; f = cam['params'][0]; cx, cy = cam['params'][1], cam['params'][2]
        pts = []
        for (u, v) in ((w*(1-crop)/2, h*(1-crop)/2), (w*(1+crop)/2, h*(1-crop)/2), (w*(1+crop)/2, h*(1+crop)/2), (w*(1-crop)/2, h*(1+crop)/2)):
            d = img['R'].T @ np.array([(u - cx)/f, (v - cy)/f, 1.0]); C = img['C']; s = (zg - C[2]) / d[2]; X = C + s*d
            pts.append((X[0] + E0, X[1] + N0))
        out[n] = pts
    return out


def main():
    tag = sys.argv[1]
    span = float(sys.argv[sys.argv.index('--span') + 1]) if '--span' in sys.argv else 2500.0
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else '/tmp'
    # --colmap-work DIR : a full-block COLMAP solve (default: the 6-frame pilot)
    work = sys.argv[sys.argv.index('--colmap-work') + 1] if '--colmap-work' in sys.argv else '/tmp/colmap_pilot'
    model = sys.argv[sys.argv.index('--colmap-model') + 1] if '--colmap-model' in sys.argv else (None if work != '/tmp/colmap_pilot' else '/tmp/colmap_pilot/sparse_pp_txt')
    ctag = 'pilot' if work == '/tmp/colmap_pilot' else tag
    clabel = 'colmap' if not os.path.exists(P('data', f'{ctag}_placedC_geo.json')) else 'placedC'
    cgeo = json.load(open(P('data', f'{ctag}_{clabel}_geo.json'))); S, W, N, E = cgeo['bbox']
    foot = colmap_footprints(work, model)
    # windows centred on the join between the two flight lines with the most
    # sidelap when there are several lines, else across the strip
    cEs = np.array([np.mean([p[0] for p in v]) for v in foot.values()])
    lines = np.sort(cEs); gaps = np.where(np.diff(lines) > 800)[0]
    if len(gaps):
        groups = np.split(lines, gaps + 1); mids = [(groups[i].mean() + groups[i+1].mean()) / 2 for i in range(len(groups) - 1)]
        cE = mids[len(mids) // 2]
    else:
        cE = cEs.mean()
    lon_c = dtmap.LON0 + cE / dtmap.MLON
    dlat = span / dtmap.MLAT; dlon = span / dtmap.MLON
    px = 800
    sources = (('OLD (per-frame + RBF)', f'{tag}', 'final', 'frameadj'),
               ('BUNDLE (close3, placed)', f'{tag}', 'placed3', 'stageA3'),
               (f'COLMAP ({clabel})', ctag, clabel, None))
    rows = []
    for f in (0.36, 0.64):
        lat = N - f * (N - S); lat0, lat1, lon0, lon1 = lat - dlat/2, lat + dlat/2, lon_c - dlon/2, lon_c + dlon/2
        row = []
        for title, t, lab, stage in sources:
            try:
                a = PP.read_window(t, lab, lat0, lon0, lat1, lon1, px)
            except Exception:
                a = np.zeros((px, px), np.uint8)
            im = Image.fromarray(PP.stretch(a)).convert('RGB'); d = ImageDraw.Draw(im)
            if stage:
                PP.outlines(t, stage, lat0, lon0, lat1, lon1, px, d, (255, 120, 0))
            else:
                E0 = (lon0 - dtmap.LON0) * dtmap.MLON; E1 = (lon1 - dtmap.LON0) * dtmap.MLON
                N0 = (lat0 - dtmap.LAT0) * dtmap.MLAT; N1 = (lat1 - dtmap.LAT0) * dtmap.MLAT
                sx = px / (E1 - E0); sy = px / (N1 - N0)
                for pts in foot.values():
                    q = [((x - E0) * sx, (N1 - y) * sy) for x, y in pts]
                    d.line(q + [q[0]], fill=(255, 120, 0), width=1)
            d.rectangle((0, 0, px, 20), fill=(0, 0, 0)); d.text((6, 4), f"{tag}  {title}   {span/1000:.1f} km", fill=(255, 255, 255))
            row.append(im)
        rows.append(row)
    sheet = Image.new('RGB', (3 * px + 20, 2 * px + 10), (0, 0, 0))
    for j, row in enumerate(rows):
        for i, im in enumerate(row):
            sheet.paste(im, (i * (px + 10), j * (px + 10)))
    p = os.path.join(out, f'compare_{tag}_{ctag}_{int(span)}.png'); sheet.save(p); print(p)


if __name__ == '__main__':
    main()
