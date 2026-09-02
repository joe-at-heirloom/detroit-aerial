#!/usr/bin/env python
"""Coarse same-film alignment: where does a COLMAP block sit relative to the old
build of the same year, when the answer may be hundreds of metres and a few degrees?

A COLMAP block is self-aligned to the catalogue positions, which for 1949 scatter
361 m around the solved cameras; the heading fitted to them was a few degrees off,
which is ~900 m at the block's ends -- far outside any fine search. The old build
of the same year is the same film, so big windows at 20 m/px phase-correlate to a
clean delta even with a few degrees of rotation, and a similarity fitted to those
offsets recovers the heading, scale and shift. Written in fieldfit's cell format so
fieldfit/fieldapply do the rest.

Usage:  ./.venv/bin/python scripts/coarsealign.py 1949 colmap final [--mpp 20] [--win 4000] [--out coarse]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import rasterio, dtmap, close2
from validate import P, reference_for


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def main():
    tag, src, ref = sys.argv[1], sys.argv[2], sys.argv[3]
    mpp = float(arg('--mpp', 20.0)); win_m = float(arg('--win', 4000.0)); out = arg('--out', 'coarse')
    g = json.load(open(P('data', f'{tag}_{src}_geo.json'))); bbox = g['bbox']
    W = int((bbox[3] - bbox[1]) * dtmap.MLON / mpp); H = int((bbox[2] - bbox[0]) * dtmap.MLAT / mpp)
    A = rasterio.open(P('mosaics', f'detroit_{tag}_{src}.tif')).read(1, out_shape=(H, W))
    B = reference_for(bbox, W, H, f'{tag}:{ref}')
    minE = (bbox[1] - dtmap.LON0) * dtmap.MLON; maxN = (bbox[2] - dtmap.LAT0) * dtmap.MLAT
    win = int(win_m / mpp); step = win // 2
    cells = []
    for y0 in range(0, H - win + 1, step):
        for x0 in range(0, W - win + 1, step):
            a = A[y0:y0 + win, x0:x0 + win]; b = B[y0:y0 + win, x0:x0 + win]
            if (a > 0).mean() < 0.7 or (b > 0).mean() < 0.7:
                continue
            dx, dy, sharp, edge = close2._pc(a.astype(np.float32), b.astype(np.float32), pad=2)
            if edge or sharp < 6.0:
                continue
            cells.append(dict(row=y0 // step, col=x0 // step, cE=minE + (x0 + win / 2) * mpp,
                              cN=maxN - (y0 + win / 2) * mpp, dE=dx * mpp, dN=-dy * mpp,
                              mag=math.hypot(dx, dy) * mpp, ratio=1.0 + min(sharp / 20.0, 2.0),
                              pegged=False, valid=1.0))
    if not cells:
        raise SystemExit("no window locked")
    m = np.array([c['mag'] for c in cells])
    print(f"{tag} {src} vs {ref}: {len(cells)} windows of {win_m:.0f} m at {mpp:.0f} m/px locked; "
          f"offset median {np.median(m):.0f}  p90 {np.percentile(m,90):.0f}  max {m.max():.0f} m", flush=True)
    json.dump(dict(tag=tag, label=src, bbox=bbox, mpp=mpp, minE=minE, maxN=maxN, shape=[H, W],
                   cells_16x3=cells, cells_32x6=cells, ref=f'{tag}:{ref}'),
              open(P('data', f'field_{tag}_{out}.json'), 'w'))
    print(f"saved data/field_{tag}_{out}.json")


if __name__ == '__main__':
    main()
