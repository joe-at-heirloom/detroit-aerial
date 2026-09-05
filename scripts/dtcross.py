#!/usr/bin/env python
"""How far do the downtown epochs disagree with EACH OTHER?

Every downtown number so far has been absolute -- each epoch measured against a
modern reference on its own. That is the right question for placing a mosaic and
the wrong one for a wipe: what you see when you drag the divider is one epoch
against another, and two epochs can each be 15 m from the truth in opposite
directions and be 30 m apart on screen.

So measure the thing the viewer actually shows. Both epochs are resampled onto
one grid, then correlated cell by cell -- the same instrument the block seams are
judged with, pointed sideways in time instead of across a flight line.

    ./.venv/bin/python scripts/dtcross.py [--mpp 1.25] [--n 5]
"""
import sys, os, json, math, argparse
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image
import rasterio
import gridval, dtmap, handadjust
Image.MAX_IMAGE_PIXELS = None

P = lambda *a: os.path.join(ROOT, *a)
TAGS = ['1949', '1956', '1961', 'now']


def layers():
    man = json.load(open(P('data', 'manifest.json')))
    return {l['id'][2:]: l for l in man['layers'] if l['id'].startswith('dt')}


def onto(l, bbox, W, H, adj=None):
    """Resample one layer onto a shared lat/lon grid. Nearest-neighbour: this is a
    measurement, and interpolation would smooth the very edges being measured.

    `adj` is the layer's hand alignment, so this measures what the viewer is
    actually serving rather than what is on disk underneath it."""
    S, W_, N, E = l['bbox']
    p = P(l['file'])
    if p.lower().endswith(('.tif', '.tiff')):
        ds = rasterio.open(p)
        im = ds.read(1)
    else:
        im = np.asarray(Image.open(p).convert('L'))
    h, w = im.shape
    s, wl, n, e = bbox
    lat = n - (np.arange(H) + 0.5) * (n - s) / H
    lon = wl + (np.arange(W) + 0.5) * (e - wl) / W
    if adj is None:
        # separable, and much the cheaper path
        sy = (N - lat) / (N - S) * h - 0.5
        sx = (lon - W_) / (E - W_) * w - 0.5
        oy = (sy >= 0) & (sy < h); ox = (sx >= 0) & (sx < w)
        out = np.zeros((H, W), np.uint8)
        if oy.any() and ox.any():
            out[np.ix_(oy, ox)] = im[np.clip(sy[oy], 0, h - 1).astype(np.int32)][
                                     :, np.clip(sx[ox], 0, w - 1).astype(np.int32)]
        return out
    LON, LAT = np.meshgrid(lon, lat)
    lon_s, lat_s = handadjust.source_lonlat(LON, LAT, (S, W_, N, E), adj)
    sx = (lon_s - W_) / (E - W_) * w - 0.5
    sy = (N - lat_s) / (N - S) * h - 0.5
    ok = (sx >= 0) & (sx < w) & (sy >= 0) & (sy < h)
    out = np.zeros((H, W), np.uint8)
    out[ok] = im[np.clip(sy[ok], 0, h - 1).astype(np.int32),
                 np.clip(sx[ok], 0, w - 1).astype(np.int32)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mpp', type=float, default=1.25)
    ap.add_argument('--n', type=int, default=5, help='grid is n x n cells')
    a = ap.parse_args()

    L = layers()
    have = [t for t in TAGS if t in L]
    # the ground every epoch actually covers, so no pair is judged on a cell one
    # of them never flew
    bb = [max(L[t]['bbox'][0] for t in have), max(L[t]['bbox'][1] for t in have),
          min(L[t]['bbox'][2] for t in have), min(L[t]['bbox'][3] for t in have)]
    W = int((bb[3] - bb[1]) * dtmap.MLON / a.mpp)
    H = int((bb[2] - bb[0]) * dtmap.MLAT / a.mpp)
    print(f"downtown overlap {bb[3]-bb[1]:.4f} x {bb[2]-bb[0]:.4f} deg "
          f"-> {W}x{H} @ {a.mpp} m/px\n")
    ADJ = handadjust.load(ROOT)
    G = {t: onto(L[t], bb, W, H, handadjust.for_layer(ADJ, 'dt' + t)) for t in have}
    for t in have:
        a = handadjust.for_layer(ADJ, 'dt' + t)
        note = ('   hand: %+.1f,%+.1f m  %+.2f deg  x%.4f'
                % (a['dE'], a['dN'], a['deg'], a['scale'])) if a else ''
        print(f"  {t:>4}  coverage {(G[t]>0).mean():.2f}{note}")
    print()

    rows = []
    for i, ta in enumerate(have):
        for tb in have[i + 1:]:
            t = gridval.prepare_cached(G[ta], G[tb], a.mpp, f'dtx_{ta}_{tb}_{a.mpp}')
            recs = gridval.grid_hier_prep(t, NY=a.n, NX=a.n,
                                          prior=gridval.Prior([], a.mpp), min_valid=0.35)
            # `ratio` is the peak's margin over the runner-up; below ~1.15 the
            # correlation has not chosen, and a pegged cell has run to its search
            # edge, which is a failure reporting itself as a large displacement.
            d = np.array([math.hypot(r['dE'], r['dN']) for r in recs
                          if 'skip' not in r and not r.get('pegged')
                          and r.get('ratio', 0) >= 1.15])
            tot = sum(1 for r in recs if 'skip' not in r)
            if not len(d):
                print(f"  {ta} vs {tb}: no cell locked ({tot} tried)"); continue
            rows.append((ta, tb, len(d), np.median(d), np.percentile(d, 90), d.max()))
            print(f"  {ta} vs {tb}:  {len(d):2d}/{tot:2d} locked   median {np.median(d):5.1f}  "
                  f"p90 {np.percentile(d,90):5.1f}  max {d.max():5.1f} m", flush=True)

    print("\n  the wipe shows a PAIR, so the pair is the number that matters")
    if rows:
        worst = max(rows, key=lambda r: r[3])
        print(f"  worst pair: {worst[0]} vs {worst[1]} at {worst[3]:.1f} m median")


if __name__ == '__main__':
    main()
