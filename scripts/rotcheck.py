#!/usr/bin/env python
"""Per-window rotation of a build against modern imagery (OpenCV ECC, euclidean).

Answers "is it off kilter": a build can have zero net rotation in its placement
field and still show local rotation if individual frames sit twisted inside the
seam tolerance. Windows of `span` m at `mpp` m/px, spread over the block; each
is registered to modern imagery with a Euclidean (shift + rotation) model.

Usage:  ./.venv/bin/python scripts/rotcheck.py 1949 placedC [--n 12] [--span 1000] [--mpp 1.0]
"""
import sys, os, json, math
import numpy as np, cv2, rasterio
from rasterio.windows import from_bounds
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'pipeline'))
from validate import P, modern_for


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def prep(a):
    a = a.astype(np.float32)
    a = cv2.GaussianBlur(a, (0, 0), 1.5)
    hp = a - cv2.GaussianBlur(a, (0, 0), 25)          # high-pass: roads and edges, not tone
    hp = (hp - hp.mean()) / (hp.std() + 1e-6)
    return hp


def main():
    tag, label = sys.argv[1], sys.argv[2]
    n = int(arg('--n', 12)); span = float(arg('--span', 1000)); mpp = float(arg('--mpp', 1.0))
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{label}.tif'))
    g = json.load(open(P('data', f'{tag}_{label}_geo.json'))); S, W, N, E = g['bbox']
    px = int(span / mpp)
    # a grid of candidate centres, keep those with coverage
    rng = np.random.RandomState(7)
    lats = np.linspace(S + 0.02, N - 0.02, 6); lons = np.linspace(W + 0.02, E - 0.02, 4)
    cands = [(la, lo) for la in lats for lo in lons]; rng.shuffle(cands)
    out = []
    for la, lo in cands:
        if len(out) >= n: break
        dlat = span / 111320.0; dlon = span / (111320.0 * math.cos(math.radians(la)))
        s, w, nn, e = la - dlat / 2, lo - dlon / 2, la + dlat / 2, lo + dlon / 2
        a = ds.read(1, window=from_bounds(w, s, e, nn, ds.transform), out_shape=(px, px), boundless=True, fill_value=0)
        if (a > 0).mean() < 0.95: continue
        m = modern_for((s, w, nn, e), px, px)
        if m.ndim == 3: m = m.mean(axis=2)
        A, B = prep(a), prep(m)
        warp = np.eye(2, 3, dtype=np.float32)
        crit = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 200, 1e-6)
        try:
            cc, warp = cv2.findTransformECC(B, A, warp, cv2.MOTION_EUCLIDEAN, crit, None, 7)
        except cv2.error:
            out.append((la, lo, None, None, None, 0.0)); continue
        th = math.degrees(math.atan2(warp[1, 0], warp[0, 0]))
        # shift at the window centre, in metres
        c = np.array([px / 2, px / 2, 1.0]); p = warp @ c
        dx, dy = (p[0] - px / 2) * mpp, (p[1] - px / 2) * mpp
        out.append((la, lo, th, dx, -dy, float(cc)))
    print(f"{tag} {label}: {len(out)} windows of {span:.0f} m at {mpp} m/px, ECC euclidean vs modern")
    ths = []
    for la, lo, th, dx, dy, cc in out:
        if th is None:
            print(f"  {la:.4f},{lo:.4f}   failed"); continue
        ths.append(th)
        print(f"  {la:.4f},{lo:.4f}   rotation {th:+7.3f} deg   shift dE {dx:+6.1f} dN {dy:+6.1f} m   corr {cc:.2f}")
    if ths:
        ths = np.array(ths)
        print(f"  rotation median {np.median(ths):+.3f} deg, spread (MAD) {1.4826*np.median(np.abs(ths-np.median(ths))):.3f} deg, "
              f"|rot| p90 {np.percentile(np.abs(ths),90):.3f} deg  (0.1 deg = 1.7 m over 1 km)")


if __name__ == '__main__':
    main()
