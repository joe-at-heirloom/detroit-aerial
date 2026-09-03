#!/usr/bin/env python
"""Street-grid bearing of a build versus modern imagery: is it off kilter?

Detroit's street grid is 90-degree symmetric, so a patch of it has a well
defined orientation mod 90. Take the image gradient, weight each pixel's
orientation by its magnitude and sum w*exp(4i*theta): the 4-theta representation
folds the two perpendicular street families onto one phase, so the sum is a
single strong vector whose argument is the grid angle. Do it for the build and
for modern imagery over the same ground, and arg(S_build * conj(S_modern))/4 is
the rotation between the two grids -- no feature matching, no search window, and
it does not care that the buildings changed, only that the streets are still
streets.

|S|/sum(w) is the coherence: how grid-like the patch is. Below ~0.05 there is no
grid to measure (parkland, river, farmland) and the angle means nothing.

Usage:  ./.venv/bin/python scripts/gridbearing.py 1949 placedC [--n 20] [--span 800] [--mpp 1.0]
"""
import sys, os, json, math
import numpy as np, cv2, rasterio
from rasterio.windows import from_bounds
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, 'scripts')); sys.path.insert(0, os.path.join(R, 'pipeline'))
from validate import P, modern_for, reference_for


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def grid_vector(a):
    """sum(w * exp(4i*theta)) and sum(w) for one patch."""
    a = a.astype(np.float32)
    a = a - cv2.GaussianBlur(a, (0, 0), 20)          # drop tone, keep structure
    n = a.shape[0]
    win = np.outer(np.hanning(n), np.hanning(n)).astype(np.float32)   # no frame edges
    a = a * win
    gx = cv2.Sobel(a, cv2.CV_32F, 1, 0, ksize=3); gy = cv2.Sobel(a, cv2.CV_32F, 0, 1, ksize=3)
    w = np.hypot(gx, gy); th = np.arctan2(gy, gx)
    keep = w > np.percentile(w, 70)                  # strongest edges only: roads, not texture
    S = np.sum(w[keep] * np.exp(4j * th[keep])); W = np.sum(w[keep])
    return S, W


def main():
    tag, label = sys.argv[1], sys.argv[2]
    n = int(arg('--n', 20)); span = float(arg('--span', 800)); mpp = float(arg('--mpp', 1.0))
    ref = arg('--ref', 'modern')
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{label}.tif'))
    g = json.load(open(P('data', f'{tag}_{label}_geo.json'))); S_, W_, N_, E_ = g['bbox']
    px = int(span / mpp)
    lats = np.linspace(S_ + 0.015, N_ - 0.015, 8); lons = np.linspace(W_ + 0.015, E_ - 0.015, 5)
    cands = [(la, lo) for la in lats for lo in lons]
    print(f"{tag} {label} vs {ref}: grid bearing over {span:.0f} m patches at {mpp} m/px")
    rows = []
    for la, lo in cands:
        if len(rows) >= n: break
        dlat = span / 111320.0; dlon = span / (111320.0 * math.cos(math.radians(la)))
        s, w, nn, e = la - dlat / 2, lo - dlon / 2, la + dlat / 2, lo + dlon / 2
        a = ds.read(1, window=from_bounds(w, s, e, nn, ds.transform), out_shape=(px, px), boundless=True, fill_value=0)
        if (a > 0).mean() < 0.98: continue
        m = reference_for((s, w, nn, e), px, px, ref) if ref != 'modern' else modern_for((s, w, nn, e), px, px)
        m = np.asarray(m); m = m.mean(axis=2) if m.ndim == 3 else m
        if (m > 0).mean() < 0.98: continue
        Sa, Wa = grid_vector(a); Sb, Wb = grid_vector(m)
        ca, cb = abs(Sa) / (Wa + 1e-9), abs(Sb) / (Wb + 1e-9)
        d = math.degrees(np.angle(Sa * np.conj(Sb))) / 4.0
        rows.append((la, lo, d, ca, cb))
        flag = '' if min(ca, cb) >= 0.05 else '   (weak grid)'
        print(f"  {la:.4f},{lo:.4f}   rotation {d:+7.3f} deg   coherence {ca:.3f}/{cb:.3f}{flag}")
    good = [r for r in rows if min(r[3], r[4]) >= 0.05]
    if good:
        d = np.array([r[2] for r in good])
        mad = 1.4826 * np.median(np.abs(d - np.median(d)))
        print(f"  --> {len(good)} patches with a real grid: rotation median {np.median(d):+.3f} deg, "
              f"MAD {mad:.3f} deg   ({abs(np.median(d))*17.5:.1f} m per km)")
    else:
        print("  --> no patch had a measurable grid")


if __name__ == '__main__':
    main()
