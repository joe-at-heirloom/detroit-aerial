#!/usr/bin/env python
"""Build a candidate raster from a jointfit.py correction (plan step 5).

The correction is composed with the block's existing placement and applied to the
COLMAP mosaic in ONE resampling -- out(p) = colmap(p + d_new(p) + d_old(p + d_new(p)))
-- never as a second warp of the already-warped served raster. Written under a new
label; the manifest is not touched.

Before any pixel is written the field is checked the way the plan asks:

  size       |d_new| over the whole raster, and separately over the ground the
             measurements covered: a cubic is free to swing where nothing held it
  gradient   the largest local scale / rotation of the new field; a continuous
             warp moves the two sides of a seam by (gradient x their separation,
             ~2 m), so this bounds the seam change directly
  foldover   det(I + grad d) must stay positive everywhere

Outside the box the measurements span, coordinates are clamped: the correction
continues flat rather than extrapolating a polynomial nobody measured.

Usage:  ./.venv/bin/python scripts/jointapply.py 1956 [--fit runs/jointfit_west.json]
            [--label placedJ] [--taper 1500] [--check-only]
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import dtmap
import jointfit as JF
import fieldapply as FA
import rbfapply
P = lambda *a: os.path.join(ROOT, *a)


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def main():
    tag = sys.argv[1]; lid = 'b' + tag
    fit = json.load(open(arg('--fit', P('runs', 'jointfit_west.json'))))
    label = arg('--label', 'placedJ'); taper = float(arg('--taper', 1500.0))
    order = fit['best'][lid]
    k = JF.basis(order, 0.0, 0.0).shape[1]
    p = np.array(fit['params'][lid]); cen = np.array(fit['centres'][lid])
    served = json.load(open(P('data', 'manifest.json')))['serve'][tag]
    old = json.load(open(P('data', f'fieldmodel_{tag}_{served}.json')))
    g = json.load(open(P('data', f'{tag}_colmap_geo.json')))
    print(f"{tag}: correction {order} ({k} params) composed with the served {served} placement "
          f"({'+'.join(m['kind'] for m in old)})", flush=True)

    # where the correction was measured: this layer's fitted cells (vs ref or film)
    run = json.load(open(P(fit['run'])))
    obs, _ = JF.cells(run, fit['ref'], fit['check'], 1.15, fit['gate'])
    pts = np.array([(o['E'], o['N']) for o in obs if lid in (o['a'], o['b'])
                    and fit['check'] not in (o['a'], o['b'])])
    from scipy.spatial import cKDTree
    tree = cKDTree(pts)

    def raw(E, N):
        return JF.basis(order, (E - cen[0]) / JF.SC, (N - cen[1]) / JF.SC) @ p

    lo = pts.min(0) - 650.0; hi = pts.max(0) + 650.0      # half a cell past the last measurement

    def d_new(E, N):
        # Beyond the measured ground the polynomial is not held by anything, so it
        # is not evaluated there: coordinates are clamped to the measurements' box
        # and the field continues flat. Continuous by construction -- blending
        # toward the NEAREST measurement jumps from cell to cell and put a 2.3%
        # gradient into the field.
        return raw(min(max(E, lo[0]), hi[0]), min(max(N, lo[1]), hi[1]))

    # --- checks on a 250 m lattice over the raster
    step = 250.0
    Es = g['minE'] + np.arange(0, g['W'] * g['mpp'], step)
    Ns = g['maxN'] - np.arange(0, g['H'] * g['mpp'], step)
    D = np.array([[d_new(E, N) for E in Es] for N in Ns])        # (ny, nx, 2)
    near = np.array([[tree.query([E, N])[0] <= 1300.0 for E in Es] for N in Ns])
    mag = np.hypot(D[..., 0], D[..., 1])
    # gradient in metres per metre
    dEdx = np.gradient(D[..., 0], step, axis=1); dEdy = -np.gradient(D[..., 0], step, axis=0)
    dNdx = np.gradient(D[..., 1], step, axis=1); dNdy = -np.gradient(D[..., 1], step, axis=0)
    grad = np.sqrt(dEdx ** 2 + dEdy ** 2 + dNdx ** 2 + dNdy ** 2)
    det = (1 + dEdx) * (1 + dNdy) - dEdy * dNdx
    rawmag = np.array([[np.hypot(*raw(E, N)) for E in Es] for N in Ns])
    print(f"  |d_new| over measured ground: median {np.median(mag[near]):.1f}  max {mag[near].max():.1f} m")
    print(f"  |d_new| beyond it (clamped):   max {mag[~near].max() if (~near).any() else 0:.1f} m "
          f"(unclamped it would reach {rawmag[~near].max() if (~near).any() else 0:.1f} m)")
    print(f"  largest local gradient {grad.max()*1000:.2f} m/km -> seam change under "
          f"{grad.max()*2.0*100:.1f} cm for a 2 m seam; min det {det.min():.5f} (foldover if <= 0)")
    ok = det.min() > 0.9 and grad.max() < 0.01
    if not ok:
        raise SystemExit("  field fails the gradient/foldover check; not applying")
    if '--check-only' in sys.argv:
        return

    old_at = FA.make_warp(old)

    def warp_at(E, N):
        a_ = d_new(E, N)
        b_ = old_at(E + a_[0], N + a_[1])
        return (a_[0] + b_[0], a_[1] + b_[1], 0.0)

    t0 = time.time()
    out = P('mosaics', f'detroit_{tag}_{label}.tif')
    tmp = os.path.join(os.environ.get('TMPDIR', '/tmp'), f'ja_{tag}.raw')
    res = rbfapply.apply(P('mosaics', f'detroit_{tag}_colmap.tif'), dict(
        minE=g['minE'], maxN=g['maxN'], W=g['W'], H=g['H'], mpp=g['mpp']),
        warp_at, out, tmp, step_m=250.0, log=lambda s: None)
    json.dump(res, open(P('data', f'{tag}_{label}_geo.json'), 'w'))
    json.dump(dict(order=order, params=p.tolist(), centre=cen.tolist(), taper=taper,
                   composed_with=f'{served}', fit=os.path.relpath(arg('--fit', P('runs', 'jointfit_west.json')), ROOT)),
              open(P('data', f'jointmodel_{tag}_{label}.json'), 'w'))
    print(f"  wrote {os.path.relpath(out, ROOT)} [{time.time()-t0:.0f}s] -- not served", flush=True)


if __name__ == '__main__':
    main()
