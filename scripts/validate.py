#!/usr/bin/env python
"""Measure a block against modern imagery -- and prove the ruler first.

Three things run here, in this order, because the project has twice been fooled
by a number that turned out to be the metric's own failure:

  1. planted-shift calibration on modern-vs-modern, masked to the historical
     footprint. Whatever this reports is the metric's error, not the mosaic's.
  2. self-consistency: plant a shift on the real historical mosaic. The measured
     displacement must move by exactly that much. A metric that locks onto the
     wrong street lattice fails this even when step 1 looks perfect.
  3. the actual measurement, on a 16 x 3 grid, reported as median / p90 / max /
     count over 25 m -- never a single average.

Usage:  ./.venv/bin/python scripts/validate.py 1961 1967 [--suffix final]
"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'pipeline'))
import gridval, dtmap
from PIL import Image
import rasterio
Image.MAX_IMAGE_PIXELS = None

MPP = 2.5
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def P(*a): return os.path.join(ROOT, *a)

# Prefer the high-resolution reference when it is present. The original
# modern_west.png is 3.5 m/px and stops short of the 1949 and 1956 blocks; the
# fetched one is 1.8 m/px, covers all four, and is resampled to this project's
# linear-latitude grid rather than left in Web Mercator.
def _ref():
    hi = P('data', 'modern_west_hi_geo.json')
    if os.path.exists(hi) and os.path.exists(P('mosaics', 'modern_west_hi.tif')):
        return json.load(open(hi))['bbox'], P('mosaics', 'modern_west_hi.tif'), True
    return (json.load(open(P('data', 'modern_west_geo.json')))['bbox'],
            P('mosaics', 'modern_west.png'), False)

_MODB, _MODPATH, _MODHI = _ref()
_MOD = None
def modern_for(bbox, W, H):
    """Resample the modern reference onto a block's exact grid.

    This used to crop to the overlap and then resize the crop to the full block
    size. Where a block extends past the modern imagery -- 1949 runs about 3 km
    further south, 1956 about 3 km further south and west -- that silently
    STRETCHED the reference across the whole block, so every comparison against it
    was measuring against a distorted map. Sample properly instead and leave the
    uncovered margin as nodata, which the masked correlation already knows how to
    ignore."""
    global _MOD
    if _MOD is None:
        if _MODHI:
            _MOD = rasterio.open(_MODPATH).read(1)
        else:
            _MOD = np.asarray(Image.open(_MODPATH).convert('L'))
    ms, mw, mn, me = _MODB; mh, mwd = _MOD.shape
    s, w, n, e = bbox
    lat = n - (np.arange(H) + 0.5) * (n - s) / H
    lon = w + (np.arange(W) + 0.5) * (e - w) / W
    sy = (mn - lat) / (mn - ms) * mh - 0.5
    sx = (lon - mw) / (me - mw) * mwd - 0.5
    oky = (sy >= 0) & (sy <= mh - 1)
    okx = (sx >= 0) & (sx <= mwd - 1)
    yc = np.clip(sy, 0, mh - 1.001); xc = np.clip(sx, 0, mwd - 1.001)
    y0 = yc.astype(np.int32); x0 = xc.astype(np.int32)
    fy = (yc - y0).astype(np.float32)[:, None]; fx = (xc - x0).astype(np.float32)[None, :]
    A = _MOD[y0][:, x0].astype(np.float32); B = _MOD[y0][:, x0 + 1].astype(np.float32)
    C = _MOD[y0 + 1][:, x0].astype(np.float32); D = _MOD[y0 + 1][:, x0 + 1].astype(np.float32)
    out = (A * (1 - fx) * (1 - fy) + B * fx * (1 - fy) + C * (1 - fx) * fy + D * fx * fy)
    out[~oky, :] = 0; out[:, ~okx] = 0
    return np.clip(out, 0, 255).astype(np.uint8)


def load(tag, suffix='rbf'):
    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json')))
    bbox = geo['bbox']
    Wp = int((bbox[3] - bbox[1]) * dtmap.MLON / MPP)
    Hp = int((bbox[2] - bbox[0]) * dtmap.MLAT / MPP)
    arr = rasterio.open(P('mosaics', f'detroit_{tag}_{suffix}.tif')).read(1, out_shape=(Hp, Wp))
    return arr, modern_for(bbox, Wp, Hp), bbox


def report(recs, label, lock=1.15):
    """Report the measurement conditioned on whether the measurement is valid.

    A cell whose correlation has no unique peak does not have a large error -- it
    has an unknown one, and quoting its number as error is how this project got
    "231.9 m worst cell" out of a mosaic that was roughly right. Cell (5,1) of 1961
    is the worked example: it reported 39.9 m at ratio 1.07, and applying that
    correction visibly splits an arterial that was already aligned.

    So: headline over cells that locked, and an explicit count of the ones that did
    not, which are neither pass nor fail but unmeasured."""
    ok = [r for r in recs if 'skip' not in r and not r.get('pegged')]
    peg = [r for r in recs if 'skip' not in r and r.get('pegged')]
    skip = [r for r in recs if 'skip' in r]
    if not ok:
        print(f"  {label}: nothing measurable"); return None
    good = [r for r in ok if r['ratio'] >= lock]
    weak = [r for r in ok if r['ratio'] < lock]
    if not good:
        print(f"  {label}: no cell locked (all {len(ok)} below ratio {lock})"); return None
    m = np.array([r['mag'] for r in good])
    print(f"  {label}: n={len(good)}  median {np.median(m):5.1f}  p90 {np.percentile(m,90):5.1f}  "
          f"max {m.max():5.1f}  >10m {(m>10).sum()}  >25m {(m>25).sum()}"
          f"   | unlocked {len(weak)}  pegged {len(peg)}  no-coverage {len(skip)}")
    strict = np.array([r['mag'] for r in ok if r['ratio'] >= 1.25])
    if len(strict):
        print(f"      confident subset (ratio>=1.25, n={len(strict)}): median "
              f"{np.median(strict):.1f}  p90 {np.percentile(strict,90):.1f}  "
              f"max {strict.max():.1f}  >25m {(strict>25).sum()}")
    return m


def main():
    tags = [a for a in sys.argv[1:] if not a.startswith('--')]
    suffix = (sys.argv[sys.argv.index('--suffix') + 1]
              if '--suffix' in sys.argv else 'final')
    tags = [t for t in tags if t != suffix]
    for tag in tags:
        t0 = time.time()
        arr, mod, bbox = load(tag, suffix)
        print(f"\n===== {tag} ({suffix}) =====  {arr.shape[1]}x{arr.shape[0]} @ {MPP} m/px  "
              f"coverage {(arr>0).mean():.2f}", flush=True)
        t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_{suffix}')
        print(f"  ridge maps ready [{time.time()-t0:.0f}s]", flush=True)

        # Calibration reuses the modern ridge maps as the "historical" side, masked
        # to the real footprint: identical content, so anything it reports is the
        # metric's own error.
        tc = dict(t)
        tc['hc'] = t['mc'] * t['hvc']; tc['hf'] = t['mf'] * t['hv']
        prior_c = gridval.Prior(gridval.coarse_field(tc, log=lambda *_: None), MPP)

        print(" 1. CALIBRATION -- modern vs modern, masked to the historical footprint")
        for dE, dN in ((0.0, 0.0), (37.5, -22.5), (120.0, 87.5)):
            recs = gridval.grid_hier_prep(tc, prior=prior_c, hshift=(dE, dN))
            u = [r for r in recs if 'skip' not in r and not r.get('pegged')]
            err = np.array([math.hypot(r['dE'] - dE, r['dN'] - dN) for r in u]) if u else np.array([9e9])
            print(f"    planted ({dE:+6.1f},{dN:+6.1f})  recovery error: median {np.median(err):4.1f}  "
                  f"max {err.max():5.1f}  wrong by >10 m: {(err>10).sum()}/{len(u)}", flush=True)

        prior = gridval.Prior(gridval.coarse_field(t), MPP)
        print(" 2. SELF-CONSISTENCY -- plant a shift on the real historical mosaic")
        base = gridval.grid_hier_prep(t, prior=prior)
        B = {(r['row'], r['col']): r for r in base if 'skip' not in r and not r.get('pegged')}
        for dE, dN in ((25.0, 0.0), (0.0, -60.0), (-75.0, 45.0)):
            recs = gridval.grid_hier_prep(t, prior=prior, hshift=(dE, dN))
            d = np.array([math.hypot(r['dE'] - B[(r['row'], r['col'])]['dE'] - dE,
                                     r['dN'] - B[(r['row'], r['col'])]['dN'] - dN)
                          for r in recs if 'skip' not in r and not r.get('pegged')
                          and (r['row'], r['col']) in B] or [9e9])
            print(f"    planted ({dE:+6.1f},{dN:+6.1f})  tracked to within: median {np.median(d):4.1f}  "
                  f"max {d.max():5.1f}  inconsistent >10 m: {(d>10).sum()}/{len(d)}", flush=True)

        print(" 3. MEASUREMENT")
        report(base, f'{tag}')
        json.dump(base, open(P('data', f'gridval_{tag}.json'), 'w'))
        print("    worst cells:  row col valid    dE     dN    mag  ratio pegged")
        for r in sorted([x for x in base if 'skip' not in x], key=lambda r: -r['mag'])[:10]:
            print(f"                  {r['row']:3d} {r['col']:3d} {r['valid']:.2f} "
                  f"{r['dE']:+6.1f} {r['dN']:+6.1f} {r['mag']:6.1f} {r['ratio']:5.2f}  {r['pegged']}")
        print(f"    [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
