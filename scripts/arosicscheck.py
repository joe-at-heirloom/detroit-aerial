#!/usr/bin/env python
"""Third-party verdict: AROSICS against USGS NAIP.

AROSICS is a published, peer-reviewed co-registration package from GFZ Potsdam
(Scheffler et al., Remote Sensing 2017). It shares no code with this project. Its
COREG_LOCAL mode lays a grid of tie points over the overlap of two rasters,
estimates a sub-pixel shift at each by phase correlation in the frequency domain,
and filters the result with its own reliability and outlier tests.

Pointing it at our final mosaic and at NAIP -- imagery from a different
organisation entirely -- gives a number that owes nothing to this project's own
correlation code, its own reference imagery, or its own idea of what "aligned"
means.

Both rasters are written to disk in a projected CRS (UTM 17N) first, because
co-registration in degrees would make a pixel a different size in x and y.

Usage:  ./.venv/bin/python scripts/arosicscheck.py 1961 [--suffix final]
                          [--km 2.4] [--n 3] [--mpp 1.0]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import rasterio
from rasterio.windows import Window
from rasterio.transform import from_origin
import dtmap, validate, naipcheck
from validate import P

OUT = '/tmp/das_arosics'


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def write_ll(path, arr, bbox):
    s, w, n, e = bbox
    H, W = arr.shape
    tr = from_origin(w, n, (e - w) / W, (n - s) / H)
    with rasterio.open(path, 'w', driver='GTiff', height=H, width=W, count=1,
                       dtype='uint8', crs='EPSG:4326', transform=tr) as ds:
        ds.write(arr, 1)


def to_utm(src, dst, res):
    from osgeo import gdal
    gdal.UseExceptions()
    gdal.Warp(dst, src, dstSRS='EPSG:32617', xRes=res, yRes=res,
              resampleAlg='bilinear', srcNodata=0, dstNodata=0)


def main():
    tag = [a for a in sys.argv[1:] if not a.startswith('--')][0]
    suffix = arg('--suffix', 'final')
    km = float(arg('--km', 2.4)); n = int(arg('--n', 3)); mpp = float(arg('--mpp', 1.0))
    # --control plants a known shift on the target before handing it over. If
    # AROSICS does not report that shift back, its verdict on the real data means
    # nothing either.
    plant = float(arg('--control', 0.0))
    os.makedirs(OUT, exist_ok=True)
    from arosics import COREG_LOCAL

    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json'))); bbox = geo['bbox']
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{suffix}.tif'))
    S, W, N, E = bbox
    span = km * 1000.0
    px = int(round(span / mpp))

    def covered(bb):
        c0 = (bb[1] - W) / (E - W) * ds.width; r0 = (N - bb[2]) / (N - S) * ds.height
        wpx = (bb[3] - bb[1]) / (E - W) * ds.width
        hpx = (bb[2] - bb[0]) / (N - S) * ds.height
        a = ds.read(1, window=Window(c0, r0, wpx, hpx), out_shape=(32, 32),
                    boundless=True, fill_value=0)
        return (a > 0).mean() > 0.98

    wins = naipcheck.windows_over(bbox, covered, n, span)
    print(f"AROSICS: {tag} ({suffix}) vs USGS NAIP, {len(wins)} areas of "
          f"{km:.1f} x {km:.1f} km at {mpp} m/px", flush=True)
    if plant:
        print(f"  CONTROL: a shift of ({plant:+.1f}, {-plant/2:+.1f}) m has been planted "
              f"on the target; AROSICS must report it back\n", flush=True)
    else:
        print(flush=True)
    allshift = []
    for k, bb in enumerate(wins):
        nap = naipcheck.fetch(bb, px, px)
        if nap is None or (nap > 0).mean() < 0.95:
            print(f"  area {k+1}: no NAIP"); continue
        c0 = (bb[1] - W) / (E - W) * ds.width; r0 = (N - bb[2]) / (N - S) * ds.height
        wpx = (bb[3] - bb[1]) / (E - W) * ds.width
        hpx = (bb[2] - bb[0]) / (N - S) * ds.height
        hist = ds.read(1, window=Window(c0, r0, wpx, hpx), out_shape=(px, px),
                       boundless=True, fill_value=0)
        if plant:
            import gridval
            hist = gridval.shift_image(hist, plant, -plant / 2, mpp)
        rll = os.path.join(OUT, f'ref_{tag}_{k}.tif')
        tll = os.path.join(OUT, f'tgt_{tag}_{k}.tif')
        write_ll(rll, nap, bb); write_ll(tll, hist, bb)
        ru = rll.replace('.tif', '_utm.tif'); tu = tll.replace('.tif', '_utm.tif')
        to_utm(rll, ru, mpp); to_utm(tll, tu, mpp)
        try:
            CRL = COREG_LOCAL(ru, tu, grid_res=100, window_size=(256, 256),
                              max_shift=40, nodata=(0, 0), min_reliability=30,
                              q=True, progress=False, CPUs=2)
            CRL.calculate_spatial_shifts()
            t = CRL.CoRegPoints_table
            ok = t[t['ABS_SHIFT'] != -9999]
            if len(ok) == 0:
                print(f"  area {k+1}: AROSICS found no reliable tie points"); continue
            a = ok['ABS_SHIFT'].to_numpy(float) * mpp
            xs = ok['X_SHIFT_M'].to_numpy(float); ys = ok['Y_SHIFT_M'].to_numpy(float)
            allshift.append(a)
            c = bb[0] + (bb[2] - bb[0]) / 2, bb[1] + (bb[3] - bb[1]) / 2
            print(f"  area {k+1} @ {c[0]:.4f},{c[1]:.4f}: {len(ok)} tie points  "
                  f"median {np.median(a):5.2f}  p90 {np.percentile(a,90):5.2f}  "
                  f"max {a.max():5.2f} m   bias X {np.median(xs):+5.2f} Y {np.median(ys):+5.2f} m",
                  flush=True)
        except Exception as ex:
            print(f"  area {k+1}: AROSICS failed: {type(ex).__name__}: {ex}", flush=True)
    if allshift:
        a = np.concatenate(allshift)
        print(f"\n  AROSICS overall: {len(a)} tie points  median {np.median(a):5.2f}  "
              f"p90 {np.percentile(a,90):5.2f}  max {a.max():5.2f} m  "
              f">10 m {(a>10).sum()}  >25 m {(a>25).sum()}")


if __name__ == '__main__':
    main()
