#!/usr/bin/env python
"""Put the downtown epochs through the same correction and the same scrutiny.

Downtown was carried over from an earlier stage of the project and never checked
against anything independent. Its claimed 3-5 m came from validation against road
centreline vectors -- a metric this project has already established has a 5-6 m
noise floor, because modern imagery scores the same against those vectors as our
imagery does. Measured against USGS NAIP instead, the downtown epochs sit 14-26 m
out.

Two distinct things are mixed into that number and only one is fixable here:

  the ground   streets and intersections genuinely offset -- a residual field
               removes this, exactly as for the blocks
  the roofs    downtown Detroit is tall, the negatives are not orthorectified, and
               relief displacement leans every building differently in every frame.
               No 2D warp can fix that. It needs a DEM.

So this corrects the ground and reports the roofs as the known limit.

Usage:  ./.venv/bin/python scripts/fixdowntown.py [1949 1956 1961] [--dry]
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image
import rasterio
from rasterio.transform import from_origin
import gridval, warpsolve, rbfapply, dtmap
import validate
from validate import P, report
Image.MAX_IMAGE_PIXELS = None
MPP = 1.25          # downtown is small; measure it finer than the blocks


def layer(tag):
    man = json.load(open(P('data', 'manifest.json')))
    for l in man['layers']:
        if l['id'] == f'dt{tag}':
            return l
    return None


def as_grid(l, mpp=MPP):
    """Sample a downtown PNG onto a linear lat/lon grid, and give it a geo dict in
    the same form the block pipeline uses."""
    S, W, N, E = l['bbox']
    im = np.asarray(Image.open(P(l['file'])).convert('L'))
    H, Wd = im.shape
    Wp = int((E - W) * dtmap.MLON / mpp); Hp = int((N - S) * dtmap.MLAT / mpp)
    lat = N - (np.arange(Hp) + 0.5) * (N - S) / Hp
    lon = W + (np.arange(Wp) + 0.5) * (E - W) / Wp
    sy = np.clip((N - lat) / (N - S) * H - 0.5, 0, H - 1).astype(np.int32)
    sx = np.clip((lon - W) / (E - W) * Wd - 0.5, 0, Wd - 1).astype(np.int32)
    arr = im[sy][:, sx]
    geo = dict(minE=(W - dtmap.LON0) * dtmap.MLON,
               maxN=(N - dtmap.LAT0) * dtmap.MLAT,
               W=Wp, H=Hp, mpp=mpp, bbox=[S, W, N, E])
    return arr, geo


def write_tif(path, arr, geo):
    S, W, N, E = geo['bbox']
    tr = from_origin(W, N, (E - W) / geo['W'], (N - S) / geo['H'])
    with rasterio.open(path, 'w', driver='GTiff', height=geo['H'], width=geo['W'],
                       count=1, dtype='uint8', crs='EPSG:4326', transform=tr,
                       tiled=True, blockxsize=512, blockysize=512,
                       compress='DEFLATE', predictor=2, BIGTIFF='YES') as ds:
        ds.write(arr, 1)
        ds.build_overviews([2, 4, 8, 16], rasterio.enums.Resampling.average)


def run(tag, dry=False):
    t0 = time.time()
    l = layer(tag)
    if l is None:
        print(f"  no downtown layer for {tag}"); return
    print(f"\n===== downtown {tag} =====", flush=True)
    arr, geo = as_grid(l)
    mod = validate.modern_for(geo['bbox'], geo['W'], geo['H'])
    print(f"  {geo['W']}x{geo['H']} @ {MPP} m/px  coverage {(arr>0).mean():.2f}", flush=True)
    src = P('mosaics', f'downtown_{tag}_src.tif')
    write_tif(src, arr, geo)
    t = gridval.prepare_cached(arr, mod, MPP, f'dt{tag}_src')
    prior = gridval.Prior(gridval.coarse_field(t, win_m=1400.0, overlap=0.4,
                                               search_m=120.0, min_ratio=1.05,
                                               log=print), MPP)
    before = gridval.grid_hier_prep(t, NY=4, NX=4, prior=prior, min_valid=0.4)
    report(before, f'downtown {tag} before')
    kept, tr, R, warp_at = warpsolve.solve(t, geo['bbox'], dtmap.MLAT, dtmap.MLON,
                                           win_m=700.0, overlap=0.4, iters=4,
                                           min_valid=0.45,
                                           lengths=(200., 350., 500., 800.))
    if kept is None or len(kept) < 12:
        print("  too little control downtown; leaving it alone"); return
    d = np.array([math.hypot(s['dE'], s['dN']) for s in kept])
    print(f"  {len(kept)} stations correcting median {np.median(d):.1f} "
          f"p90 {np.percentile(d,90):.1f} max {d.max():.1f} m", flush=True)
    if dry:
        return
    out = P('mosaics', f'downtown_{tag}_fix.tif')
    rbfapply.apply(src, geo, warp_at, out, f'/tmp/dt_{tag}.raw',
                   step_m=max(40.0, R.length / 8), log=lambda s: None)
    geo2 = dict(geo); json.dump(geo2, open(P('data', f'dt{tag}_fix_geo.json'), 'w'))
    arr2 = rasterio.open(out).read(1)
    t2 = gridval.prepare_cached(arr2, mod, MPP, f'dt{tag}_fix_{int(t0)}')
    prior2 = gridval.Prior(gridval.coarse_field(t2, win_m=1400.0, overlap=0.4,
                                                search_m=120.0, min_ratio=1.05,
                                                log=lambda *_: None), MPP)
    after = gridval.grid_hier_prep(t2, NY=4, NX=4, prior=prior2, min_valid=0.4)
    report(before, f'downtown {tag} before')
    report(after, f'downtown {tag} after ')
    print(f"  wrote {out}  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    dry = '--dry' in sys.argv
    tags = [a for a in sys.argv[1:] if not a.startswith('--')] or ['1949', '1956', '1961']
    for tg in tags:
        run(tg, dry)
