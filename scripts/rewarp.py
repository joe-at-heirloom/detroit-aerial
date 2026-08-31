#!/usr/bin/env python
"""Re-solve a block's local warp on control that cannot alias.

Two modes:
  --source rbf   iterate on top of the existing pass-one mosaic (fast, no rebuild)
  --source raw   rebuild the pre-warp composite from the frames and warp it once
                 (slower, but one resampling instead of two)

Either way the field comes from pipeline/warpsolve.py, which never runs a
correlation search wider than half a city block.

Usage:  ./.venv/bin/python scripts/rewarp.py 1961 1967 [--source raw|rbf] [--dry]
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, stations, rbfwarp, rbfapply, warpsolve, dtmap
import mosaic_sim, apply_anchor
from PIL import Image
import rasterio
from rasterio.transform import from_origin
from validate import load, modern_for, report, MPP, P
Image.MAX_IMAGE_PIXELS = None

SP = ('/private/tmp/claude-501/-Users-joelint-Documents-Apps/'
      'b1877f49-3048-4d7e-83b9-c2a051af5bc9/scratchpad/detroit')
NATIVE = {'1961': 'native.json', '1967': 'native1967.json',
          '1949': 'native1949.json', '1956': 'native1956.json'}
BUNDLE = {'1961': ('sim_bundle.json', 'sim_anchor.json'),
          '1967': ('sim_bundle.json', 'sim_anchor.json'),
          '1949': ('sim_bundle_4956.json', 'sim_anchor_4956.json'),
          '1956': ('sim_bundle_4956.json', 'sim_anchor_4956.json')}


def build_prewarp(tag, log=print):
    """Rebuild the un-warped composite that pass one threw away."""
    bf, af = BUNDLE[tag]
    b = json.load(open(P('data', bf)))[tag]
    an = json.load(open(P('data', af)))[tag]
    c = json.load(open(P('data', f'calib{tag}.json')))
    nat = json.load(open(P('data', NATIVE[tag])))
    ppm = b['ppm_full']

    def to_dt(e, n):
        lat = c['lat0'] + n / c['MLAT']; lon = c['lon0'] + e / c['MLON']
        return ((lon - dtmap.LON0) * dtmap.MLON, (lat - dtmap.LAT0) * dtmap.MLAT)

    sol = {}
    for r in b['frames']:
        if r not in nat or not os.path.exists(f"{SP}/fullres/{r}.jpg"):
            continue
        e, n = to_dt(*b['pos'][r]); w, h = nat[r]
        sol[r] = dict(ok=True, e=e, n=n, rot=b['rot'][r], dE=0.0, dN=0.0,
                      gw=w / ppm, gh=h / ppm, ratio=2.0)
    fixed = apply_anchor.corrected(sol, an['anchor'], an['cE'], an['cN'])
    mpp = 1.0 / ppm
    log(f"  compositing {len(fixed)} frames at {mpp:.3f} m/px")
    os.environ['MOSAIC_RAW'] = f'/tmp/full_{tag}.raw'
    m = mosaic_sim.build(fixed, f"{SP}/fullres", ppm, mpp, chunk=768, log=lambda s: None)
    la_top = dtmap.LAT0 + m['maxN'] / dtmap.MLAT
    lo_left = dtmap.LON0 + m['minE'] / dtmap.MLON
    geo = dict(minE=float(m['minE']), maxN=float(m['maxN']), W=m['W'], H=m['H'], mpp=mpp,
               bbox=[dtmap.LAT0 + (m['maxN'] - m['H'] * mpp) / dtmap.MLAT, lo_left, la_top,
                     dtmap.LON0 + (m['minE'] + m['W'] * mpp) / dtmap.MLON])
    pre = P('mosaics', f'detroit_{tag}_pre.tif')
    tr = from_origin(lo_left, la_top, mpp / dtmap.MLON, mpp / dtmap.MLAT)
    with rasterio.open(pre, 'w', driver='GTiff', height=m['H'], width=m['W'], count=1,
                       dtype='uint8', crs='EPSG:4326', transform=tr, tiled=True,
                       blockxsize=512, blockysize=512, compress='DEFLATE', predictor=2,
                       num_threads='ALL_CPUS', BIGTIFF='YES') as ds:
        for y in range(0, m['H'], 2048):
            yb = min(m['H'], y + 2048)
            ds.write(np.asarray(m['arr'][y:yb]), 1,
                     window=rasterio.windows.Window(0, y, m['W'], yb - y))
        ds.build_overviews([2, 4, 8, 16, 32, 64], rasterio.enums.Resampling.average)
    del m
    for f in (f'/tmp/full_{tag}.raw',):
        try: os.remove(f)
        except OSError: pass
    json.dump(geo, open(P('data', f'{tag}_pre_geo.json'), 'w'))
    log(f"  wrote {pre}  {geo['W']}x{geo['H']}")
    return geo


def run(tag, source='rbf', dry=False):
    t0 = time.time()
    print(f"\n===== {tag}  (source: {source}) =====", flush=True)
    if source == 'raw':
        if not os.path.exists(P('data', f'{tag}_pre_geo.json')):
            build_prewarp(tag)
        suffix = 'pre'
    else:
        suffix = 'rbf'
    arr, mod, bbox = load(tag, suffix)
    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json')))
    print(f"  {arr.shape[1]}x{arr.shape[0]} @ {MPP} m/px  coverage {(arr>0).mean():.2f}",
          flush=True)
    t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_{suffix}')
    prior0 = gridval.Prior(gridval.coarse_field(t, log=lambda *_: None), MPP)
    before = gridval.grid_hier_prep(t, prior=prior0)
    report(before, 'before')

    kept, tr, R, warp_at, held_out = warpsolve.solve(
        t, bbox, dtmap.MLAT, dtmap.MLON, win_m=2000.0, overlap=0.4, iters=5,
        lengths=(400., 600., 800., 1200., 1800.))
    if kept is None:
        print("  no usable field"); return
    d = np.array([math.hypot(s['dE'], s['dN']) for s in kept])
    print(f"  final control: {len(kept)} stations correcting median {np.median(d):.1f} "
          f"p90 {np.percentile(d,90):.1f} max {d.max():.1f} m", flush=True)
    json.dump(kept, open(P('data', f'stations_v2_{tag}.json'), 'w'))
    json.dump(dict(kind='poly1_rbf', trend=tr, length=R.length),
              open(P('data', f'rbffit_v2_{tag}.json'), 'w'))
    R.dump(P('data', f'rbfres_v2_{tag}.json'))
    if dry:
        print(f"  [dry] {time.time()-t0:.0f}s"); return

    out_tif = P('mosaics', f'detroit_{tag}_v2.tif')
    g = dict(minE=geo['minE'], maxN=geo['maxN'], W=geo['W'], H=geo['H'], mpp=geo['mpp'])
    out = rbfapply.apply(P('mosaics', f'detroit_{tag}_{suffix}.tif'), g, warp_at,
                         out_tif, f'/tmp/v2_{tag}.raw',
                         step_m=max(50.0, R.length / 8), log=lambda s: None)
    json.dump(out, open(P('data', f'{tag}_v2_geo.json'), 'w'))
    print(f"  wrote {out_tif} [{time.time()-t0:.0f}s]", flush=True)

    arr2, mod2, bbox2 = load(tag, 'v2')
    t2 = gridval.prepare_cached(arr2, mod2, MPP, f'{tag}_v2_{int(t0)}')
    prior2 = gridval.Prior(gridval.coarse_field(t2, log=lambda *_: None), MPP)
    after = gridval.grid_hier_prep(t2, prior=prior2)
    report(before, 'before')
    report(after, 'after ')
    json.dump(after, open(P('data', f'gridval_{tag}_v2.json'), 'w'))
    print(f"  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    src = sys.argv[sys.argv.index('--source') + 1] if '--source' in sys.argv else 'rbf'
    dry = '--dry' in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    args = [a for a in args if a not in ('raw', 'rbf')]
    for tg in args:
        run(tg, src, dry)
