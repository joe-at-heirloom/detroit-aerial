#!/usr/bin/env python
"""Rebuild a block end to end, correcting the frames themselves.

  pre-warp composite -> per-frame (dE, dN, rot) against modern imagery
                     -> re-composite -> residual RBF field -> final mosaic

The per-frame step is the one that matters for the worst cells. A continuous warp
cannot represent the step change in error that happens where two negatives meet,
so it averages across the seam and leaves half the disagreement on either side.
Correcting each frame before compositing removes the step; the RBF afterwards only
has the smooth within-frame part left to do.

Usage:  ./.venv/bin/python scripts/rebuild.py 1961 [--rots] [--skip-frames]
                                             [--frame-iters 3]
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, stations, rbfwarp, rbfapply, warpsolve, frameadjust, dtmap
import mosaic_sim, apply_anchor
from PIL import Image
import rasterio
from rasterio.transform import from_origin
from validate import load, modern_for, report, MPP, P
from rewarp import build_prewarp, SP, BUNDLE, NATIVE
Image.MAX_IMAGE_PIXELS = None


def placements(tag):
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
    return apply_anchor.corrected(sol, an['anchor'], an['cE'], an['cN']), ppm


def composite(tag, sol, ppm, label, log=print):
    mpp = 1.0 / ppm
    os.environ['MOSAIC_RAW'] = f'/tmp/comp_{tag}.raw'
    m = mosaic_sim.build(sol, f"{SP}/fullres", ppm, mpp, chunk=768, log=lambda s: None)
    la_top = dtmap.LAT0 + m['maxN'] / dtmap.MLAT
    lo_left = dtmap.LON0 + m['minE'] / dtmap.MLON
    geo = dict(minE=float(m['minE']), maxN=float(m['maxN']), W=m['W'], H=m['H'], mpp=mpp,
               bbox=[dtmap.LAT0 + (m['maxN'] - m['H'] * mpp) / dtmap.MLAT, lo_left, la_top,
                     dtmap.LON0 + (m['minE'] + m['W'] * mpp) / dtmap.MLON])
    path = P('mosaics', f'detroit_{tag}_{label}.tif')
    tr = from_origin(lo_left, la_top, mpp / dtmap.MLON, mpp / dtmap.MLAT)
    with rasterio.open(path, 'w', driver='GTiff', height=m['H'], width=m['W'], count=1,
                       dtype='uint8', crs='EPSG:4326', transform=tr, tiled=True,
                       blockxsize=512, blockysize=512, compress='DEFLATE', predictor=2,
                       num_threads='ALL_CPUS', BIGTIFF='YES') as ds:
        for y in range(0, m['H'], 2048):
            yb = min(m['H'], y + 2048)
            ds.write(np.asarray(m['arr'][y:yb]), 1,
                     window=rasterio.windows.Window(0, y, m['W'], yb - y))
        ds.build_overviews([2, 4, 8, 16, 32, 64], rasterio.enums.Resampling.average)
    del m
    try: os.remove(f'/tmp/comp_{tag}.raw')
    except OSError: pass
    json.dump(geo, open(P('data', f'{tag}_{label}_geo.json'), 'w'))
    log(f"  composited {label}: {geo['W']}x{geo['H']} @ {mpp:.3f} m/px")
    return geo


def run(tag, rots=False, skip_frames=False, frame_iters=3):
    t0 = time.time()
    print(f"\n===== rebuild {tag} =====", flush=True)
    sol, ppm = placements(tag)
    print(f"  {len(sol)} frames", flush=True)

    if not os.path.exists(P('data', f'{tag}_pre_geo.json')):
        composite(tag, sol, ppm, 'pre')
    arr, mod, bbox = load(tag, 'pre')
    t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_pre')
    prior = gridval.Prior(gridval.coarse_field(t), MPP)
    before = gridval.grid_hier_prep(t, prior=prior)
    report(before, 'pre-warp composite')

    if not skip_frames:
        geo = json.load(open(P('data', f'{tag}_pre_geo.json')))
        minE, maxN = geo['minE'], geo['maxN']
        angles = tuple(np.arange(-0.8, 0.81, 0.4)) if rots else (0.0,)
        # Iterate the per-frame solve. A frame is rendered from the negative with
        # whatever correction it currently carries, and matched against the modern
        # ridge map -- which never changes -- so this costs nothing but the solve:
        # no re-compositing until the end. It matters because a frame that failed
        # to lock on the first pass did so with a prior tens of metres off, and
        # borrowing its neighbours' answer is strictly worse than solving its own.
        for it in range(frame_iters):
            pr = prior if it == 0 else None      # later passes start from ~0
            ang = angles if (rots and it == frame_iters - 1) else (0.0,)
            print(f"  per-frame solve, pass {it+1}/{frame_iters} "
                  f"({len(ang)} angle(s))", flush=True)
            fx = frameadjust.solve_frames(sol, f"{SP}/fullres", t, minE, maxN, MPP,
                                          prior=pr, rots=ang)
            if not fx:
                print("    nothing locked; stopping", flush=True); break
            d = np.array([math.hypot(v['dE'], v['dN']) for v in fx.values()])
            print(f"    {len(fx)}/{len(sol)} frames locked; this pass moves them "
                  f"median {np.median(d):.1f}  p90 {np.percentile(d,90):.1f}  "
                  f"max {d.max():.1f} m", flush=True)
            fx = frameadjust.regularise(fx, sol)
            for r, v in fx.items():
                sol[r]['dE'] += v['dE']; sol[r]['dN'] += v['dN']
                sol[r]['rot'] += v['rot']
            if np.median(d) < 1.5:
                print("    frames converged", flush=True); break
        json.dump({r: dict(dE=sol[r]['dE'], dN=sol[r]['dN'], rot=sol[r]['rot'])
                   for r in sol}, open(P('data', f'frameadj_{tag}.json'), 'w'))
        tot = np.array([math.hypot(sol[r]['dE'], sol[r]['dN']) for r in sol])
        print(f"  total per-frame correction: median {np.median(tot):.1f}  "
              f"p90 {np.percentile(tot,90):.1f}  max {tot.max():.1f} m", flush=True)
        composite(tag, sol, ppm, 'fa')
        arr, mod, bbox = load(tag, 'fa')
        t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_fa')
        prior = gridval.Prior(gridval.coarse_field(t), MPP)
        mid = gridval.grid_hier_prep(t, prior=prior)
        report(mid, 'after per-frame     ')
        src, geo = 'fa', json.load(open(P('data', f'{tag}_fa_geo.json')))
    else:
        src, geo = 'pre', json.load(open(P('data', f'{tag}_pre_geo.json')))

    kept, tr, R, warp_at = warpsolve.solve(t, bbox, dtmap.MLAT, dtmap.MLON,
                                           win_m=2000.0, overlap=0.4, iters=5,
                                           lengths=(250., 400., 600., 800., 1200.))
    if kept is None:
        print("  no residual field"); return
    json.dump(kept, open(P('data', f'stations_final_{tag}.json'), 'w'))
    json.dump(dict(kind='poly1_rbf', trend=tr, length=R.length),
              open(P('data', f'rbffit_final_{tag}.json'), 'w'))
    R.dump(P('data', f'rbfres_final_{tag}.json'))
    g = dict(minE=geo['minE'], maxN=geo['maxN'], W=geo['W'], H=geo['H'], mpp=geo['mpp'])
    out = rbfapply.apply(P('mosaics', f'detroit_{tag}_{src}.tif'), g, warp_at,
                         P('mosaics', f'detroit_{tag}_final.tif'),
                         f'/tmp/fin_{tag}.raw', step_m=max(50.0, R.length / 8),
                         log=lambda s: None)
    json.dump(out, open(P('data', f'{tag}_final_geo.json'), 'w'))
    arr2, mod2, _ = load(tag, 'final')
    t2 = gridval.prepare_cached(arr2, mod2, MPP, f'{tag}_final')
    prior2 = gridval.Prior(gridval.coarse_field(t2, log=lambda *_: None), MPP)
    after = gridval.grid_hier_prep(t2, prior=prior2)
    json.dump(after, open(P('data', f'gridval_{tag}.json'), 'w'))
    report(before, 'pre-warp composite  ')
    if not skip_frames:
        report(mid, 'after per-frame     ')
    report(after, 'FINAL               ')
    print(f"  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    rots = '--rots' in sys.argv
    sk = '--skip-frames' in sys.argv
    it = int(sys.argv[sys.argv.index('--frame-iters') + 1]) if '--frame-iters' in sys.argv else 3
    for tg in [a for a in sys.argv[1:] if not a.startswith('--')
               and a != str(it)]:
        run(tg, rots, sk, it)
