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
                                                  [--ref 1961] [--src]

`--src` re-runs from `downtown_<tag>_src.tif`, the uncorrected raster this script
wrote the first time. Without it a second run reads whatever the manifest points
at -- which, once a correction exists, is the CORRECTED mosaic, so the new warp
lands on top of the old one and the reported "before" is the previous run's
"after". Always pass it when redoing an epoch.
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
import validate, naipcheck
from validate import P, report
Image.MAX_IMAGE_PIXELS = None
MPP = 1.25          # downtown is small; measure it finer than the blocks


def layer(tag):
    man = json.load(open(P('data', 'manifest.json')))
    for l in man['layers']:
        if l['id'] == f'dt{tag}':
            return l
    return None


def from_src(tag, mpp=MPP):
    """The uncorrected raster on the shared grid, with the geometry it was cut on."""
    src = P('mosaics', f'downtown_{tag}_src.tif')
    geo = P('data', f'dt{tag}_fix_geo.json')
    if not (os.path.exists(src) and os.path.exists(geo)):
        raise SystemExit(f'--src needs downtown_{tag}_src.tif and dt{tag}_fix_geo.json')
    g = json.load(open(geo))
    if abs(g['mpp'] - mpp) > 1e-9:
        raise SystemExit(f"src grid is {g['mpp']} m/px, this run wants {mpp}")
    return rasterio.open(src).read(1), g


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


def reference(geo, ref):
    """NAIP by default; another downtown epoch when chaining.

    Downtown 1956 against modern imagery is the worst case in the whole project --
    the core was comprehensively rebuilt between 1956 and now, so most of what the
    correlation sees has no counterpart. Its control came back correcting 43 m and
    left the block at 32 m, which is what wrong control looks like. Downtown 1961 is
    five years away and now good to ~4 m, so chaining through it asks a question the
    imagery can actually answer."""
    if ref in (None, 'naip'):
        return naipcheck.fetch(geo['bbox'], geo['W'], geo['H'])
    src = P('mosaics', f'downtown_{ref}_fix.tif')
    if not os.path.exists(src):
        raise SystemExit(f'downtown {ref} has no corrected mosaic yet')
    ds = rasterio.open(src)
    a = ds.read(1, out_shape=(geo['H'], geo['W']))
    nap = naipcheck.fetch(geo['bbox'], geo['W'], geo['H'])
    if nap is None:
        return a
    return np.where(a > 0, a, nap).astype(np.uint8)


class Constant:
    """A Prior that returns the same shift everywhere."""

    def __init__(self, dE, dN):
        self.v = (float(dE), float(dN))

    def at(self, y, x):
        return self.v


def global_prior(t, search_m=250.0, min_ratio=1.30, min_mag=30.0, log=print):
    """One rigid shift for the whole scene, solved before the residual field.

    The fine search is capped at +/-40 m and walks about 120 m over three
    re-centrings. A mosaic further out than that cannot be reached -- and it does
    not fail loudly. Detroit's street lattice repeats every ~97.5 m, so the search
    locks onto the wrong member of it and returns a confident small number.
    Downtown 1956 has been reported at ~32 m for exactly that reason. Solved on
    the whole raster with a 250 m search it is **170 m** out, and all sixteen
    cells of a 4x4 grid then agree with that to within +/-20 m.

    A coarse *field* has no room downtown -- the scene is 3.4 km across and the
    field wants 6 km windows of mile-grid arterials, which is why `run` disables
    it. One global shift is a different question, and the whole raster is a single
    window with plenty of structure to answer it.

    Returns None when the peak is not convincing or the scene is already inside
    the fine search's reach, in which case a zero prior is the right answer.
    """
    H, W = t['shape']
    r = match_whole(t, search_m)
    if r is None:
        log("  global: no peak"); return None
    mag = math.hypot(r['dE'], r['dN'])
    ok = r['coarse_ratio'] >= min_ratio and mag >= min_mag
    log(f"  global: dE {r['dE']:+.1f} dN {r['dN']:+.1f} ({mag:.0f} m) "
        f"coarse ratio {r['coarse_ratio']:.2f} -> {'using it' if ok else 'ignored'}")
    return Constant(r['dE'], r['dN']) if ok else None


def match_whole(t, search_m):
    H, W = t['shape']
    return gridval.match_hier_prep(t, 0, H, 0, W, coarse_search=search_m)


def run(tag, dry=False, ref='naip', use_src=False):
    t0 = time.time()
    l = layer(tag)
    if l is None:
        print(f"  no downtown layer for {tag}"); return
    print(f"\n===== downtown {tag} (reference: {ref}"
          f"{', from src' if use_src else ''}) =====", flush=True)
    arr, geo = from_src(tag) if use_src else as_grid(l)
    # Downtown sits east of the modern reference raster, which was fetched for the
    # four West Detroit flight blocks and stops at -83.11; downtown runs -83.06 to
    # -83.02. Handing it that reference gives an all-zero modern side, which is
    # exactly what happened -- every cell reported "no peak" and the cause looked
    # like a matching failure. NAIP covers it, and being independent of Esri is a
    # bonus rather than a problem.
    mod = reference(geo, ref)
    if mod is None or (mod > 0).mean() < 0.5:
        print("  no reference imagery available over downtown"); return
    print(f"  {geo['W']}x{geo['H']} @ {MPP} m/px  coverage {(arr>0).mean():.2f}", flush=True)
    src = P('mosaics', f'downtown_{tag}_src.tif')
    if not use_src:
        write_tif(src, arr, geo)
    t = gridval.prepare_cached(arr, mod, MPP, f'dt{tag}_src_{ref}')
    # No coarse stage downtown. The regional coarse field is built on 6 km windows
    # of mile-grid arterials; the whole downtown raster is 3.4 km across, so there
    # is no room for one, and every attempt pegs at its search edge. It is not
    # needed either -- downtown is already within about 25 m, comfortably inside a
    # +/-40 m fine search, which is what the coarse stage exists to deliver.
    # A coarse field has no room here, but one rigid shift for the whole scene
    # does -- and without it an epoch further out than ~120 m is measured against
    # the wrong street lattice rather than reported as unreachable.
    pz = prior = global_prior(t) or gridval.Prior([], MPP)
    before = gridval.grid_hier_prep(t, NY=4, NX=4, prior=prior, min_valid=0.4)
    report(before, f'downtown {tag} before')
    kept, tr, R, warp_at, held_out = warpsolve.solve(t, geo['bbox'], dtmap.MLAT, dtmap.MLON,
                                           win_m=700.0, overlap=0.4, iters=4,
                                           min_valid=0.45, prior=pz,
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
    t2 = gridval.prepare_cached(arr2, mod, MPP, f'dt{tag}_fix_{ref}_{int(t0)}')
    # the corrected raster should need no prior; if it still does, the fix did not take
    prior2 = gridval.Prior([], MPP)
    after = gridval.grid_hier_prep(t2, NY=4, NX=4, prior=prior2, min_valid=0.4)
    report(before, f'downtown {tag} before')
    report(after, f'downtown {tag} after ')
    print(f"  wrote {out}  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    dry = '--dry' in sys.argv
    use_src = '--src' in sys.argv
    ref = sys.argv[sys.argv.index('--ref') + 1] if '--ref' in sys.argv else 'naip'
    tags = [a for a in sys.argv[1:] if not a.startswith('--') and a != ref] \
        or ['1949', '1956', '1961']
    for tg in tags:
        run(tg, dry, ref, use_src)
