#!/usr/bin/env python
"""Every layer against every other layer, on one grid, with one instrument.

`dtcross.py` asks the question a wipe asks -- how far apart are the two layers on
screen -- but only for downtown, and by reading whole rasters into memory. This
asks it for every pair of layers in a manifest group, which is what the alignment
plan's measurement system is (ALIGNMENT_PLAN.md, step 2):

  1. Each layer is read ONCE onto a shared lat/lon grid (windowed GDAL reads with
     area averaging, so a 6 GB NAIP year never sits in memory whole), the hand
     alignment in data/adjust.json applied exactly as the tile server applies it,
     and ridge-filtered at both of gridval's bands. The ridge maps are cached on
     disk by source signature; a changed raster or adjustment rebuilds its maps.
  2. Each PAIR gets gridval's regional coarse field on 6 km arterial windows --
     a 1-2 km cell has no unique arterial peak and would alias a whole block --
     then the alias-proof fine match on fixed geographic cells about that prior.
  3. Each pair is re-measured with a planted field (translation + 1.5 m/km scale
     + 1 mrad rotation) on one side. The measurement must move by exactly that
     field; how far it does not is the instrument's error for that pair, and it
     is reported next to the number rather than assumed.

Nothing here fits or moves anything.

Usage:  ./.venv/bin/python scripts/crossmatrix.py west [--mpp 2.0] [--cell 1300] [--jobs 8]
            [--layers b1961,l2016,...] [--no-adjust] [--no-selfcheck] [--out FILE]
        ./.venv/bin/python scripts/crossmatrix.py downtown --mpp 1.25 --cell 800 --min-ratio 1.3
"""
import sys, os, json, math, time, argparse, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import gridval, handadjust, dtmap
Image = None

P = lambda *a: os.path.join(ROOT, *a)
CACHE = P('runs', 'cache')
# the modern reference validate.py measures against; not a viewer layer, but
# every absolute number in HANDOFF is relative to it
EXTRA = {'west': [dict(id='esrihi', label='Esri hi', file='mosaics/modern_west_hi.tif',
                       geo='data/modern_west_hi_geo.json')]}
PLANT = dict(dE=30.0, dN=-20.0, scale=1.5e-3, rot=1.0e-3)   # metres, -, radians
TRACK_MAX = 5.0     # metres: a cell that misses its planted field by more is dropped


def load_layers(group, only=None, extra=True):
    man = json.load(open(P('data', 'manifest.json')))
    L = [dict(l) for l in man['layers'] if l['group'] == group]
    if extra:
        for x in EXTRA.get(group, []):
            if os.path.exists(P(x['file'])):
                L.append(dict(id=x['id'], label=x['label'], file=x['file'],
                              bbox=json.load(open(P(x['geo'])))['bbox']))
    if only:
        L = [l for l in L if l['id'] in only]
    return L


def grid_for(layers, mpp):
    s = min(l['bbox'][0] for l in layers); w = min(l['bbox'][1] for l in layers)
    n = max(l['bbox'][2] for l in layers); e = max(l['bbox'][3] for l in layers)
    W = int(round((e - w) * dtmap.MLON / mpp)); H = int(round((n - s) * dtmap.MLAT / mpp))
    W -= W % 2; H -= H % 2            # the coarse maps are a 2x downsample
    return dict(bbox=[s, w, n, e], W=W, H=H, mpp=mpp)


def _gray(a):
    if a.ndim == 2:
        return a
    g = (0.299 * a[0] + 0.587 * a[1] + 0.114 * a[2])
    g[(a[:3] == 0).all(0)] = 0
    return np.clip(np.round(g), 0, 255).astype(np.uint8)


def read_onto(l, g, adj):
    """The layer on the grid, 0 = nodata. GeoTIFFs are placed by their own
    transform (what the tile server reads); plain images by the manifest bbox."""
    import rasterio
    from rasterio.warp import reproject, Resampling
    from rasterio.transform import from_bounds
    s, w, n, e = g['bbox']; H, W = g['H'], g['W']
    out = np.zeros((H, W), np.uint8)
    path = P(l['file'])
    if adj is not None:
        # hand-aligned: sample through handadjust's inverse, as dtcross and the
        # tile server do. Only downtown carries adjustments and it is small.
        src, (S_, W_, N_, E_) = _full(path, l)
        h, w_ = src.shape
        for y0 in range(0, H, 1024):
            y1 = min(H, y0 + 1024)
            lat = n - (np.arange(y0, y1) + 0.5) * (n - s) / H
            lon = w + (np.arange(W) + 0.5) * (e - w) / W
            LON, LAT = np.meshgrid(lon, lat)
            lo, la = handadjust.source_lonlat(LON, LAT, tuple(l['bbox']), adj)
            sx = (lo - W_) / (E_ - W_) * w_ - 0.5; sy = (N_ - la) / (N_ - S_) * h - 0.5
            ok = (sx >= 0) & (sx < w_) & (sy >= 0) & (sy < h)
            blk = np.zeros(ok.shape, np.uint8)
            blk[ok] = src[sy[ok].astype(np.int32), sx[ok].astype(np.int32)]
            out[y0:y1] = blk
        return out
    if path.lower().endswith(('.tif', '.tiff')):
        with rasterio.open(path) as ds:
            nb = min(ds.count, 3)
            step = 1024
            for y0 in range(0, H, step):
                y1 = min(H, y0 + step)
                tr = from_bounds(w, n - y1 * (n - s) / H, e, n - y0 * (n - s) / H, W, y1 - y0)
                dst = np.zeros((nb, y1 - y0, W), np.uint8)
                reproject(source=rasterio.band(ds, list(range(1, nb + 1))), destination=dst,
                          src_transform=ds.transform, src_crs=ds.crs, dst_transform=tr,
                          dst_crs=ds.crs, resampling=Resampling.average,
                          src_nodata=0, dst_nodata=0, num_threads=2)
                out[y0:y1] = _gray(dst if nb > 1 else dst[0])
        return out
    src, (S_, W_, N_, E_) = _full(path, l)
    tr_src = from_bounds(W_, S_, E_, N_, src.shape[1], src.shape[0])
    reproject(source=src, destination=out, src_transform=tr_src, src_crs='EPSG:4326',
              dst_transform=from_bounds(w, s, e, n, W, H), dst_crs='EPSG:4326',
              resampling=Resampling.average, src_nodata=0, dst_nodata=0)
    return out


def _full(path, l):
    global Image
    if path.lower().endswith(('.tif', '.tiff')):
        import rasterio
        with rasterio.open(path) as ds:
            b = ds.bounds
            a = ds.read(list(range(1, min(ds.count, 3) + 1)))
            return _gray(a if a.shape[0] > 1 else a[0]), (b.bottom, b.left, b.top, b.right)
    if Image is None:
        from PIL import Image as _I; _I.MAX_IMAGE_PIXELS = None; Image = _I
    return np.asarray(Image.open(path).convert('L')), tuple(l['bbox'])


def signature(l, g, adj):
    st = os.stat(P(l['file']))
    return dict(file=l['file'], size=st.st_size, mtime=int(st.st_mtime),
                bbox=l['bbox'], adj=adj, grid=g, v=2)


def build_maps(args):
    """Worker: ridge maps for one layer, from cache when the signature matches."""
    l, g, adj, cdir = args
    sig = signature(l, g, adj)
    meta = os.path.join(cdir, l['id'] + '.json')
    if os.path.exists(meta) and json.load(open(meta)) == sig:
        return l['id'], 'cached', 0.0
    t0 = time.time()
    I = read_onto(l, g, adj)
    mpp = g['mpp']
    f, v = gridval.ridge_full(I.astype(np.float32), mpp, (18.0, 46.0), 320.0, strip=3072)
    Id = gridval._ds(I.astype(np.float32), 2)
    c, vc = gridval.ridge_full(Id, mpp * 2, gridval.COARSE_BAND, gridval.COARSE_FLAT, strip=3072)
    for nm, a in (('f', f), ('v', v), ('c', c), ('vc', vc)):
        np.save(os.path.join(cdir, f"{l['id']}.{nm}.npy"),
                a.astype(np.float32) if nm in ('f', 'c') else a.astype(bool))
    np.save(os.path.join(cdir, f"{l['id']}.img.npy"), I)
    json.dump(sig, open(meta, 'w'))
    return l['id'], f'built, coverage {float((I > 0).mean()):.2f}', time.time() - t0


def maps(cdir, lid):
    return {nm: np.load(os.path.join(cdir, f'{lid}.{nm}.npy'), mmap_mode='r')
            for nm in ('f', 'v', 'c', 'vc')}


def overlap_px(la, lb, g):
    """Pixel box of the two layers' bbox intersection, even-aligned for ds=2."""
    s, w, n, e = g['bbox']; H, W = g['H'], g['W']
    S = max(la['bbox'][0], lb['bbox'][0]); Wl = max(la['bbox'][1], lb['bbox'][1])
    N = min(la['bbox'][2], lb['bbox'][2]); E = min(la['bbox'][3], lb['bbox'][3])
    if S >= N or Wl >= E:
        return None
    y0 = int((n - N) / (n - s) * H) // 2 * 2; y1 = min(H, -(-int(math.ceil((n - S) / (n - s) * H)) // 2) * 2)
    x0 = int((Wl - w) / (e - w) * W) // 2 * 2; x1 = min(W, -(-int(math.ceil((E - w) / (e - w) * W)) // 2) * 2)
    return y0, y1, x0, x1


def plant_at(y, x, cy, cx, mpp):
    """The planted field at pixel (y, x): translation + scale + rotation about the
    pair's centre, in metres E/N."""
    dx = (x - cx) * mpp; dy = -(y - cy) * mpp
    s, r = PLANT['scale'], PLANT['rot']
    return (PLANT['dE'] + s * dx - r * dy, PLANT['dN'] + s * dy + r * dx)


def measure_pair(args):
    a, b, g, cdir, cell_m, min_valid, selfcheck = args
    t0 = time.time()
    box = overlap_px(a, b, g)
    rec = dict(a=a['id'], b=b['id'])
    if box is None:
        return dict(rec, skip='no overlap')
    y0, y1, x0, x1 = box
    A, B = maps(cdir, a['id']), maps(cdir, b['id'])
    t = dict(mpp=g['mpp'], ds=2, shape=(y1 - y0, x1 - x0),
             hf=A['f'][y0:y1, x0:x1], hv=A['v'][y0:y1, x0:x1],
             hc=A['c'][y0 // 2:y1 // 2, x0 // 2:x1 // 2], hvc=A['vc'][y0 // 2:y1 // 2, x0 // 2:x1 // 2],
             mf=B['f'][y0:y1, x0:x1], mv=B['v'][y0:y1, x0:x1],
             mc=B['c'][y0 // 2:y1 // 2, x0 // 2:x1 // 2], mvc=B['vc'][y0 // 2:y1 // 2, x0 // 2:x1 // 2])
    both = np.asarray(t['hv'][::8, ::8]) & np.asarray(t['mv'][::8, ::8])
    rec['overlap_km2'] = float(both.sum() * (8 * g['mpp']) ** 2 / 1e6)
    if rec['overlap_km2'] < 1.0:
        return dict(rec, skip='overlap under 1 km2')
    # coarse field, accepted only when most windows lock (gridval.prior_for's rule)
    pts = gridval.coarse_field(t, log=lambda *_: None)
    nwin = getattr(gridval.coarse_field, 'last_n_windows', 0) or 0
    use = len(pts) >= 12 and len(pts) >= 0.6 * max(nwin, 1)
    prior = gridval.Prior(pts if use else [], g['mpp'])
    V = np.array([[p['dE'], p['dN']] for p in pts]) if pts else np.zeros((0, 2))
    rec['coarse'] = dict(windows=nwin, locked=len(pts), used=use,
                         median_dE=float(np.median(V[:, 0])) if len(V) else None,
                         median_dN=float(np.median(V[:, 1])) if len(V) else None,
                         max_mag=float(np.hypot(V[:, 0], V[:, 1]).max()) if len(V) else None)
    H, W = t['shape']; cp = max(8, int(round(cell_m / g['mpp'])))
    cy, cx = H / 2, W / 2
    s, w, n, e = g['bbox']
    cells = []
    for yy in range(0, H - cp // 2, cp):
        for xx in range(0, W - cp // 2, cp):
            ya, yb, xa, xb = yy, min(H, yy + cp), xx, min(W, xx + cp)
            va = float(np.asarray(t['hv'][ya:yb:4, xa:xb:4]).mean())
            vb = float(np.asarray(t['mv'][ya:yb:4, xa:xb:4]).mean())
            my, mx = (ya + yb) / 2, (xa + xb) / 2
            c = dict(lat=n - (y0 + my) * (n - s) / g['H'], lon=w + (x0 + mx) * (e - w) / g['W'],
                     va=round(va, 3), vb=round(vb, 3))
            if min(va, vb) < min_valid:
                c['skip'] = 'coverage'; cells.append(c); continue
            pr = prior.at(my, mx) if use else None
            r = gridval.match_hier_prep(t, ya, yb, xa, xb, prior=pr)
            if r is None:
                c['skip'] = 'nopeak'; cells.append(c); continue
            c.update({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
            if selfcheck:
                pE, pN = plant_at(my, mx, cy, cx, g['mpp'])
                # match_hier_prep adds hshift to a supplied prior itself
                r2 = gridval.match_hier_prep(t, ya, yb, xa, xb, prior=pr, hshift=(pE, pN))
                if r2 is not None and not r2['pegged']:
                    c['track'] = round(math.hypot(r2['dE'] - r['dE'] - pE, r2['dN'] - r['dN'] - pN), 3)
            cells.append(c)
    rec['cells'] = cells
    rec['secs'] = round(time.time() - t0, 1)
    return rec


def summarise(rec, min_ratio):
    if 'skip' in rec:
        return None
    tried = [c for c in rec['cells'] if 'skip' not in c or c['skip'] == 'nopeak']
    lk = [c for c in rec['cells'] if 'skip' not in c and not c['pegged'] and c['ratio'] >= min_ratio]
    tr = np.array([c['track'] for c in lk if 'track' in c])
    # a cell that did not follow its planted field was measuring something other
    # than the ground (downtown against Today: roofs); its number is not evidence
    ok = [c for c in lk if c.get('track', 0.0) <= TRACK_MAX]
    d = np.array([c['mag'] for c in ok])
    if not len(d):
        return dict(n=0, tried=len(tried), untracked=len(lk) - len(ok))
    return dict(n=len(d), tried=len(tried), untracked=len(lk) - len(ok), median=float(np.median(d)),
                p90=float(np.percentile(d, 90)), max=float(d.max()),
                over10=int((d > 10).sum()), over25=int((d > 25).sum()),
                bias_dE=float(np.median([c['dE'] for c in ok])),
                bias_dN=float(np.median([c['dN'] for c in ok])),
                track_median=float(np.median(tr)) if len(tr) else None,
                track_p90=float(np.percentile(tr, 90)) if len(tr) else None,
                track_over5=int((tr > 5).sum()) if len(tr) else None)


def resummarise(path):
    d = json.load(open(path))
    for rec in d['pairs']:
        rec['summary'] = summarise(rec, d['min_ratio'])
    json.dump(d, open(path, 'w'))


def main():
    if len(sys.argv) == 3 and sys.argv[1] == '--resummarise':
        return resummarise(sys.argv[2])
    ap = argparse.ArgumentParser()
    ap.add_argument('group')
    ap.add_argument('--mpp', type=float, default=2.0)
    ap.add_argument('--cell', type=float, default=1300.0, help='fine cell size, metres')
    ap.add_argument('--min-valid', type=float, default=0.45)
    ap.add_argument('--min-ratio', type=float, default=1.15)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--map-jobs', type=int, default=3, help='ridge builds are memory-heavy')
    ap.add_argument('--layers')
    ap.add_argument('--no-adjust', action='store_true')
    ap.add_argument('--no-selfcheck', action='store_true')
    ap.add_argument('--out')
    a = ap.parse_args()
    from concurrent.futures import ProcessPoolExecutor

    L = load_layers(a.group, a.layers.split(',') if a.layers else None)
    g = grid_for(L, a.mpp)
    ADJ = {} if a.no_adjust else handadjust.load(ROOT)
    tagadj = 'noadj' if a.no_adjust else 'adj'
    cdir = os.path.join(CACHE, f"xm_{a.group}_{a.mpp}")
    os.makedirs(cdir, exist_ok=True)
    print(f"{a.group}: {len(L)} layers, grid {g['W']}x{g['H']} @ {a.mpp} m/px, cache {cdir}", flush=True)
    for l in L:
        adj = handadjust.for_layer(ADJ, l['id'])
        l['adj'] = adj
        if adj:
            print(f"  {l['id']}: hand alignment {adj}", flush=True)
    # an adjusted layer's maps must not be confused with the same layer unadjusted
    for l in L:
        if l['adj'] is not None:
            l['id'] = l['id'] + '+adj'
    t0 = time.time()
    with ProcessPoolExecutor(a.map_jobs) as ex:
        for lid, what, secs in ex.map(build_maps, [(dict(l), g, l['adj'], cdir) for l in L]):
            print(f"  maps {lid:10s} {what}  [{secs:.0f}s]", flush=True)
    print(f"  ridge maps ready [{time.time()-t0:.0f}s]\n", flush=True)

    pairs = list(itertools.combinations(L, 2))
    out = dict(group=a.group, grid=g, cell_m=a.cell, min_valid=a.min_valid, min_ratio=a.min_ratio,
               plant=PLANT, adjust=tagadj, when=time.strftime('%Y-%m-%d %H:%M'),
               layers=[dict(id=l['id'], label=l['label'], file=l['file'], bbox=l['bbox'], adj=l['adj'])
                       for l in L], pairs=[])
    with ProcessPoolExecutor(a.jobs) as ex:
        args = [(x, y, g, cdir, a.cell, a.min_valid, not a.no_selfcheck) for x, y in pairs]
        for rec in ex.map(measure_pair, args):
            rec['summary'] = summarise(rec, a.min_ratio)
            out['pairs'].append(rec)
            s = rec['summary']
            if s is None:
                print(f"  {rec['a']:>9} vs {rec['b']:<9} {rec['skip']}", flush=True)
            elif s['n'] == 0:
                print(f"  {rec['a']:>9} vs {rec['b']:<9} nothing locked of {s['tried']}", flush=True)
            else:
                tk = (f"   track {s['track_median']:.1f}/{s['track_p90']:.1f}"
                      if s['track_median'] is not None else '')
                print(f"  {rec['a']:>9} vs {rec['b']:<9} {s['n']:3d}/{s['tried']:3d}  "
                      f"median {s['median']:5.1f}  p90 {s['p90']:5.1f}  max {s['max']:5.1f}  "
                      f">25 {s['over25']:2d}  bias {s['bias_dE']:+5.1f},{s['bias_dN']:+5.1f}"
                      f"{tk}  [{rec['secs']:.0f}s]", flush=True)
    path = a.out or P('runs', f"crossmatrix_{a.group}_{tagadj}.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(out, open(path, 'w'))
    print(f"\nwrote {path}  [{time.time()-t0:.0f}s]")


if __name__ == '__main__':
    main()
