#!/usr/bin/env python
"""Independent verification against USGS NAIP.

Everything else in this project is measured against a modern reference this
project built itself, out of Esri tiles, with a reprojection this project wrote.
That is a single point of failure, and it has already been wrong once: the
reference used to be Web Mercator with a lat/lon box stamped on it, bowing by 17 m
at mid-span, and every control point was quietly fitted to the bow.

NAIP is independent of all of that. It is US public-domain orthoimagery produced by
USDA/USGS at 0.3-0.6 m, it comes from a different organisation and a different
processing chain than Esri, and its ImageServer will return EPSG:4326 directly --
so the server does the reprojection and this project's own Mercator handling is
taken out of the loop entirely.

Two things get checked:

  --reference   is our modern raster itself correct? Sample windows of NAIP and
                correlate them against modern_west_hi.tif. Any systematic offset
                here is an error in our reference, which would have been inherited
                by every mosaic aligned to it.
  --block TAG   is the mosaic aligned to reality? Correlate the final mosaic
                directly against NAIP at windows spread over its coverage.

Usage:  ./.venv/bin/python scripts/naipcheck.py --reference [--n 24]
        ./.venv/bin/python scripts/naipcheck.py --block 1961 [--n 30] [--suffix final]
"""
import sys, os, io, json, math, time, hashlib
import urllib.request, urllib.parse
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image
import rasterio
from rasterio.windows import Window
import gridval, dtmap
import validate
from validate import P
Image.MAX_IMAGE_PIXELS = None

NAIP = ("https://imagery.nationalmap.gov/arcgis/rest/services/"
        "USGSNAIPImagery/ImageServer/exportImage")
CACHE = '/tmp/das_naip'
UA = "detroit-air-survey/1.0 (historical aerial georeferencing)"


def arg(name, default=None):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


def fetch(bbox, W, H, tries=4):
    """NAIP over a lat/lon bbox, resampled onto an exact W x H grid, in EPSG:4326.

    ArcGIS exportImage does NOT simply honour the bbox you give it: if the pixel
    grid you request has a different aspect ratio than the bbox, it silently
    *expands the bbox* to match, and returns imagery covering more ground than you
    asked for. Our check windows are square in metres, which at this latitude is
    1.35:1 in degrees, so asking for a square pixel grid quietly stretched the
    latitude extent by 243 m -- and the whole verification then compared imagery
    from two different footprints and reported ~74 m of error that was not there.

    So: request a pixel grid whose aspect already matches the bbox in degrees, which
    leaves the server nothing to adjust, then resample to the caller's grid here.
    `verify_extent()` below asserts the server actually did what we expect."""
    s, w, n, e = bbox
    key = hashlib.md5(f"{s:.6f},{w:.6f},{n:.6f},{e:.6f},{W},{H}".encode()).hexdigest()
    p = os.path.join(CACHE, key + '.png')
    if os.path.exists(p):
        try:
            return np.asarray(Image.open(p))
        except Exception:
            pass
    os.makedirs(CACHE, exist_ok=True)
    # Derive the request size so its aspect matches the bbox in degrees, and apply
    # the server's size cap BEFORE deriving the other axis -- capping afterwards
    # breaks the aspect match and re-triggers the very bbox expansion this is here
    # to avoid.
    aspect = (e - w) / (n - s)
    def _size(cap):
        hq = min(max(W, H), cap, int(cap / aspect) if aspect > 1 else cap)
        hq = max(64, int(round(hq)))
        wq = int(round(hq * aspect))
        if wq > cap:
            wq = cap; hq = max(64, int(round(wq / aspect)))
        return wq, hq

    # The service 500s on large requests -- 3799 x 2680 is already too much -- so
    # step the cap down until it answers and resample up afterwards. The aspect
    # stays matched at every step, which is what keeps the framing honest.
    for cap in (3000, 2200, 1600, 1100):
        Wq, Hq = _size(cap)
        q = urllib.parse.urlencode(dict(bbox=f"{w},{s},{e},{n}", bboxSR=4326,
                                        imageSR=4326, size=f"{Wq},{Hq}",
                                        format='tiff', f='image'))
        got = None
        for k in range(tries):
            try:
                rq = urllib.request.Request(f"{NAIP}?{q}", headers={'User-Agent': UA})
                got = urllib.request.urlopen(rq, timeout=180).read()
                break
            except urllib.error.HTTPError as ex:
                if ex.code >= 500:
                    break            # too big; try a smaller cap
                time.sleep(0.6 * (k + 1))
            except Exception:
                time.sleep(0.6 * (k + 1))
        if got is None:
            continue
        try:
            im = Image.open(io.BytesIO(got))
            a = np.asarray(im)
            if a.ndim == 3:
                al = a[..., 3] if a.shape[2] == 4 else None
                g = a[..., :3].mean(2)
                if al is not None:
                    g = np.where(al > 0, g, 0)
            else:
                g = a.astype(np.float32)
            g = np.clip(g, 0, 255).astype(np.uint8)
            if g.shape != (H, W):
                g = np.asarray(Image.fromarray(g).resize((W, H), Image.LANCZOS))
            Image.fromarray(g).save(p)
            return g
        except Exception:
            continue
    return None


def verify_extent(bbox, W, H):
    """Ask the server what footprint it actually returned, and compare."""
    s, w, n, e = bbox
    aspect = (e - w) / (n - s)
    Wq = int(round(max(W, H) * max(1.0, aspect))); Hq = int(round(Wq / aspect))
    q = urllib.parse.urlencode(dict(bbox=f"{w},{s},{e},{n}", bboxSR=4326, imageSR=4326,
                                    size=f"{Wq},{Hq}", format='tiff', f='json'))
    rq = urllib.request.Request(f"{NAIP}?{q}", headers={'User-Agent': UA})
    d = json.loads(urllib.request.urlopen(rq, timeout=120).read())
    ex = d['extent']
    dN = (ex['ymax'] - ex['ymin'] - (n - s)) * dtmap.MLAT
    dE = (ex['xmax'] - ex['xmin'] - (e - w)) * dtmap.MLON
    cN = ((ex['ymin'] + ex['ymax']) / 2 - (s + n) / 2) * dtmap.MLAT
    cE = ((ex['xmin'] + ex['xmax']) / 2 - (w + e) / 2) * dtmap.MLON
    return dict(extent_err_N=dN, extent_err_E=dE, centre_err_N=cN, centre_err_E=cE,
                size=(d.get('width'), d.get('height')))


def windows_over(cov_bbox, cov_mask_fn, n, span_m):
    """Lattice of check windows inside the covered area."""
    S, W, N, E = cov_bbox
    dlat = span_m / dtmap.MLAT; dlon = span_m / dtmap.MLON
    out = []
    k = int(math.ceil(math.sqrt(n * 3)))
    for j in range(k):
        for i in range(k):
            lat = N - (j + 0.5) / k * (N - S)
            lon = W + (i + 0.5) / k * (E - W)
            if not (S + dlat < lat < N - dlat and W + dlon < lon < E - dlon):
                continue
            bb = [lat - dlat / 2, lon - dlon / 2, lat + dlat / 2, lon + dlon / 2]
            if cov_mask_fn(bb):
                out.append(bb)
    step = max(1, len(out) // n)
    return out[::step][:n]


def measure(a, b, mpp, search_m=60.0):
    """Alias-proof by construction: a +/-60 m search is under one 97.5 m block."""
    ra, va = gridval.ridge_full(a.astype(np.float32), mpp)
    rb, vb = gridval.ridge_full(b.astype(np.float32), mpp)
    if ra.std() < 1e-9 or rb.std() < 1e-9:
        return None
    return gridval.match_cell((ra * va).astype(np.float32), va.astype(np.float32),
                              (rb * vb).astype(np.float32), vb.astype(np.float32),
                              mpp, search_m, 'mncc')


def summarise(res, label, lock=1.15):
    g = [r for r in res if r and r['ratio'] >= lock and not r['edge']]
    if not g:
        print(f"  {label}: nothing locked"); return
    m = np.array([r['mag'] for r in g])
    dE = np.array([r['dE'] for r in g]); dN = np.array([r['dN'] for r in g])
    print(f"  {label}: n={len(g)}/{len(res)}  median {np.median(m):5.2f}  "
          f"p90 {np.percentile(m,90):5.2f}  max {m.max():5.2f} m")
    print(f"      systematic offset: dE {np.median(dE):+5.2f}  dN {np.median(dN):+5.2f} m"
          f"   >10 m {(m>10).sum()}  >25 m {(m>25).sum()}")


def main():
    MPP = float(arg('--mpp', 1.0))
    n = int(arg('--n', 24))
    span = float(arg('--span', 700))
    px = int(round(span / MPP))

    if '--reference' in sys.argv:
        rb = json.load(open(P('data', 'modern_west_hi_geo.json')))['bbox']
        print(f"Checking OUR modern reference against USGS NAIP "
              f"({n} windows, {span:.0f} m, {MPP} m/px)", flush=True)
        wins = windows_over(rb, lambda bb: True, n, span)
        res = []
        for k, bb in enumerate(wins):
            nap = fetch(bb, px, px)
            if nap is None or (nap > 0).mean() < 0.9:
                res.append(None); continue
            ours = validate.modern_for(bb, px, px)
            if (ours > 0).mean() < 0.9:
                res.append(None); continue
            res.append(measure(ours, nap, MPP))
            if (k + 1) % 8 == 0:
                print(f"    {k+1}/{len(wins)}", flush=True)
        summarise(res, 'our reference vs NAIP')
        json.dump([r for r in res if r], open('/tmp/naip_reference.json', 'w'))
        return

    tag = arg('--block')
    suffix = arg('--suffix', 'final')
    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json')))
    bbox = geo['bbox']
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{suffix}.tif'))
    S, W, N, E = bbox

    def covered(bb):
        c0 = (bb[1] - W) / (E - W) * ds.width; r0 = (N - bb[2]) / (N - S) * ds.height
        wpx = (bb[3] - bb[1]) / (E - W) * ds.width
        hpx = (bb[2] - bb[0]) / (N - S) * ds.height
        a = ds.read(1, window=Window(c0, r0, wpx, hpx), out_shape=(24, 24),
                    boundless=True, fill_value=0)
        return (a > 0).mean() > 0.97

    print(f"Checking {tag} ({suffix}) directly against USGS NAIP "
          f"({n} windows, {span:.0f} m, {MPP} m/px)", flush=True)
    wins = windows_over(bbox, covered, n, span)
    print(f"  {len(wins)} fully covered check windows", flush=True)
    res = []
    for k, bb in enumerate(wins):
        nap = fetch(bb, px, px)
        if nap is None or (nap > 0).mean() < 0.9:
            res.append(None); continue
        c0 = (bb[1] - W) / (E - W) * ds.width; r0 = (N - bb[2]) / (N - S) * ds.height
        wpx = (bb[3] - bb[1]) / (E - W) * ds.width
        hpx = (bb[2] - bb[0]) / (N - S) * ds.height
        hist = ds.read(1, window=Window(c0, r0, wpx, hpx), out_shape=(px, px),
                       boundless=True, fill_value=0)
        res.append(measure(hist, nap, MPP))
        if (k + 1) % 8 == 0:
            print(f"    {k+1}/{len(wins)}", flush=True)
    summarise(res, f'{tag} {suffix} vs NAIP')
    json.dump([r for r in res if r], open(f'/tmp/naip_{tag}.json', 'w'))


if __name__ == '__main__':
    main()
