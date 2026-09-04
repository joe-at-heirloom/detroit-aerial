#!/usr/bin/env python
"""Fetch a modern reference raster covering every block.

Source: Esri World Imagery (the basemap this project already credits).

The important part is the reprojection. Tiles arrive in Web Mercator, where the
north-south scale grows as 1/cos(latitude); this project's grid is linear in
latitude. Across the 33 km the blocks span, that difference is about 0.5%, or
160 m end to end -- more than an order of magnitude larger than the accuracy being
chased. Stamping a lat/lon bounding box onto a Mercator mosaic would therefore
build that stretch into the reference and into every control point measured
against it, so each output row is resampled from its own true Mercator latitude.

Usage:  ./.venv/bin/python scripts/fetchmodern.py [--z 16] [--out modern_west_hi]
"""
import sys, os, math, json, time, io
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import urllib.request
from PIL import Image
import rasterio
from rasterio.transform import from_origin
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import dtmap
Image.MAX_IMAGE_PIXELS = None
def P(*a): return os.path.join(ROOT, *a)

URL = ("https://services.arcgisonline.com/ArcGIS/rest/services/"
       "World_Imagery/MapServer/tile/{z}/{y}/{x}")
UA = "detroit-air-survey/1.0 (historical aerial georeferencing; contact via repo)"
CACHE = '/tmp/das_tiles'
TS = 256


def arg(n, d):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def lon2x(lon, z): return (lon + 180.0) / 360.0 * (1 << z)
def lat2y(lat, z):
    r = math.radians(lat)
    return (1.0 - math.log(math.tan(r) + 1.0 / math.cos(r)) / math.pi) / 2.0 * (1 << z)
def y2lat(y, z):
    n = math.pi * (1 - 2 * y / (1 << z))
    return math.degrees(math.atan(math.sinh(n)))


def fetch(z, x, y, tries=4):
    p = os.path.join(CACHE, str(z), str(x), f'{y}.jpg')
    if os.path.exists(p) and os.path.getsize(p) > 0:
        try:
            return np.asarray(Image.open(p).convert('L'))
        except Exception:
            pass
    os.makedirs(os.path.dirname(p), exist_ok=True)
    for k in range(tries):
        try:
            rq = urllib.request.Request(URL.format(z=z, x=x, y=y),
                                        headers={'User-Agent': UA})
            with urllib.request.urlopen(rq, timeout=30) as r:
                b = r.read()
            im = Image.open(io.BytesIO(b)).convert('L')
            with open(p, 'wb') as f:
                f.write(b)
            return np.asarray(im)
        except Exception:
            time.sleep(0.4 * (k + 1))
    return None


def main():
    z = int(arg('--z', 16))
    name = arg('--out', 'modern_west_hi')
    # --bbox S,W,N,E covers ground no block reaches yet. Locating a roll of film
    # that overlaps nothing already solved needs a reference wider than the blocks,
    # and the collection holds ten 1949 rolls over the whole county.
    if '--bbox' in sys.argv:
        S, W, N, E = [float(x) for x in arg('--bbox', '').split(',')]
        print(f"bbox {S:.4f}..{N:.4f} N, {W:.4f}..{E:.4f} W", flush=True)
        return _build(z, name, S, W, N, E)
    # union of every block's extent, padded
    boxes = []
    for t in ('1949', '1956', '1961', '1967'):
        for suf in ('pre', 'rbf', 'final', 'v2'):
            f = P('data', f'{t}_{suf}_geo.json')
            if os.path.exists(f):
                boxes.append(json.load(open(f))['bbox']); break
    pad = 0.004
    S = min(b[0] for b in boxes) - pad; W = min(b[1] for b in boxes) - pad
    N = max(b[2] for b in boxes) + pad; E = max(b[3] for b in boxes) + pad
    print(f"union {S:.4f}..{N:.4f} N, {W:.4f}..{E:.4f} W", flush=True)
    return _build(z, name, S, W, N, E)


def _build(z, name, S, W, N, E):
    tx0 = int(math.floor(lon2x(W, z))); tx1 = int(math.ceil(lon2x(E, z)))
    ty0 = int(math.floor(lat2y(N, z))); ty1 = int(math.ceil(lat2y(S, z)))
    nx, ny = tx1 - tx0, ty1 - ty0
    print(f"z{z}: {nx} x {ny} = {nx*ny} tiles", flush=True)

    merc = np.zeros((ny * TS, nx * TS), np.uint8)
    todo = [(x, y) for y in range(ty0, ty1) for x in range(tx0, tx1)]
    done = [0]; missing = [0]
    def work(t):
        x, y = t
        a = fetch(z, x, y)
        if a is None:
            missing[0] += 1
        else:
            merc[(y - ty0) * TS:(y - ty0 + 1) * TS, (x - tx0) * TS:(x - tx0 + 1) * TS] = a
        done[0] += 1
        if done[0] % 400 == 0:
            print(f"  {done[0]}/{len(todo)} tiles", flush=True)
    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(work, todo))
    print(f"  fetched {len(todo)-missing[0]}/{len(todo)} ({missing[0]} missing)", flush=True)

    # target grid: linear in lat/lon, at the tiles' own ground resolution
    mpp = 156543.03392 * math.cos(math.radians((S + N) / 2)) / (1 << z)
    Wp = int((E - W) * dtmap.MLON / mpp); Hp = int((N - S) * dtmap.MLAT / mpp)
    print(f"  resampling to {Wp} x {Hp} @ {mpp:.2f} m/px (linear lat/lon)", flush=True)
    lon = W + (np.arange(Wp) + 0.5) * (E - W) / Wp
    sx = lon2x(lon, z) - tx0
    sx = np.clip(sx * TS, 0, merc.shape[1] - 1.001)
    x0 = sx.astype(np.int32); fx = (sx - x0).astype(np.float32)[None, :]
    out = np.zeros((Hp, Wp), np.uint8)
    band = 2048
    for j0 in range(0, Hp, band):
        j1 = min(Hp, j0 + band)
        lat = N - (np.arange(j0, j1) + 0.5) * (N - S) / Hp
        sy = np.array([lat2y(v, z) for v in lat]) - ty0
        sy = np.clip(sy * TS, 0, merc.shape[0] - 1.001)
        y0 = sy.astype(np.int32); fy = (sy - y0).astype(np.float32)[:, None]
        A = merc[y0][:, x0].astype(np.float32); B = merc[y0][:, x0 + 1].astype(np.float32)
        C = merc[y0 + 1][:, x0].astype(np.float32); D = merc[y0 + 1][:, x0 + 1].astype(np.float32)
        out[j0:j1] = np.clip(A * (1 - fx) * (1 - fy) + B * fx * (1 - fy)
                             + C * (1 - fx) * fy + D * fx * fy, 0, 255).astype(np.uint8)
    del merc
    tif = P('mosaics', f'{name}.tif')
    tr = from_origin(W, N, (E - W) / Wp, (N - S) / Hp)
    with rasterio.open(tif, 'w', driver='GTiff', height=Hp, width=Wp, count=1,
                       dtype='uint8', crs='EPSG:4326', transform=tr, tiled=True,
                       blockxsize=512, blockysize=512, compress='DEFLATE', predictor=2,
                       num_threads='ALL_CPUS', BIGTIFF='YES') as ds:
        for y in range(0, Hp, 2048):
            yb = min(Hp, y + 2048)
            ds.write(out[y:yb], 1, window=rasterio.windows.Window(0, y, Wp, yb - y))
        ds.build_overviews([2, 4, 8, 16, 32], rasterio.enums.Resampling.average)
    json.dump(dict(bbox=[S, W, N, E], W=Wp, H=Hp, mpp=mpp, zoom=z),
              open(P('data', f'{name}_geo.json'), 'w'))
    print(f"  wrote {tif} ({os.path.getsize(tif)/1e6:.0f} MB)", flush=True)


if __name__ == '__main__':
    main()
