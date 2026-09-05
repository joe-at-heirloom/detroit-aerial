#!/usr/bin/env python
"""Fetch an already-orthorectified imagery layer onto this project's grid.

Two kinds of source, one output: a 3-band GeoTIFF in EPSG:4326 with the same
linear lat/lon transform and square-metre pixels every other layer uses, plus
data/layer_<name>_geo.json, which manifest.py picks up and serve.py reads.

  arcgis  an ArcGIS ImageServer (the City of Detroit's 1998 / 2005 / 2010
          orthos). Exported chunk by chunk in EPSG:4326 at square-degree pixels,
          so the requested pixel grid has exactly the bbox's aspect -- an
          ImageServer silently widens a bbox whose aspect does not match, which
          once put 243 m of phantom error into this project (HANDOFF, ruled out
          #7). Every returned extent is checked against the request anyway, and
          the chunk is then warped onto the project grid.
  naip    USDA NAIP through Microsoft Planetary Computer's STAC (2012-2022 over
          Detroit). Quarter-quad COGs in UTM 17N, warped onto the project grid.

The extent is the union of the four served West Detroit blocks, padded, so the
wipe reaches every corner of the film. Outside a source's own coverage the layer
is black, which serve.py already draws as transparent.

  ./.venv/bin/python scripts/fetchlayer.py 1998 --arcgis https://egis.detroitmi.gov/image/rest/services/Imagery/1998_Aerial_Imagery/ImageServer --mpp 1.0
  ./.venv/bin/python scripts/fetchlayer.py 2018 --naip 2018 --mpp 0.6
  [--bbox S,W,N,E]  a small extent for a trial run
"""
import sys, os, json, math, time, io, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from rasterio.warp import reproject, transform_bounds, Resampling
from rasterio.crs import CRS
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import dtmap
def P(*a): return os.path.join(ROOT, *a)

UA = {'User-Agent': 'detroit-air-survey/1.0 (historical aerial georeferencing; contact via repo)'}
CHUNK = 4000            # output pixels per side per request / warp band
PAD = 0.004             # degrees, same as fetchmodern


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def target_extent():
    if '--bbox' in sys.argv:
        return [float(x) for x in arg('--bbox').split(',')]
    man = json.load(open(P('data', 'manifest.json')))
    boxes = [l['bbox'] for l in man['layers'] if l['group'] == 'west' and l['id'].startswith('b')]
    return [min(b[0] for b in boxes) - PAD, min(b[1] for b in boxes) - PAD,
            max(b[2] for b in boxes) + PAD, max(b[3] for b in boxes) + PAD]


def grid(S, W, N, E, mpp):
    Wp = int((E - W) * dtmap.MLON / mpp); Hp = int((N - S) * dtmap.MLAT / mpp)
    return Wp, Hp, from_origin(W, N, (E - W) / Wp, (N - S) / Hp)


def http(url, data=None, tries=5, timeout=300):
    for k in range(tries):
        try:
            rq = urllib.request.Request(url, data=data, headers=dict(UA, **({'Content-Type': 'application/json'} if data else {})))
            with urllib.request.urlopen(rq, timeout=timeout) as r:
                return r.read()
        except Exception as exc:
            if k == tries - 1:
                raise
            time.sleep(3 * (k + 1))


# ------------------------------------------------------------------ arcgis
def arcgis_info(base):
    """The service's extent in lat/lon (None if its CRS cannot be parsed) and the
    largest request it accepts, in pixels."""
    d = json.loads(http(base + '?f=json'))
    mx = (int(d.get('maxImageWidth') or 4000), int(d.get('maxImageHeight') or 4000))
    try:
        e = d['extent']; sr = e['spatialReference']
        crs = CRS.from_wkt(sr['wkt']) if 'wkt' in sr else CRS.from_epsg(sr.get('latestWkid') or sr['wkid'])
        w, s, e_, n = transform_bounds(crs, CRS.from_epsg(4326), e['xmin'], e['ymin'], e['xmax'], e['ymax'])
        return [s, w, n, e_], mx
    except Exception as exc:
        print(f"  service extent unknown ({exc}); requesting everything", flush=True)
        return None, mx


def arcgis_chunk(base, S, W, N, E, mpp):
    """One chunk, exported at square-degree pixels, returned with its transform."""
    d = mpp / dtmap.MLAT                               # degrees per pixel, both axes
    w = max(1, int(round((E - W) / d))); h = max(1, int(round((N - S) / d)))
    q = dict(bbox=f"{W},{S},{E},{N}", bboxSR=4326, imageSR=4326, size=f"{w},{h}",
             format='tiff', pixelType='U8', noData=0, interpolation='RSP_BilinearInterpolation', f='image')
    if '--where' in sys.argv:
        # a mosaic dataset holding several vintages (USGS's NAIP service carries a
        # Year attribute): pick one, e.g. --where "Year=2022"
        q['mosaicRule'] = json.dumps(dict(mosaicMethod='esriMosaicAttribute', where=arg('--where'),
                                          ascending=True, mosaicOperation='MT_FIRST'))
    for k in range(4):
        data = http(base + '/exportImage?' + urllib.parse.urlencode(q))
        if data[:4] in (b'II*\x00', b'MM\x00*'):
            break
        if k == 3:
            raise RuntimeError(f"not a TIFF for {W},{S},{E},{N} at {w}x{h}: {data[:200]!r}")
        time.sleep(5 * (k + 1))                       # an error page: usually transient
    with MemoryFile(data) as mf:
        with mf.open() as ds:
            b = ds.bounds
            # the extent the server actually returned must be the one asked for
            tol = 2 * d
            if abs(b.left - W) > tol or abs(b.right - E) > tol or abs(b.bottom - S) > tol or abs(b.top - N) > tol:
                raise RuntimeError(f"server widened the extent: asked {W},{S},{E},{N} got {b}")
            arr = ds.read(list(range(1, min(3, ds.count) + 1)))
            if arr.shape[0] == 1:
                arr = np.repeat(arr, 3, axis=0)
            return arr, ds.transform, ds.crs


def build_arcgis(name, base, mpp, S, W, N, E):
    Wp, Hp, tr = grid(S, W, N, E, mpp)
    cover, (mxw, mxh) = arcgis_info(base)
    if '--maxreq' in sys.argv:
        # some servers cap the RESPONSE, not the pixel count: imagery.michigan.gov
        # returns an error page for a 4-band TIFF over ~25 MB, so ask for less
        # (2400,2400 there) than its advertised 15000 x 4100
        rw, rh = [int(x) for x in arg('--maxreq').split(',')]
        mxw, mxh = min(mxw, rw), min(mxh, rh)
    # a chunk is requested at square-degree pixels, so it is MLAT/MLON (~1.35x)
    # wider in request pixels than in output pixels; keep inside the service's cap
    cx = min(CHUNK, int((mxw - 8) * dtmap.MLON / dtmap.MLAT)); cy = min(CHUNK, mxh - 8)
    print(f"  {Wp} x {Hp} @ {mpp} m/px in {cx} x {cy} chunks; service covers {cover}", flush=True)
    dst = out_dataset(name, Wp, Hp, tr)
    jobs = [(i, j) for j in range(0, Hp, cy) for i in range(0, Wp, cx)]
    done = [0]; skipped = [0]; t0 = time.time()

    def work(ij):
        i, j = ij
        cw = min(cx, Wp - i); ch = min(cy, Hp - j)
        cW = W + (E - W) * i / Wp; cE = W + (E - W) * (i + cw) / Wp
        cN = N - (N - S) * j / Hp; cS = N - (N - S) * (j + ch) / Hp
        if cover and (cE <= cover[1] or cW >= cover[3] or cN <= cover[0] or cS >= cover[2]):
            skipped[0] += 1; done[0] += 1
            return
        # ask for a hair more than the chunk so the warp's edge pixels have neighbours
        m = 2 * mpp / dtmap.MLAT
        arr, str_, scrs = arcgis_chunk(base, cS - m, cW - m, cN + m, cE + m, mpp)
        out = np.zeros((3, ch, cw), np.uint8)
        ctr = from_origin(cW, cN, (E - W) / Wp, (N - S) / Hp)
        reproject(arr, out, src_transform=str_, src_crs=scrs, dst_transform=ctr,
                  dst_crs=CRS.from_epsg(4326), resampling=Resampling.bilinear, src_nodata=None)
        with wlock:
            dst.write(out, window=rasterio.windows.Window(i, j, cw, ch))
        done[0] += 1
        if done[0] % 10 == 0:
            print(f"  {done[0]}/{len(jobs)} chunks ({skipped[0]} outside coverage) {time.time()-t0:.0f}s", flush=True)

    import threading
    wlock = threading.Lock()
    with ThreadPoolExecutor(max_workers=int(arg('--workers', 3))) as ex:
        list(ex.map(work, jobs))
    finish(dst, name, S, W, N, E, Wp, Hp, mpp, dict(kind='arcgis', url=base))


# -------------------------------------------------------------------- naip
STAC = 'https://planetarycomputer.microsoft.com/api/stac/v1/search'
TOKEN = 'https://planetarycomputer.microsoft.com/api/sas/v1/token/naip'


def naip_items(year, S, W, N, E):
    body = json.dumps(dict(collections=['naip'], bbox=[W, S, E, N],
                           datetime=f"{year}-01-01/{year}-12-31", limit=500)).encode()
    feats = json.loads(http(STAC, data=body))['features']
    return [(f['id'], f['assets']['image']['href'], f['bbox']) for f in feats]


# The quarter-quads are cloud-optimised GeoTIFFs on Azure: ~150 MB each at 1 m,
# ~500 MB at 0.6 m, 34 of them over the extent, and this machine's link is about
# 10 MB/s. Reading only the overlapping windows through GDAL's range requests was
# tried and ran ten times slower per byte than a plain download, so the quads are
# streamed whole to a cache (deleted by the caller once the layer is written) on
# several connections. The SAS token lasts an hour and is refreshed per file.
def build_naip(name, year, mpp, S, W, N, E):
    import shutil
    items = [it for it in naip_items(year, S, W, N, E)
             if it[2][1] < N and it[2][3] > S and it[2][0] < E and it[2][2] > W]
    if not items:
        raise SystemExit(f"no NAIP {year} over the extent")
    Wp, Hp, tr = grid(S, W, N, E, mpp)
    print(f"  {Wp} x {Hp} @ {mpp} m/px; {len(items)} quarter-quads for {year}", flush=True)
    cache = arg('--cache', f"/tmp/das_naip/{year}")
    os.makedirs(cache, exist_ok=True)
    state = dict(tok=None, at=0)

    def token():
        if time.time() - state['at'] > 30 * 60:
            state['tok'] = json.loads(http(TOKEN))['token']; state['at'] = time.time()
        return state['tok']

    def fetch(it):
        iid, href, _ = it
        p = os.path.join(cache, iid + '.tif')
        if os.path.exists(p) and os.path.getsize(p) > 1_000_000:
            return p
        for k in range(5):
            try:
                rq = urllib.request.Request(href + '?' + token(), headers=UA)
                with urllib.request.urlopen(rq, timeout=600) as r, open(p + '.part', 'wb') as f:
                    shutil.copyfileobj(r, f, 1 << 20)
                os.replace(p + '.part', p)
                return p
            except Exception as exc:
                if k == 4:
                    raise
                time.sleep(5 * (k + 1))
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=int(arg('--workers', 6))) as ex:
        paths = list(ex.map(fetch, items))
    gb = sum(os.path.getsize(p) for p in paths) / 1e9
    print(f"  {len(paths)} quads, {gb:.1f} GB, {time.time()-t0:.0f}s ({gb*1e3/max(1,time.time()-t0):.1f} MB/s)", flush=True)

    dst = out_dataset(name, Wp, Hp, tr)
    srcs = [(rasterio.open(p), b) for p, (_, _, b) in zip(paths, items)]
    dlon = (E - W) / Wp; dlat = (N - S) / Hp
    for j in range(0, Hp, CHUNK):
        ch = min(CHUNK, Hp - j)
        bN = N - dlat * j; bS = N - dlat * (j + ch)
        out = np.zeros((3, ch, Wp), np.uint8)

        def one(src):
            ds, (w_, s_, e_, n_) = src
            if n_ <= bS or s_ >= bN or e_ <= W or w_ >= E:
                return None
            # only the columns this quad can touch, with a pixel of margin
            i0 = max(0, int((w_ - W) / dlon) - 1); i1 = min(Wp, int(math.ceil((e_ - W) / dlon)) + 1)
            tmp = np.zeros((3, ch, i1 - i0), np.uint8)
            reproject(rasterio.band(ds, [1, 2, 3]), tmp,
                      dst_transform=from_origin(W + i0 * dlon, bN, dlon, dlat), dst_crs=CRS.from_epsg(4326),
                      resampling=Resampling.bilinear, src_nodata=0, dst_nodata=0, num_threads=2)
            return i0, tmp
        with ThreadPoolExecutor(max_workers=4) as ex:
            for r in ex.map(one, srcs):
                if r is None:
                    continue
                i0, tmp = r
                # paste only where the quad has data, so the black margin of one
                # quad never overwrites its neighbour
                m = tmp.max(axis=0) > 0
                sub = out[:, :, i0:i0 + tmp.shape[2]]
                sub[:, m] = tmp[:, m]
        dst.write(out, window=rasterio.windows.Window(0, j, Wp, ch))
        print(f"  rows {j+ch}/{Hp} {time.time()-t0:.0f}s", flush=True)
    for ds, _ in srcs:
        ds.close()
    finish(dst, name, S, W, N, E, Wp, Hp, mpp, dict(kind='naip', year=year, items=[i[0] for i in items]))


# ------------------------------------------------------------------ output
def out_dataset(name, Wp, Hp, tr):
    tif = P('mosaics', f'layer_{name}.tif')
    return rasterio.open(tif, 'w', driver='GTiff', height=Hp, width=Wp, count=3, dtype='uint8',
                         crs='EPSG:4326', transform=tr, tiled=True, blockxsize=512, blockysize=512,
                         compress='DEFLATE', predictor=2, photometric='RGB',
                         num_threads='ALL_CPUS', BIGTIFF='YES')


def finish(dst, name, S, W, N, E, Wp, Hp, mpp, source):
    print("  building overviews", flush=True)
    dst.build_overviews([2, 4, 8, 16, 32, 64], Resampling.average)
    dst.close()
    tif = P('mosaics', f'layer_{name}.tif')
    json.dump(dict(bbox=[S, W, N, E], W=Wp, H=Hp, mpp=mpp, label=arg('--label', name), source=source),
              open(P('data', f'layer_{name}_geo.json'), 'w'))
    print(f"  wrote {tif} ({os.path.getsize(tif)/1e6:.0f} MB)", flush=True)


def main():
    name = sys.argv[1]
    S, W, N, E = target_extent()
    mpp = float(arg('--mpp', 1.0))
    print(f"layer {name}: {S:.4f}..{N:.4f} N, {W:.4f}..{E:.4f} W", flush=True)
    if '--arcgis' in sys.argv:
        build_arcgis(name, arg('--arcgis').rstrip('/'), mpp, S, W, N, E)
    elif '--naip' in sys.argv:
        build_naip(name, arg('--naip'), mpp, S, W, N, E)
    else:
        raise SystemExit('--arcgis URL or --naip YEAR')


if __name__ == '__main__':
    main()
