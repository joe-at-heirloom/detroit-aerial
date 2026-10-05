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

  tiles   a Web Mercator tile cache (ArcGIS MapServer tile layer or any XYZ
          template with {z} {y} {x}): the City's MiSAIL caches, Michigan Tech's
          1951. Tiles at one zoom (--zoom, else the first as fine as --mpp) are
          mosaicked per chunk and warped onto the project grid.

  ./.venv/bin/python scripts/fetchlayer.py 1998 --arcgis https://egis.detroitmi.gov/image/rest/services/Imagery/1998_Aerial_Imagery/ImageServer --mpp 1.0
  ./.venv/bin/python scripts/fetchlayer.py 2018 --naip 2018 --mpp 0.6
  ./.venv/bin/python scripts/fetchlayer.py 1951 --tiles "https://tiles.arcgis.com/.../MapServer/tile/{z}/{y}/{x}" --mpp 0.6
  [--bbox S,W,N,E]  a small extent for a trial run
  [--clip]          store only the part of the extent the source covers (from its
                    service extent, or --cover S,W,N,E): Oakland County's years
                    reach only ~4 km past 8 Mile, so a full-extent file is 90 % black
  [--label "2010 NAIP"] [--prio -1] [--credit "USDA NAIP"]  see finish()
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
STRETCH = {}
BANDS = {}


def band_order(base, S, W, N, E):
    """Which raw bands are red, green and blue. The State's NAIP services do not
    agree -- 2005 is NIR,R,G,B and 2009 is B,G,R,NIR, where 2010-2020 are
    R,G,B,NIR -- and reading bands 1-3 blind made 2005 false-colour and swapped
    2009's red and blue. The service's own rendering knows: correlate a small
    raw export against its PNG at the extent's centre and take the best band for
    each channel. --bands 1,2,3 (0-based bandIds) overrides."""
    if '--bands' in sys.argv:
        BANDS[base] = arg('--bands'); return
    from PIL import Image
    la, lo = (S + N) / 2, (W + E) / 2; h = 0.004
    q = dict(bbox=f"{lo-h*1.35},{la-h},{lo+h*1.35},{la+h}", bboxSR=4326, imageSR=4326, size='270,200', f='image')
    try:
        raw = http(base + '/exportImage?' + urllib.parse.urlencode(dict(q, format='tiff', pixelType='U8')))
        png = http(base + '/exportImage?' + urllib.parse.urlencode(dict(q, format='png')))
        with MemoryFile(raw) as mf, mf.open() as ds:
            a = ds.read().astype(float)
        p = np.asarray(Image.open(io.BytesIO(png)).convert('RGB')).astype(float)
    except Exception as exc:
        print(f"  band order not checked ({exc})", flush=True); return
    if a.shape[0] < 3 or a.std() == 0:
        return
    best = [int(np.nanargmax([np.corrcoef(p[..., c].ravel(), a[b].ravel())[0, 1] for b in range(a.shape[0])]))
            for c in range(3)]
    if best != [0, 1, 2] and len(set(best)) == 3:
        BANDS[base] = ','.join(map(str, best))
        print(f"  service renders R,G,B from raw bands {best} (0-based): using bandIds {BANDS[base]}", flush=True)


def arcgis_info(base):
    """The service's extent in lat/lon (None if its CRS cannot be parsed) and the
    largest request it accepts, in pixels."""
    d = json.loads(http(base + '?f=json'))
    mx = (int(d.get('maxImageWidth') or 4000), int(d.get('maxImageHeight') or 4000))
    # A 16-bit service (Oakland County's 2008, values 10k-63k) asked for pixelType
    # U8 saturates to a white sheet. Have the server stretch it: a 0.5 % percent
    # clip over the dataset's own statistics (not per request, so chunks match).
    # Its native 4-band U16 also overflows the response limit, hence bandIds.
    if d.get('pixelType') not in (None, 'U8'):
        # output 1..255, not 0..255: the darkest picture must not become nodata
        STRETCH[base] = json.dumps(dict(rasterFunction='Stretch', outputPixelType='U8',
                                        rasterFunctionArguments=dict(StretchType=6, MinPercent=0.5,
                                                                     MaxPercent=0.5, DRA=False,
                                                                     Min=1, Max=255)))
        print(f"  {d['pixelType']} service: server-side percent-clip stretch to 8 bits", flush=True)
    try:
        e = d.get('extent') or d['fullExtent']; sr = e['spatialReference']      # MapServers say fullExtent
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
    if base.endswith('/MapServer'):
        return mapserver_chunk(base, q, W, S, E, N, d)
    if base in STRETCH:
        q['renderingRule'] = STRETCH[base]; q['bandIds'] = '0,1,2'
    elif base in BANDS:
        q['bandIds'] = BANDS[base]
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


_ms_checked = []


def mapserver_chunk(base, q, W, S, E, N, d):
    """A dynamic MapServer (the City of Taylor's yearly aerials) has /export, not
    /exportImage. Its f=json answer names an output file the City's server does not
    actually serve (404), so the picture comes from f=image; the extent the server
    draws is checked once, on the first chunk, through f=json -- the request
    keeps the bbox's aspect in degrees, which is what makes it honour the bbox.
    Transparent pixels are nodata."""
    from PIL import Image
    if '--native' in sys.argv:
        return mapserver_native(base, W, S, E, N, d)
    q = dict(bbox=q['bbox'], bboxSR=4326, imageSR=4326, size=q['size'], format='png32',
             transparent='true', dpi=96)
    if '--layers' in sys.argv:
        q['layers'] = arg('--layers')
    if not _ms_checked:
        e = json.loads(http(base + '/export?' + urllib.parse.urlencode(dict(q, f='json'))))['extent']
        tol = 2 * d
        if abs(e['xmin'] - W) > tol or abs(e['xmax'] - E) > tol or abs(e['ymin'] - S) > tol or abs(e['ymax'] - N) > tol:
            raise RuntimeError(f"server widened the extent: asked {W},{S},{E},{N} got {e}")
        _ms_checked.append(1)
    for k in range(4):
        data = http(base + '/export?' + urllib.parse.urlencode(dict(q, f='image')))
        if data[:8] == b'\x89PNG\r\n\x1a\n':
            break
        if k == 3:
            raise RuntimeError(f"not a PNG for {W},{S},{E},{N}: {data[:200]!r}")
        time.sleep(5 * (k + 1))
    a = np.asarray(Image.open(io.BytesIO(data)).convert('RGBA'))
    rgb = np.transpose(a[..., :3], (2, 0, 1)).copy()
    rgb[:, a[..., 3] == 0] = 0
    tr = from_origin(W, N, (E - W) / a.shape[1], (N - S) / a.shape[0])
    return rgb, tr, CRS.from_epsg(4326)


_ms_crs = {}


def mapserver_native(base, W, S, E, N, d):
    """--native: some services draw nothing when asked to reproject (Taylor's 2004
    returns an empty PNG in EPSG:4326 and the full picture in its own state
    plane). Ask in the service's own coordinates for the box around the chunk,
    with square pixels there, and let the chunk's warp do the reprojection."""
    from PIL import Image
    if base not in _ms_crs:
        sr = json.loads(http(base + '?f=json'))['spatialReference']
        _ms_crs[base] = (sr, CRS.from_wkt(sr['wkt']) if 'wkt' in sr else CRS.from_epsg(sr.get('latestWkid') or sr['wkid']))
    sr, crs = _ms_crs[base]
    x0, y0, x1, y1 = transform_bounds(CRS.from_epsg(4326), crs, W, S, E, N, densify_pts=21)
    unit = crs.linear_units_factor[1] if crs.is_projected else 1.0      # metres per unit
    px = d * dtmap.MLAT / unit                                          # same ground pixel, in native units
    w = max(1, int(round((x1 - x0) / px))); h = max(1, int(round((y1 - y0) / px)))
    srj = json.dumps({'wkt': sr['wkt']} if 'wkt' in sr else {'wkid': sr.get('latestWkid') or sr['wkid']})
    q = dict(bbox=f"{x0},{y0},{x1},{y1}", bboxSR=srj, imageSR=srj, size=f"{w},{h}", format='png32',
             transparent='true', dpi=96, f='image')
    for k in range(4):
        data = http(base + '/export?' + urllib.parse.urlencode(q))
        if data[:8] == b'\x89PNG\r\n\x1a\n':
            break
        if k == 3:
            raise RuntimeError(f"not a PNG for {W},{S},{E},{N}: {data[:200]!r}")
        time.sleep(5 * (k + 1))
    a = np.asarray(Image.open(io.BytesIO(data)).convert('RGBA'))
    rgb = np.transpose(a[..., :3], (2, 0, 1)).copy()
    rgb[:, a[..., 3] == 0] = 0
    return rgb, from_origin(x0, y1, (x1 - x0) / a.shape[1], (y1 - y0) / a.shape[0]), crs


def build_arcgis(name, base, mpp, S, W, N, E):
    Wp, Hp, tr = grid(S, W, N, E, mpp)
    cover, (mxw, mxh) = arcgis_info(base)
    if not base.endswith('/MapServer') and base not in STRETCH:
        band_order(base, S, W, N, E)
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
    g = dict(bbox=[S, W, N, E], W=Wp, H=Hp, mpp=mpp, label=arg('--label', name), source=source)
    # Two layers of one year (the City's 2010 and NAIP 2010): the year view puts
    # the higher --prio on top; --credit is the publisher named in its readout.
    if '--prio' in sys.argv:
        g['prio'] = float(arg('--prio'))
    if '--credit' in sys.argv:
        g['credit'] = arg('--credit')
    # --part-of 2002: more ground of the same picture, fetched separately (the
    # strip east of the West block to downtown). The year view draws it as one
    # layer with its parent -- one outline, one caption -- and the scene tabs leave
    # it out of their rails.
    if '--part-of' in sys.argv:
        g['part_of'] = 'l' + arg('--part-of')
    json.dump(g, open(P('data', f'layer_{name}_geo.json'), 'w'))
    print(f"  wrote {tif} ({os.path.getsize(tif)/1e6:.0f} MB)", flush=True)


# ------------------------------------------------------------------- tiles
# A tile cache (an ArcGIS MapServer/tile layer, or any XYZ template with {z} {y}
# {x}) on the standard Web Mercator grid: the City's MiSAIL 2020 / 2024 caches,
# Michigan Tech's 1951 layer. Each output chunk gathers the tiles under it at one
# zoom, mosaics them in EPSG:3857 and warps onto the project grid. A missing tile
# (404, or fully transparent) is nodata, i.e. black, like everywhere else.
HALF = 20037508.342789244


def tile_cover(tmpl):
    """The cache's full extent in lat/lon, if it is an ArcGIS MapServer."""
    if '/MapServer/tile/' not in tmpl and '/ImageServer/tile/' not in tmpl:
        return None
    base = tmpl.split('/tile/')[0]
    try:
        d = json.loads(http(base + '?f=json'))
        e = d.get('fullExtent') or d['extent']; sr = e['spatialReference']
        crs = CRS.from_epsg(sr.get('latestWkid') or sr['wkid'])
        w, s, e_, n = transform_bounds(crs, CRS.from_epsg(4326), e['xmin'], e['ymin'], e['xmax'], e['ymax'])
        return [s, w, n, e_]
    except Exception as exc:
        print(f"  cache extent unknown ({exc})", flush=True)
        return None


def get_tile(url):
    """bytes, or None when the cache has no tile there."""
    for k in range(5):
        try:
            rq = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(rq, timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (400, 404, 410):
                return None
            if k == 4:
                raise
            if exc.code in (403, 429, 503):
                # a small municipal server rate-limits (Taylor's answered 403
                # after ~1000 tiles): back right off rather than hammer it
                time.sleep(30 * (k + 1))
                continue
        except Exception:
            if k == 4:
                raise
        time.sleep(2 * (k + 1))


def tile_rgb(data):
    from PIL import Image
    im = Image.open(io.BytesIO(data))
    if im.mode in ('RGBA', 'LA', 'P', 'PA'):
        im = im.convert('RGBA')
        a = np.asarray(im)
        rgb = a[..., :3].copy()
        rgb[a[..., 3] == 0] = 0
    else:
        rgb = np.asarray(im.convert('RGB')).copy()
    return np.transpose(rgb, (2, 0, 1))


def build_tiles(name, tmpl, mpp, S, W, N, E):
    lat0 = math.radians((S + N) / 2)
    z = int(arg('--zoom', 0))
    if not z:
        # the coarsest zoom that is at least as fine as the output pixel
        z = next(z for z in range(10, 23) if 2 * HALF / 256 / 2 ** z * math.cos(lat0) <= mpp)
    tm = 2 * HALF / 2 ** z                                # tile size in metres (Mercator)
    Wp, Hp, tr = grid(S, W, N, E, mpp)
    cover = tile_cover(tmpl)
    cx = cy = 2048
    print(f"  {Wp} x {Hp} @ {mpp} m/px from z{z} tiles "
          f"({tm/256*math.cos(lat0):.2f} m/px); cache covers {cover}", flush=True)
    dst = out_dataset(name, Wp, Hp, tr)
    jobs = [(i, j) for j in range(0, Hp, cy) for i in range(0, Wp, cx)]
    done = [0]; ntiles = [0]; t0 = time.time()
    import threading
    wlock = threading.Lock()
    inner = ThreadPoolExecutor(max_workers=int(arg('--tile-workers', 12)))

    def mx(lon): return math.radians(lon) * 6378137.0
    def my(lat): return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2)) * 6378137.0

    def work(ij):
        i, j = ij
        cw = min(cx, Wp - i); ch = min(cy, Hp - j)
        cW = W + (E - W) * i / Wp; cE = W + (E - W) * (i + cw) / Wp
        cN = N - (N - S) * j / Hp; cS = N - (N - S) * (j + ch) / Hp
        if cover and (cE <= cover[1] or cW >= cover[3] or cN <= cover[0] or cS >= cover[2]):
            done[0] += 1
            return
        x0 = int((mx(cW) + HALF) // tm); x1 = int((mx(cE) + HALF) // tm)
        y0 = int((HALF - my(cN)) // tm); y1 = int((HALF - my(cS)) // tm)
        keys = [(x, y) for y in range(y0, y1 + 1) for x in range(x0, x1 + 1)]
        datas = list(inner.map(lambda k: get_tile(tmpl.replace('{z}', str(z)).replace('{x}', str(k[0]))
                                                       .replace('{y}', str(k[1]))), keys))
        if not any(datas):
            done[0] += 1
            return
        mos = np.zeros((3, (y1 - y0 + 1) * 256, (x1 - x0 + 1) * 256), np.uint8)
        for (x, y), d in zip(keys, datas):
            if d:
                try:
                    t = tile_rgb(d)
                except Exception:
                    continue
                mos[:, (y - y0) * 256:(y - y0 + 1) * 256, (x - x0) * 256:(x - x0 + 1) * 256] = t[:, :256, :256]
        str_ = from_origin(-HALF + x0 * tm, HALF - y0 * tm, tm / 256, tm / 256)
        out = np.zeros((3, ch, cw), np.uint8)
        ctr = from_origin(cW, cN, (E - W) / Wp, (N - S) / Hp)
        reproject(mos, out, src_transform=str_, src_crs=CRS.from_epsg(3857), dst_transform=ctr,
                  dst_crs=CRS.from_epsg(4326), resampling=Resampling.bilinear, src_nodata=0, dst_nodata=0)
        with wlock:
            dst.write(out, window=rasterio.windows.Window(i, j, cw, ch))
            ntiles[0] += sum(1 for d in datas if d)
        done[0] += 1
        if done[0] % 10 == 0:
            print(f"  {done[0]}/{len(jobs)} chunks, {ntiles[0]} tiles {time.time()-t0:.0f}s", flush=True)

    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(work, jobs))
    inner.shutdown()
    finish(dst, name, S, W, N, E, Wp, Hp, mpp, dict(kind='tiles', url=tmpl, zoom=z))


# -------------------------------------------------------------------- cogs
# A set of tiled GeoTIFFs on a web server, listed by a tile index (a shapefile,
# zipped or not, with a `url` field, as NOAA's Digital Coast publishes) or by a
# text file of URLs: NOAA's 2008 Great Lakes border ortho (489 tiles over the
# extent, 62 MB each at 0.3 m) and its 2025 / 2026 Detroit shoreline sets. Each
# tile is read straight off the server at the overview nearest --mpp, so a 0.6 m
# layer moves a quarter of the bytes, then warped into an output row band.
def cog_list(src, S, W, N, E):
    from shapely.geometry import box
    ext = box(W, S, E, N)
    if src.lower().endswith(('.shp', '.zip')):
        import geopandas as gpd
        g = gpd.read_file(src).to_crs(4326)
        g = g[g.intersects(ext)]
        col = next(c for c in g.columns if c.lower() == 'url')
        return [(u, geom.bounds) for u, geom in zip(g[col], g.geometry)]
    urls = [u.strip() for u in (http(src).decode() if src.startswith('http') else open(src).read()).split()
            if u.strip().lower().endswith('.tif')]

    def head(u):
        with rasterio.open('/vsicurl/' + u) as ds:
            return u, transform_bounds(ds.crs, CRS.from_epsg(4326), *ds.bounds)
    with ThreadPoolExecutor(max_workers=8) as ex:
        out = list(ex.map(head, urls))
    return [(u, b) for u, b in out if box(*b).intersects(ext)]


def build_cogs(name, src, mpp, S, W, N, E):
    os.environ.setdefault('GDAL_DISABLE_READDIR_ON_OPEN', 'EMPTY_DIR')
    os.environ.setdefault('CPL_VSIL_CURL_ALLOWED_EXTENSIONS', '.tif')
    os.environ.setdefault('GDAL_HTTP_MAX_RETRY', '5')
    os.environ.setdefault('GDAL_HTTP_RETRY_DELAY', '3')
    items = cog_list(src, S, W, N, E)
    if not items:
        raise SystemExit("no tiles over the extent")
    Wp, Hp, tr = grid(S, W, N, E, mpp)
    print(f"  {Wp} x {Hp} @ {mpp} m/px from {len(items)} tiles", flush=True)
    dst = out_dataset(name, Wp, Hp, tr)
    dlon = (E - W) / Wp; dlat = (N - S) / Hp
    cache = {}
    t0 = time.time()

    def load(u):
        """The tile at the overview nearest mpp: (array, transform, crs)."""
        if u in cache:
            return cache[u]
        for k in range(4):
            try:
                with rasterio.open('/vsicurl/' + u) as ds:
                    res = abs(ds.transform.a)
                    f = max(1, int(mpp / res * 1.0001))
                    oh, ow = ds.height // f, ds.width // f
                    a = ds.read([1, 2, 3], out_shape=(3, oh, ow), resampling=Resampling.average)
                    if ds.count >= 5 or any(ci.name == 'alpha' for ci in ds.colorinterp):
                        m = ds.read(ds.count, out_shape=(oh, ow), resampling=Resampling.nearest)
                        a[:, m == 0] = 0
                    t = ds.transform * rasterio.Affine.scale(ds.width / ow, ds.height / oh)
                    # the raster's own extent: a tile index's footprint can be a
                    # few tens of metres short of it, which once cut a seam
                    # between every pair of tiles
                    bb = transform_bounds(ds.crs, CRS.from_epsg(4326), *ds.bounds)
                    cache[u] = (a, t, ds.crs, bb)
                    return cache[u]
            except Exception as exc:
                if k == 3:
                    print(f"  FAILED {u}: {exc}", flush=True)
                    cache[u] = None
                    return None
                time.sleep(5 * (k + 1))

    nbytes = [0]
    for j in range(0, Hp, CHUNK):
        ch = min(CHUNK, Hp - j)
        bN = N - dlat * j; bS = N - dlat * (j + ch)
        mine = [it for it in items if it[1][3] + 0.002 > bS and it[1][1] - 0.002 < bN]
        # drop tiles no longer needed so the cache holds about two bands
        for u in [u for u in cache if u not in {m[0] for m in mine}]:
            del cache[u]
        out = np.zeros((3, ch, Wp), np.uint8)
        with ThreadPoolExecutor(max_workers=int(arg('--workers', 8))) as ex:
            loaded = list(ex.map(lambda it: (it, load(it[0])), mine))

        def one(x):
            _, got = x
            if got is None:
                return None
            a, t, crs, (w_, s_, e_, n_) = got
            if n_ <= bS or s_ >= bN:
                return None
            i0 = max(0, int((w_ - W) / dlon) - 2); i1 = min(Wp, int(math.ceil((e_ - W) / dlon)) + 2)
            if i1 <= i0:
                return None
            tmp = np.zeros((3, ch, i1 - i0), np.uint8)
            reproject(a, tmp, src_transform=t, src_crs=crs,
                      dst_transform=from_origin(W + i0 * dlon, bN, dlon, dlat), dst_crs=CRS.from_epsg(4326),
                      resampling=Resampling.bilinear, src_nodata=0, dst_nodata=0)
            return i0, tmp
        with ThreadPoolExecutor(max_workers=4) as ex:
            for r in ex.map(one, loaded):
                if r is None:
                    continue
                i0, tmp = r
                m = tmp.max(axis=0) > 0
                sub = out[:, :, i0:i0 + tmp.shape[2]]
                sub[:, m] = tmp[:, m]
        dst.write(out, window=rasterio.windows.Window(0, j, Wp, ch))
        print(f"  rows {j+ch}/{Hp}, {len(mine)} tiles {time.time()-t0:.0f}s", flush=True)
    finish(dst, name, S, W, N, E, Wp, Hp, mpp, dict(kind='cogs', url=src, tiles=len(items)))


def clip_to(ext, cover):
    """--clip: store only where the source has coverage, not the whole extent."""
    if not cover:
        return ext
    S, W, N, E = ext
    S2, W2, N2, E2 = max(S, cover[0]), max(W, cover[1]), min(N, cover[2]), min(E, cover[3])
    if S2 >= N2 or W2 >= E2:
        raise SystemExit(f"source coverage {cover} misses the extent {ext}")
    return [S2, W2, N2, E2]


def main():
    name = sys.argv[1]
    S, W, N, E = target_extent()
    mpp = float(arg('--mpp', 1.0))
    if '--clip' in sys.argv:
        cover = ([float(x) for x in arg('--cover').split(',')] if '--cover' in sys.argv else
                 arcgis_info(arg('--arcgis').rstrip('/'))[0] if '--arcgis' in sys.argv else
                 tile_cover(arg('--tiles')) if '--tiles' in sys.argv else None)
        S, W, N, E = clip_to([S, W, N, E], cover)
    print(f"layer {name}: {S:.4f}..{N:.4f} N, {W:.4f}..{E:.4f} W", flush=True)
    if '--arcgis' in sys.argv:
        build_arcgis(name, arg('--arcgis').rstrip('/'), mpp, S, W, N, E)
    elif '--naip' in sys.argv:
        build_naip(name, arg('--naip'), mpp, S, W, N, E)
    elif '--tiles' in sys.argv:
        build_tiles(name, arg('--tiles'), mpp, S, W, N, E)
    elif '--cogs' in sys.argv:
        build_cogs(name, arg('--cogs'), mpp, S, W, N, E)
    else:
        raise SystemExit('--arcgis URL, --naip YEAR, --tiles TEMPLATE or --cogs INDEX')


if __name__ == '__main__':
    main()
