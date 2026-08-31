"""Resample a mosaic through a locally-weighted displacement field.

The field is evaluated on a coarse lattice (the RBF has a ~1 km length scale, so
a 200 m lattice is far finer than the signal) and bilinearly interpolated per
pixel. Sampling is inverse: output(p) = source(p + d(p))."""
import numpy as np, math, os, json
import rasterio
from rasterio.transform import from_origin

LAT0, LON0 = 42.3340, -83.0450
_p = math.radians(LAT0)
MLAT = 111132.92 - 559.82 * math.cos(2 * _p) + 1.175 * math.cos(4 * _p)
MLON = 111412.84 * math.cos(_p) - 93.5 * math.cos(3 * _p)


def field(warp_at, minE, maxN, W, H, mpp, step_m=200.0):
    """Sample the displacement field onto a lattice covering the raster."""
    nx = max(2, int(W * mpp / step_m) + 2)
    ny = max(2, int(H * mpp / step_m) + 2)
    gx = np.linspace(0, W * mpp, nx)
    gy = np.linspace(0, H * mpp, ny)
    dE = np.zeros((ny, nx), np.float32)
    dN = np.zeros((ny, nx), np.float32)
    for j, oy in enumerate(gy):
        for i, ox in enumerate(gx):
            e = minE + ox
            n = maxN - oy
            d = warp_at(e, n)
            dE[j, i] = d[0]
            dN[j, i] = d[1]
    return gx, gy, dE, dN


def _bilerp(grid, gx, gy, xs, ys):
    fx = np.clip(np.interp(xs, gx, np.arange(len(gx))), 0, len(gx) - 1.001)
    fy = np.clip(np.interp(ys, gy, np.arange(len(gy))), 0, len(gy) - 1.001)
    x0 = fx.astype(np.int32); y0 = fy.astype(np.int32)
    tx = (fx - x0)[None, :]; ty = (fy - y0)[:, None]
    g = grid
    return (g[y0][:, x0] * (1 - tx) * (1 - ty) + g[y0][:, x0 + 1] * tx * (1 - ty)
            + g[y0 + 1][:, x0] * (1 - tx) * ty + g[y0 + 1][:, x0 + 1] * tx * ty)


def apply(src_tif, geo, warp_at, out_tif, raw_path, chunk=1024, log=print):
    g = json.load(open(geo)) if isinstance(geo, str) else geo
    minE, maxN, W, H, mpp = g['minE'], g['maxN'], g['W'], g['H'], g['mpp']
    ds = rasterio.open(src_tif)
    gx, gy, DE, DN = field(warp_at, minE, maxN, W, H, mpp)
    log(f"  displacement lattice {DE.shape[1]}x{DE.shape[0]}  "
        f"dE {DE.min():+.0f}..{DE.max():+.0f}  dN {DN.min():+.0f}..{DN.max():+.0f} m")
    out = np.memmap(raw_path, dtype=np.uint8, mode='w+', shape=(H, W)); out[:] = 0
    xs_m = (np.arange(W) + 0.5) * mpp
    for y0 in range(0, H, chunk):
        y1 = min(H, y0 + chunk)
        ys_m = (np.arange(y0, y1) + 0.5) * mpp
        de = _bilerp(DE, gx, gy, xs_m, ys_m)
        dn = _bilerp(DN, gx, gy, xs_m, ys_m)
        # inverse sample: output at p reads source at p + d(p)
        sx = (xs_m[None, :] + de) / mpp
        sy = (ys_m[:, None] - dn) / mpp
        x0 = int(max(0, np.floor(sx.min()) - 2)); x1 = int(min(W, np.ceil(sx.max()) + 2))
        r0 = int(max(0, np.floor(sy.min()) - 2)); r1 = int(min(H, np.ceil(sy.max()) + 2))
        if x1 <= x0 or r1 <= r0:
            continue
        band = ds.read(1, window=rasterio.windows.Window(x0, r0, x1 - x0, r1 - r0))
        ix = np.clip(sx - x0, 0, band.shape[1] - 1.001)
        iy = np.clip(sy - r0, 0, band.shape[0] - 1.001)
        xa = ix.astype(np.int32); ya = iy.astype(np.int32)
        fx = (ix - xa).astype(np.float32); fy = (iy - ya).astype(np.float32)
        v = (band[ya, xa] * (1 - fx) * (1 - fy) + band[ya, xa + 1] * fx * (1 - fy)
             + band[ya + 1, xa] * (1 - fx) * fy + band[ya + 1, xa + 1] * fx * fy)
        inb = ((sx >= 0) & (sx < W - 1) & (sy >= 0) & (sy < H - 1))
        out[y0:y1] = np.where(inb, v, 0).astype(np.uint8)
        if (y0 // chunk) % 12 == 0:
            log(f"    row {y0}/{H}")
    out.flush()
    la_top = LAT0 + maxN / MLAT
    lo_left = LON0 + minE / MLON
    tr = from_origin(lo_left, la_top, mpp / MLON, mpp / MLAT)
    with rasterio.open(out_tif, 'w', driver='GTiff', height=H, width=W, count=1,
                       dtype='uint8', crs='EPSG:4326', transform=tr,
                       tiled=True, blockxsize=512, blockysize=512,
                       compress='DEFLATE', predictor=2, num_threads='ALL_CPUS',
                       BIGTIFF='YES') as dst:
        for y in range(0, H, 2048):
            yb = min(H, y + 2048)
            dst.write(np.asarray(out[y:yb]), 1,
                      window=rasterio.windows.Window(0, y, W, yb - y))
        dst.build_overviews([2, 4, 8, 16, 32, 64], rasterio.enums.Resampling.average)
    del out
    try: os.remove(raw_path)
    except OSError: pass
    return dict(bbox=[LAT0 + (maxN - H * mpp) / MLAT, lo_left, la_top,
                      LON0 + (minE + W * mpp) / MLON], W=W, H=H, mpp=mpp,
                minE=minE, maxN=maxN)
