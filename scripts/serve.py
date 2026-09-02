#!/usr/bin/env python3
"""Local XYZ tile server for the Detroit Air Survey.

Reads windows straight out of the source rasters on demand, so there is no tile
pyramid on disk and no resolution ceiling -- you can zoom to the native 0.63 m/px
of the negatives. Sources are EPSG:4326 with a linear lat/lon transform; tiles are
Web Mercator, so each tile is resampled row-by-row (latitude is the only nonlinear
axis at this scale).

    python3 scripts/serve.py [--port 8770]
"""
import argparse, io, json, math, os, sys, threading
from functools import lru_cache
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler

import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = json.load(open(os.path.join(ROOT, "data", "manifest.json")))
LAYERS = {l["id"]: l for l in MANIFEST["layers"]}
TILE = 256

try:
    import rasterio
    HAVE_RIO = True
except ImportError:
    HAVE_RIO = False


# ---------------------------------------------------------------- source access
class Source:
    """Uniform read(window)->ndarray over either a GeoTIFF or a plain image.

    Neither rasterio datasets nor PIL image handles are thread-safe, and this is a
    threading server, so every worker thread gets its own handle."""

    def __init__(self, path, bbox):
        self.path = path
        self.S, self.W, self.N, self.E = bbox
        self.is_tif = path.lower().endswith((".tif", ".tiff"))
        self._tl = threading.local()
        if self.is_tif and not HAVE_RIO:
            raise RuntimeError("rasterio required for %s" % path)
        h = self._handle()
        self.w, self.h = (h.width, h.height) if self.is_tif else h.size

    def _handle(self):
        h = getattr(self._tl, "h", None)
        if h is None:
            h = rasterio.open(self.path) if self.is_tif else Image.open(self.path)
            self._tl.h = h
        return h

    def read(self, px0, py0, px1, py1, out_w, out_h):
        px0 = max(0, min(self.w - 1, px0)); px1 = max(px0 + 1, min(self.w, px1))
        py0 = max(0, min(self.h - 1, py0)); py1 = max(py0 + 1, min(self.h, py1))
        h = self._handle()
        if self.is_tif:
            win = rasterio.windows.Window(px0, py0, px1 - px0, py1 - py0)
            return h.read(1, window=win, out_shape=(out_h, out_w),
                          resampling=rasterio.enums.Resampling.average)
        crop = h.crop((px0, py0, px1, py1))
        if crop.size != (out_w, out_h):
            crop = crop.resize((out_w, out_h), Image.LANCZOS)
        return np.asarray(crop.convert("L") if crop.mode not in ("L", "RGB") else crop)


_sources = {}
_lock = threading.Lock()


def source(layer_id):
    with _lock:
        if layer_id not in _sources:
            l = LAYERS[layer_id]
            _sources[layer_id] = Source(os.path.join(ROOT, l["file"]), l["bbox"])
        return _sources[layer_id]


# ---------------------------------------------------------------- tile geometry
def tile_bounds(z, x, y):
    n = 2 ** z
    lon0 = x / n * 360.0 - 180.0
    lon1 = (x + 1) / n * 360.0 - 180.0
    lat0 = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
    lat1 = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * (y + 1) / n))))
    return lat1, lon0, lat0, lon1          # S, W, N, E


def render_tile(layer_id, z, x, y):
    l = LAYERS[layer_id]
    src = source(layer_id)
    tS, tW, tN, tE = tile_bounds(z, x, y)
    if tE <= src.W or tW >= src.E or tN <= src.S or tS >= src.N:
        return None                        # tile lies outside the mosaic

    # longitude is linear in both frames -> a single column mapping
    fx0 = (tW - src.W) / (src.E - src.W) * src.w
    fx1 = (tE - src.W) / (src.E - src.W) * src.w
    px0, px1 = int(math.floor(fx0)), int(math.ceil(fx1))

    # latitude: tile rows are Mercator, source rows are linear in degrees
    def merc(lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    mN, mS = merc(tN), merc(tS)
    rows = np.arange(TILE) + 0.5
    lat = np.array([math.degrees(2 * math.atan(math.exp(mN + (mS - mN) * r / TILE)) - math.pi / 2)
                    for r in rows])
    fy = (src.N - lat) / (src.N - src.S) * src.h
    py0, py1 = int(math.floor(fy.min())), int(math.ceil(fy.max()))
    if px1 <= px0 or py1 <= py0:
        return None

    band = src.read(px0, py0, px1, py1, TILE, max(1, py1 - py0))
    if band is None or band.size == 0:
        return None
    if band.ndim == 3:
        band = band[..., :3]
    # resample rows onto the Mercator grid
    idx = np.clip(((fy - py0) / max(1, (py1 - py0)) * band.shape[0]).astype(int),
                  0, band.shape[0] - 1)
    out = band[idx]
    if out.ndim == 2:
        rgb = np.dstack([out, out, out])
    else:
        rgb = out
    alpha = np.where(rgb.max(axis=2) > 0, 255, 0).astype(np.uint8)
    rgba = np.dstack([rgb.astype(np.uint8), alpha])
    im = Image.fromarray(rgba, "RGBA")
    if im.size != (TILE, TILE):
        im = im.resize((TILE, TILE), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=False)
    return buf.getvalue()


# ---------------------------------------------------------------- road overlay
_ROADS = None


def roads():
    """City of Detroit centrelines (+ OSM beyond the city limits), in lat/lon."""
    global _ROADS
    if _ROADS is None:
        p = os.path.join(ROOT, "data", "overlay_segs.npy")
        seg = np.load(p)                       # metres in the dtmap frame
        LAT0, LON0 = 42.3340, -83.0450
        phi = math.radians(LAT0)
        MLAT = 111132.92 - 559.82 * math.cos(2 * phi) + 1.175 * math.cos(4 * phi)
        MLON = 111412.84 * math.cos(phi) - 93.5 * math.cos(3 * phi)
        _ROADS = np.column_stack([
            LON0 + seg[:, 0] / MLON, LAT0 + seg[:, 1] / MLAT,
            LON0 + seg[:, 2] / MLON, LAT0 + seg[:, 3] / MLAT]).astype(np.float64)
    return _ROADS


def render_roads(z, x, y):
    from PIL import ImageDraw
    if z < 14:                              # whole-city view: too many segments to be useful
        return None
    S, W, N, E = tile_bounds(z, x, y)
    pad = (E - W) * 0.2
    R = roads()
    m = ~(((R[:, 0] < W - pad) & (R[:, 2] < W - pad)) |
          ((R[:, 0] > E + pad) & (R[:, 2] > E + pad)) |
          ((R[:, 1] < S - pad) & (R[:, 3] < S - pad)) |
          ((R[:, 1] > N + pad) & (R[:, 3] > N + pad)))
    sel = R[m]
    im = Image.new("RGBA", (TILE, TILE), (0, 0, 0, 0))
    if len(sel) == 0:
        buf = io.BytesIO(); im.save(buf, "PNG"); return buf.getvalue()

    def merc(lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    mN, mS = merc(N), merc(S)
    d = ImageDraw.Draw(im)
    X0 = (sel[:, 0] - W) / (E - W) * TILE
    X1 = (sel[:, 2] - W) / (E - W) * TILE
    Y0 = (mN - np.array([merc(v) for v in sel[:, 1]])) / (mN - mS) * TILE
    Y1 = (mN - np.array([merc(v) for v in sel[:, 3]])) / (mN - mS) * TILE
    wdt = 1 if z < 16 else 2
    for a, b, c, e in zip(X0, Y0, X1, Y1):
        d.line([(a, b), (c, e)], fill=(0, 229, 255, 235), width=wdt)
    buf = io.BytesIO()
    im.save(buf, "PNG")
    return buf.getvalue()


BLANK = None


def blank():
    global BLANK
    if BLANK is None:
        buf = io.BytesIO()
        Image.new("RGBA", (TILE, TILE), (0, 0, 0, 0)).save(buf, "PNG")
        BLANK = buf.getvalue()
    return BLANK


# ---------------------------------------------------------------- http
class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=ROOT, **kw)

    def log_message(self, *a):
        pass

    def do_GET(self):
        if self.path.startswith("/t/"):
            try:
                _, _, rest = self.path.partition("/t/")
                lid, z, x, yy = rest.split("/")
                y = int(yy.split(".")[0])
                if lid == "__roads__":
                    png = render_roads(int(z), int(x), y)
                else:
                    png = render_tile(lid, int(z), int(x), y)
            except Exception as exc:
                sys.stderr.write("tile error %s: %s\n" % (self.path, exc))
                png = None
            data = png or blank()
            self.send_response(200)
            self.send_header("Content-Type", "image/png")
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "max-age=600")
            self.end_headers()
            self.wfile.write(data)
            return
        if self.path in ("/", "/index.html"):
            self.path = "/viewer/index.html"
        return super().do_GET()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8770)
    a = ap.parse_args()
    print("Detroit Air Survey")
    print("  root   %s" % ROOT)
    print("  layers %s" % ", ".join(LAYERS))
    print("  http://localhost:%d" % a.port)
    ThreadingHTTPServer(("127.0.0.1", a.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
