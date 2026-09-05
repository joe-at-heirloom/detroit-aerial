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
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
import handadjust
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
        # the film mosaics are one band; the modern orthos (1998 on) are three
        self.bands = [1, 2, 3] if self.is_tif and h.count >= 3 else 1

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
            if self.bands == 1:
                return h.read(1, window=win, out_shape=(out_h, out_w),
                              resampling=rasterio.enums.Resampling.average)
            a = h.read(self.bands, window=win, out_shape=(3, out_h, out_w),
                       resampling=rasterio.enums.Resampling.average)
            return np.ascontiguousarray(np.transpose(a, (1, 2, 0)))
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

    # A tile that hangs over the edge of the mosaic may only paint the part of
    # itself that lands on real ground. Clamping the read window WITHOUT also
    # narrowing the output rectangle -- which this did once -- resampled the few
    # columns that exist across the full 256, so the outermost strip of every
    # block was stretched and displaced by up to a whole tile. So: work out which
    # columns and rows of the OUTPUT are drawable first, and paint only those.

    # longitude is linear in both frames -> a single column mapping
    lo0 = max(tW, src.W); lo1 = min(tE, src.E)
    if lo1 <= lo0:
        return None
    oc0 = max(0, int(math.floor((lo0 - tW) / (tE - tW) * TILE)))
    oc1 = min(TILE, int(math.ceil((lo1 - tW) / (tE - tW) * TILE)))
    if oc1 <= oc0:
        return None
    lo0e = tW + (tE - tW) * oc0 / TILE; lo1e = tW + (tE - tW) * oc1 / TILE
    px0 = max(0, int(math.floor((lo0e - src.W) / (src.E - src.W) * src.w)))
    px1 = min(src.w, int(math.ceil((lo1e - src.W) / (src.E - src.W) * src.w)))

    # latitude: tile rows are Mercator, source rows are linear in degrees
    def merc(lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    mN, mS = merc(tN), merc(tS)
    rows = np.arange(TILE) + 0.5
    lat = np.array([math.degrees(2 * math.atan(math.exp(mN + (mS - mN) * r / TILE)) - math.pi / 2)
                    for r in rows])
    fy = (src.N - lat) / (src.N - src.S) * src.h
    okr = (fy >= 0) & (fy < src.h)
    if not okr.any():
        return None
    or0 = int(np.argmax(okr)); or1 = TILE - int(np.argmax(okr[::-1]))
    py0 = max(0, int(math.floor(fy[or0])))
    py1 = min(src.h, int(math.ceil(fy[or1 - 1])))
    if px1 <= px0 or py1 <= py0:
        return None

    band = src.read(px0, py0, px1, py1, oc1 - oc0, max(1, py1 - py0))
    if band is None or band.size == 0:
        return None
    if band.ndim == 3:
        band = band[..., :3]
    # resample rows onto the Mercator grid, over the drawable rows only
    idx = np.clip(((fy[or0:or1] - py0) / max(1, (py1 - py0)) * band.shape[0]).astype(int),
                  0, band.shape[0] - 1)
    sub = band[idx]
    if sub.ndim == 2:
        sub = np.dstack([sub, sub, sub])
    sub = sub[:, :oc1 - oc0, :3]

    rgb = np.zeros((TILE, TILE, 3), np.uint8)
    rgb[or0:or0 + sub.shape[0], oc0:oc0 + sub.shape[1]] = sub
    alpha = np.zeros((TILE, TILE), np.uint8)
    alpha[or0:or0 + sub.shape[0], oc0:oc0 + sub.shape[1]] = np.where(sub.max(axis=2) > 0, 255, 0)
    rgba = np.dstack([rgb, alpha])
    im = Image.fromarray(rgba, "RGBA")
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=False)
    return buf.getvalue()


# ---------------------------------------------------------------- adjustments
# A hand alignment: one similarity transform per layer, in ground metres about the
# layer's own centre. It is applied here, at serve time, rather than baked into the
# raster -- so it is instant, reversible, and readable as four numbers. Baking it
# into the GeoTIFF is a separate, deliberate step (scripts/bake.py), which is what
# the offline measurement tools read.
#
# Sign convention, fixed here and matched by the viewer: `deg` is CLOCKWISE on
# screen, which is what a person dragging a rotate handle means, and the same
# direction CSS `rotate()` turns. In the E/N frame that is a negative mathematical
# angle, which is where the minus below comes from.
ADJUST_PATH = handadjust.path(ROOT)
IDENTITY = handadjust.IDENTITY
LAT0, LON0, MLAT, MLON = (handadjust.LAT0, handadjust.LON0,
                          handadjust.MLAT, handadjust.MLON)
_adj_cache = {"mtime": None, "val": {}}


def adjustments():
    """Reload when the file changes, so a save takes effect without a restart."""
    try:
        m = os.path.getmtime(ADJUST_PATH)
    except OSError:
        _adj_cache.update(mtime=None, val={})
        return {}
    if m != _adj_cache["mtime"]:
        _adj_cache.update(mtime=m, val=handadjust.load(ROOT))
    return _adj_cache["val"]


def adjust_for(layer_id):
    return handadjust.for_layer(adjustments(), layer_id)


def render_tile_adj(layer_id, z, x, y, adj):
    """Tile render through a similarity transform.

    The unadjusted path is separable -- longitude is linear in both frames, so it
    maps columns, and latitude maps rows -- which is why it is fast. Rotation
    couples the axes, so this one builds the full 2D inverse map instead. It is
    256x256 of numpy either way and the cost does not show."""
    l = LAYERS[layer_id]
    src = source(layer_id)
    tS, tW, tN, tE = tile_bounds(z, x, y)

    lon = tW + (np.arange(TILE) + 0.5) * (tE - tW) / TILE

    def merc(lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    mN, mS = merc(tN), merc(tS)
    rows = (np.arange(TILE) + 0.5) / TILE
    lat = np.array([math.degrees(2 * math.atan(math.exp(mN + (mS - mN) * r)) - math.pi / 2)
                    for r in rows])
    LON, LAT = np.meshgrid(lon, lat)
    lon_s, lat_s = handadjust.source_lonlat(LON, LAT, (src.S, src.W, src.N, src.E), adj)

    fx = (lon_s - src.W) / (src.E - src.W) * src.w
    fy = (src.N - lat_s) / (src.N - src.S) * src.h
    inside = (fx >= 0) & (fx < src.w) & (fy >= 0) & (fy < src.h)
    if not inside.any():
        return None

    px0 = max(0, int(math.floor(fx[inside].min()))); px1 = min(src.w, int(math.ceil(fx[inside].max())) + 1)
    py0 = max(0, int(math.floor(fy[inside].min()))); py1 = min(src.h, int(math.ceil(fy[inside].max())) + 1)
    if px1 <= px0 or py1 <= py0:
        return None
    # never pull more than a few tiles' worth of source, however zoomed out we are
    ow = max(1, min(px1 - px0, TILE * 3)); oh = max(1, min(py1 - py0, TILE * 3))
    band = src.read(px0, py0, px1, py1, ow, oh)
    if band is None or band.size == 0:
        return None
    if band.ndim == 3:
        band = band[..., :3]
    if band.ndim == 2:
        band = np.dstack([band, band, band])

    ix = np.clip(((fx - px0) / (px1 - px0) * band.shape[1]).astype(np.int32), 0, band.shape[1] - 1)
    iy = np.clip(((fy - py0) / (py1 - py0) * band.shape[0]).astype(np.int32), 0, band.shape[0] - 1)
    rgb = band[iy, ix][..., :3]
    rgb = np.where(inside[..., None], rgb, 0).astype(np.uint8)
    alpha = np.where(inside & (rgb.max(axis=2) > 0), 255, 0).astype(np.uint8)
    im = Image.fromarray(np.dstack([rgb, alpha]), "RGBA")
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


# ---------------------------------------------------------------- places
# Pins are researched by hand and the blurbs are written by hand, but the
# coordinate should not be: the viewer knows exactly which ground the probe is
# open on, so it posts here and this appends the stub. The file is the one piece
# of hand-written work in the repo that is not regenerable, so it is written to a
# temporary file and renamed -- a half-written places.json would lose all of it.
PLACES_PATH = os.path.join(ROOT, "data", "places.json")
KINDS = {"housing", "school", "factory", "retail", "office", "civic", "park",
         "transport", "street", "leisure", "hospital", "worship", "farm",
         "landform", "other"}
TO_KINDS = {"demolished", "cleared", "closed", "rebuilt", "filled", "channelised"}
FIELDS = ["id", "lat", "lon", "title", "from", "to", "to_kind", "kind", "rank",
          "blurb", "sources", "group"]


_places_lock = threading.Lock()
_adjust_lock = threading.Lock()


def slug(title, taken):
    base = "".join(c if c.isalnum() else "-" for c in title.lower()).strip("-")
    while "--" in base:
        base = base.replace("--", "-")
    base = base[:44].strip("-") or "place"
    out, n = base, 2
    while out in taken:
        out, n = "%s-%d" % (base, n), n + 1
    return out


def year(v):
    """A year is an int, or a string the viewer can still sort ("1950s"); nothing else."""
    if v in (None, ""):
        return None
    if isinstance(v, bool):
        raise ValueError("year")
    if isinstance(v, int):
        return v if 1700 <= v <= 2100 else None
    v = str(v).strip()[:12]
    return v if any(c.isdigit() for c in v) else None


def set_adjust(req):
    """Record a hand alignment. Written whole, so a save is one atomic swap."""
    lid = str(req.get("layer", ""))
    if lid not in LAYERS:
        raise ValueError("unknown layer %r" % lid)
    a = {}
    for k, lo, hi in (("dE", -5000.0, 5000.0), ("dN", -5000.0, 5000.0),
                      ("deg", -180.0, 180.0), ("scale", 0.5, 2.0)):
        v = float(req.get(k, IDENTITY[k]))
        if not (lo <= v <= hi):
            raise ValueError("%s out of range: %g" % (k, v))
        a[k] = round(v, 4)
    with _adjust_lock:
        doc = {"layers": {}}
        if os.path.exists(ADJUST_PATH):
            try:
                doc = json.load(open(ADJUST_PATH))
                doc.setdefault("layers", {})
            except Exception:
                pass
        if a == IDENTITY:
            doc["layers"].pop(lid, None)      # identity is absence, not a row of zeros
        else:
            doc["layers"][lid] = a
        tmp = ADJUST_PATH + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(doc, fh, indent=1, sort_keys=True)
            fh.write("\n")
        os.replace(tmp, ADJUST_PATH)
    return {"layer": lid, "adjust": a}


def add_place(req):
    lat, lon = float(req["lat"]), float(req["lon"])
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        raise ValueError("coordinate out of range")
    title = str(req.get("title", "")).strip()
    if not title:
        raise ValueError("title is required")
    group = str(req.get("group", "")).strip()
    if group not in MANIFEST["groups"]:
        raise ValueError("unknown group %r" % group)
    kind = str(req.get("kind") or "other")
    if kind not in KINDS:
        raise ValueError("unknown kind %r" % kind)
    to_kind = req.get("to_kind") or None
    if to_kind is not None and to_kind not in TO_KINDS:
        raise ValueError("unknown to_kind %r" % to_kind)
    rank = int(req.get("rank") or 2)
    if rank not in (1, 2, 3):
        raise ValueError("rank must be 1, 2 or 3")
    note = str(req.get("note", "")).strip()

    with _places_lock:
        doc = json.load(open(PLACES_PATH))
        place = dict(id=slug(title, {p["id"] for p in doc["places"]}),
                     lat=round(lat, 5), lon=round(lon, 5), title=title,
                     to_kind=to_kind, kind=kind, rank=rank,
                     blurb=note, sources=[], group=group)
        place["from"] = year(req.get("from"))
        place["to"] = year(req.get("to"))
        doc["places"].append({k: place[k] for k in FIELDS})
        tmp = PLACES_PATH + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(doc, fh, indent=1, ensure_ascii=False)
            fh.write("\n")
        os.replace(tmp, PLACES_PATH)
    return place


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
                    adj = adjust_for(lid)
                    png = (render_tile_adj(lid, int(z), int(x), y, adj) if adj
                           else render_tile(lid, int(z), int(x), y))
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
        if self.path.split("?")[0] == "/api/adjust":
            self._json({"layers": adjustments()})
            return
        if self.path in ("/", "/index.html"):
            self.path = "/viewer/index.html"
        return super().do_GET()

    def _json(self, obj, code=200):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        # The server binds 127.0.0.1, so this is reachable only from this machine.
        handler = {"/api/place": add_place, "/api/adjust": set_adjust}.get(self.path)
        if handler is None:
            self.send_error(404)
            return
        try:
            n = int(self.headers.get("Content-Length") or 0)
            if not 0 < n <= 65536:
                raise ValueError("bad body length")
            self._json(handler(json.loads(self.rfile.read(n))))
        except Exception as exc:
            self._json({"error": str(exc)}, 400)


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
