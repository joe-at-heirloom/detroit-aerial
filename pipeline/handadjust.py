"""A hand alignment: one similarity transform per layer, about the layer's centre.

Kept in `data/adjust.json` and applied at read time -- by the tile server when it
cuts tiles, and by the measurement scripts when they sample a layer -- rather than
baked into the raster. So it is instant, reversible, and four readable numbers.

The geometry lives here and nowhere else. This project's recurring failure is two
implementations of the same transform drifting apart until a measurement quietly
answers a different question from the thing being served, so the viewer, the tile
server and `dtcross.py` all come through this file.

Sign convention: `deg` is CLOCKWISE ON SCREEN -- what a person turning a handle
means, and the direction CSS `rotate()` turns. In the E/N frame that is a negative
mathematical angle, which is the minus in `source_lonlat`.
"""
import json, math, os

IDENTITY = dict(dE=0.0, dN=0.0, deg=0.0, scale=1.0)

# the project's one metric frame (pipeline/dtmap.py holds the same constants)
LAT0, LON0 = 42.3340, -83.0450
_phi = math.radians(LAT0)
MLAT = 111132.92 - 559.82 * math.cos(2 * _phi) + 1.175 * math.cos(4 * _phi)
MLON = 111412.84 * math.cos(_phi) - 93.5 * math.cos(3 * _phi)


def path(root):
    return os.path.join(root, "data", "adjust.json")


def load(root):
    try:
        return json.load(open(path(root))).get("layers", {})
    except Exception:
        return {}


def for_layer(layers, layer_id):
    """The layer's transform, or None when it is the identity -- so callers can
    keep their fast path and only pay for the general one when it is asked for."""
    a = (layers or {}).get(layer_id)
    if not a:
        return None
    a = {k: float(a.get(k, IDENTITY[k])) for k in IDENTITY}
    if a["scale"] <= 0:
        a["scale"] = 1.0
    if (abs(a["dE"]) < 1e-6 and abs(a["dN"]) < 1e-6
            and abs(a["deg"]) < 1e-9 and abs(a["scale"] - 1.0) < 1e-9):
        return None
    return a


def source_lonlat(lon, lat, bbox, adj):
    """Where to READ from, for content asked for at (lon, lat).

    Forward, the transform sends a source point p to R*s*(p - c) + c + d; this is
    its inverse, because rendering asks the opposite question. `bbox` is the
    layer's own (S, W, N, E) and gives the pivot c.
    """
    S, W, N, E = bbox
    Ec = ((W + E) / 2.0 - LON0) * MLON
    Nc = ((S + N) / 2.0 - LAT0) * MLAT
    u = (lon - LON0) * MLON - Ec - adj["dE"]
    v = (lat - LAT0) * MLAT - Nc - adj["dN"]
    th = math.radians(adj["deg"])          # inverse of the screen-clockwise rotation
    ct, st = math.cos(th), math.sin(th)
    s = adj["scale"]
    Es = (u * ct - v * st) / s + Ec
    Ns = (u * st + v * ct) / s + Nc
    return LON0 + Es / MLON, LAT0 + Ns / MLAT
