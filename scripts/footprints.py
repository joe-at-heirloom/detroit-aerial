#!/usr/bin/env python
"""Where each layer actually has ground, and what the year view shows per year.

A layer's bbox is a rectangle; its coverage is not. The film blocks are ragged,
the City of Detroit's orthos stop at the city limits, Oakland County's stop at
8 Mile. The year view (viewer, "Through the years") draws the newest layer at or
before the chosen year on top and lets older layers show wherever they reach
further, each outlined thinly and captioned with its year -- so it needs the real
outlines, and for each year, which part of which layer is left visible.

Footprints come from the rasters themselves: a coarse read (~15 m cells) of
where any band is non-zero, which is exactly what serve.py draws as opaque,
cleaned of specks and pinholes and vectorised. All geometry is done in ground
metres about dtmap's origin and written back as lon/lat.

    ./.venv/bin/python scripts/footprints.py          # after manifest.py
writes data/timeline.json
"""
import os, sys, json, math, time
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_bounds
from rasterio import features
from scipy import ndimage
from shapely.geometry import shape, mapping, Polygon, MultiPolygon, Point
from shapely.ops import unary_union
from shapely import affinity
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import dtmap
def P(*a): return os.path.join(ROOT, *a)

CELL = 15.0          # metres per mask cell
MIN_PART = 0.05e6    # m^2: smaller islands of coverage are specks, not coverage
MIN_HOLE = 0.25e6    # m^2: smaller holes are dark pixels inside imagery
SIMPLIFY = 12.0      # m
SLIVER = 40.0        # m: a visible remainder thinner than 2x this is an edge artefact
TODAY = 9999


def to_m(g):
    return affinity.affine_transform(g, [dtmap.MLON, 0, 0, dtmap.MLAT,
                                         -dtmap.LON0 * dtmap.MLON, -dtmap.LAT0 * dtmap.MLAT])


def to_ll(g):
    return affinity.affine_transform(g, [1 / dtmap.MLON, 0, 0, 1 / dtmap.MLAT,
                                         dtmap.LON0, dtmap.LAT0])


def year_of(l):
    if 'year' in l:
        return int(l['year'])
    try:
        return int(str(l['label'])[:4])
    except ValueError:
        return TODAY                        # "Today"


HOSTS = [('egis.detroitmi.gov', 'City of Detroit'), ('qvkbeam7Wirps6zC', 'City of Detroit'),
         ('imagery.michigan.gov', 'State of Michigan'), ('nationalmap.gov', 'USGS'),
         ('oakgov.com', 'Oakland County'), ('planetarycomputer', 'USDA NAIP'),
         ('wayback', 'Esri Wayback'), ('RPhrOu9XQzI31xTa', 'Michigan Tech')]


def credit(l):
    """Who made the picture, in a few words, for the year view's readout."""
    if l.get('credit'):
        return l['credit']
    if l['id'].startswith(('b', 'dt')) and l['id'] != 'dtnow':
        return 'DTE film'
    if l['id'] in ('modwest', 'dtnow'):
        return 'Esri World Imagery'
    g = P('data', f"layer_{l['id'][1:]}_geo.json")
    if os.path.exists(g):
        d = json.load(open(g))
        if d.get('credit'):
            return d['credit']
        url = str(d.get('source', {}).get('url', ''))
        if 'NAIP' in url or d.get('source', {}).get('kind') == 'naip':
            return 'USDA NAIP'
        if 'NAPP' in url:
            return 'USGS NAPP'
        for k, v in HOSTS:
            if k in url:
                return v
    return ''


def coverage_mask(l):
    S, W, N, E = l['bbox']
    w = max(8, int(round((E - W) * dtmap.MLON / CELL)))
    h = max(8, int(round((N - S) * dtmap.MLAT / CELL)))
    f = P(l['file'])
    if f.lower().endswith(('.tif', '.tiff')):
        with rasterio.open(f) as ds:
            a = ds.read(list(range(1, min(3, ds.count) + 1)), out_shape=(min(3, ds.count), h, w),
                        resampling=Resampling.nearest)
        m = a.max(axis=0) > 0
    else:
        im = Image.open(f)
        im.draft('RGB', (w, h))
        a = np.asarray(im.convert('L').resize((w, h), Image.NEAREST))
        m = a > 0
    return m, from_bounds(W, S, E, N, w, h)


def clean(m):
    m = ndimage.binary_opening(m, iterations=2)
    m = ndimage.binary_closing(m, iterations=3)
    # fill holes below MIN_HOLE: label the background, refill small enclosed parts
    bg, n = ndimage.label(~m)
    if n:
        sizes = ndimage.sum(np.ones_like(bg), bg, range(1, n + 1)) * CELL * CELL
        edge = set(np.unique(np.concatenate([bg[0], bg[-1], bg[:, 0], bg[:, -1]]))) - {0}
        small = [i + 1 for i, s in enumerate(sizes) if s < MIN_HOLE and (i + 1) not in edge]
        if small:
            m = m | np.isin(bg, small)
    return m


def footprint(l):
    if l.get('url'):
        # A remote layer has no raster here: its outline is the boxes it is known
        # to cover, snapped outward to the edges of the tiles it is drawn from --
        # Leaflet draws whole tiles, so that is the ground it actually paints.
        from shapely.geometry import box
        z = l.get('minNativeZoom', 16); n = 2 ** z
        def snap(c):
            S, W, N, E = c
            x0 = math.floor((W + 180) / 360 * n); x1 = math.ceil((E + 180) / 360 * n)
            ty = lambda lat: (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n
            y0 = math.floor(ty(N)); y1 = math.ceil(ty(S))
            lat = lambda y: math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
            return box(x0 / n * 360 - 180, lat(y1), x1 / n * 360 - 180, lat(y0))
        return to_m(unary_union([snap(c) for c in l.get('cover') or [l['bbox']]]))
    m, tr = coverage_mask(l)
    m = clean(m)
    polys = [shape(g) for g, v in features.shapes(m.astype(np.uint8), mask=m, transform=tr) if v == 1]
    if not polys:
        return None
    g = to_m(unary_union(polys)).buffer(0)
    parts = [p for p in getattr(g, 'geoms', [g]) if p.area >= MIN_PART]
    if not parts:
        return None
    g = unary_union(parts).simplify(SIMPLIFY, preserve_topology=True)
    return g


def rounded(geom_ll):
    def r(c):
        return [[round(x, 5), round(y, 5)] for x, y in c]
    gj = mapping(geom_ll)
    if gj['type'] == 'Polygon':
        gj['coordinates'] = [r(ring) for ring in gj['coordinates']]
    elif gj['type'] == 'MultiPolygon':
        gj['coordinates'] = [[r(ring) for ring in poly] for poly in gj['coordinates']]
    return gj


def label_points(region):
    """Captions sit inside each sizeable piece, near its north-west corner, where
    a map frame's caption would go. One per piece over a third of a square km."""
    out = []
    for p in sorted(getattr(region, 'geoms', [region]), key=lambda q: -q.area):
        if p.area < 0.33e6:
            continue
        inset = p.buffer(-180)
        if inset.is_empty:
            inset = p.buffer(-60)
        if inset.is_empty:
            pt = p.representative_point()
        else:
            x0, _, _, y1 = p.bounds
            best = None
            for q in getattr(inset, 'geoms', [inset]):
                for x, y in q.exterior.coords:
                    d = (x - x0) + (y1 - y)
                    if best is None or d < best[0]:
                        best = (d, x, y)
            pt = Point(best[1], best[2])
        ll = to_ll(pt)
        out.append([round(ll.y, 5), round(ll.x, 5)])
        if len(out) >= 4:
            break
    return out


def build(only=None):
    """The timeline over the manifest's layers, or over `only` (a set of layer
    ids). The static site carries fewer layers than the local server -- the
    film, and remote tiles for Today and the Wayback years -- and a year's
    regions are only right if they were cut against the layers actually there,
    so build_static.py calls this with exactly the ids it serves."""
    man = json.load(open(P('data', 'manifest.json')))
    layers = [l for l in man['layers'] if only is None or l['id'] in only]
    order = {l['id']: i for i, l in enumerate(layers)}
    geo = {}
    info = {}
    for l in layers:
        g = footprint(l)
        if g is None:
            print(f"  {l['id']}: no coverage found")
            continue
        geo[l['id']] = g
        prio = l.get('prio', 0)
        gj = P('data', f"layer_{l['id'][1:]}_geo.json") if l['id'].startswith('l') else None
        if gj and os.path.exists(gj):
            prio = json.load(open(gj)).get('prio', 0)
        info[l['id']] = dict(year=year_of(l), label=str(l['label']), group=l['group'], prio=prio,
                             credit=credit(l),
                             km2=round(g.area / 1e6, 1), footprint=rounded(to_ll(g)))
        if l.get('part_of'):
            info[l['id']]['part_of'] = l['part_of']
        print(f"  {l['id']:>8} {info[l['id']]['label']:>6}  {g.area/1e6:7.1f} km2  "
              f"{sum(len(p.exterior.coords) for p in getattr(g, 'geoms', [g]))} vertices", flush=True)

    # a part (more ground of the same picture, fetched separately) sits in the
    # stack exactly where its parent does
    for i, v in info.items():
        par = v.get('part_of')
        if par in info:
            v['prio'] = info[par]['prio']
            order[i] = order[par] + 0.5
    ids = sorted(info, key=lambda i: (info[i]['year'], info[i]['prio'], order[i]))
    years = sorted({info[i]['year'] for i in ids})
    steps = []
    for y in years:
        # newest first: each layer keeps what no newer layer (at or before y) covers
        stack = [i for i in ids if info[i]['year'] <= y][::-1]
        covered = None
        shown, vis_of = [], {}
        for i in stack:
            vis = geo[i] if covered is None else geo[i].difference(covered)
            vis = vis.buffer(-SLIVER).buffer(SLIVER).intersection(geo[i])
            covered = geo[i] if covered is None else unary_union([covered, geo[i]])
            if vis.is_empty or vis.area < MIN_PART:
                continue
            vis = vis.simplify(SIMPLIFY, preserve_topology=True)
            vis_of[i] = vis
            b = to_ll(vis).bounds
            shown.append(dict(id=i, bounds=[round(b[1], 5), round(b[0], 5), round(b[3], 5), round(b[2], 5)],
                              km2=round(vis.area / 1e6, 1)))
        # outlines and captions: one per picture, a parent and its parts together
        outlines, seen = [], {}
        for e in shown:
            key = info[e['id']].get('part_of') or e['id']
            seen.setdefault(key, []).append(vis_of[e['id']])
        for key, parts in seen.items():
            u = parts[0] if len(parts) == 1 else (            # close the hairline where parts meet
                unary_union(parts).buffer(1, join_style=2).buffer(-1, join_style=2).simplify(SIMPLIFY))
            # a part whose parent is not in this year's stack still captions as itself
            oid = key if key in info else next(e['id'] for e in shown
                                               if (info[e['id']].get('part_of') or e['id']) == key)
            outlines.append(dict(id=oid, region=rounded(to_ll(u)), labels=label_points(u)))
        label = 'Today' if y == TODAY else str(y)
        steps.append(dict(year=y, label=label, stack=shown, outlines=outlines))
        print(f"  {label:>6}: " + ', '.join(f"{info[s['id']]['label']}/{s['id']} {s['km2']}" for s in shown), flush=True)

    return dict(layers={i: info[i] for i in ids}, years=steps, today=TODAY)


def main():
    t0 = time.time()
    out = build()
    json.dump(out, open(P('data', 'timeline.json'), 'w'), separators=(',', ':'))
    print(f"  wrote data/timeline.json ({os.path.getsize(P('data', 'timeline.json'))/1e3:.0f} KB) "
          f"in {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
