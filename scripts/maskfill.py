#!/usr/bin/env python
"""Turn a layer's WHITE fill into nodata.

Some ImageServers fill outside their coverage with white instead of honouring
noData=0: the City of Detroit's 2010 ortho is 17 % white over our extent (all of
Dearborn), 2005 about 4 %; Wayne County's tile caches are white outside the
county, Oakland's 2023/2025 white outside theirs. serve.py draws only black as
transparent, so that white hid every older year in the year view and showed as a
white sheet in the wipe.

The fill is not clean white: it is JPEG-noisy near-white with grey seams where
the server's own tiles meet, and an anti-aliased fringe where it meets the
picture. A first version killed only pixels white on every band (>= 250) and
left the noise behind as a speckled sheet over the year beneath (2010 over
Dearborn) and thin white slivers along the edges. So:

  1. A coarse grid (16 px cells) marks cells that are almost entirely near-white
     (every band >= 235, bands within 10 of each other) or already nodata.
     Connected groups of them over MIN_KM2 that reach nodata or the raster's
     edge are fill, grown by one cell (a big white roof reaches neither). Counting
     nodata as fill-compatible makes white attached to nodata count, and makes a
     rerun on an already-masked file find the same regions.
  2. At full resolution, inside that mask: near-white pixels become 0, and then
     anything only a few pixels thick (specks, seams, slivers) -- what a 5x5
     opening of the coverage mask removes -- becomes 0 too. A picture's real
     edge is a large shape and survives the opening; a bright roof is a few
     cells at most, never a fill region, and is left alone.

    ./.venv/bin/python scripts/maskfill.py 2010 [2005 ...]
rewrites mosaics/layer_<name>.tif (via a temporary file) and its overviews.
"""
import os, sys, time
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window
from scipy import ndimage
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def P(*a): return os.path.join(ROOT, *a)

C = 16            # coarse cell, px
NEAR = 235        # every band at or above this ...
SPREAD = 10       # ... and within this of each other: near-white, neutral
MIN_KM2 = 0.05
OPEN = 5          # px: thinner than this inside the fill is noise
M = 8             # px of overlap between strips, for the opening


def whiteish(a):
    mn = a.min(axis=0)
    return (mn >= NEAR) & ((a.max(axis=0) - mn) <= SPREAD)


def one(name):
    t0 = time.time()
    src_p = P('mosaics', f'layer_{name}.tif')
    with rasterio.open(src_p) as ds:
        H, W = ds.height, ds.width
        mpp = abs(ds.transform.e) * 111000          # rough metres per pixel, for the area rule
        h, w = (H + C - 1) // C, (W + C - 1) // C
        fw = np.zeros((h, w), np.float32)           # near-white fraction per cell
        fz = np.zeros((h, w), np.float32)           # nodata fraction per cell
        for j in range(0, h, 64):
            jh = min(64, h - j)
            a = ds.read([1, 2, 3], window=Window(0, j * C, W, min(jh * C, H - j * C)))
            pad = ((0, jh * C - a.shape[1]), (0, w * C - a.shape[2]))
            wt = np.pad(whiteish(a), pad); zt = np.pad(a.max(axis=0) == 0, pad, constant_values=True)
            fw[j:j + jh] = wt.reshape(jh, C, w, C).mean(axis=(1, 3))
            fz[j:j + jh] = zt.reshape(jh, C, w, C).mean(axis=(1, 3))
        cells = (fw + fz) > 0.9
        lab, n = ndimage.label(cells)
        fill = np.zeros_like(cells)
        if n:
            sizes = ndimage.sum(np.ones_like(lab), lab, range(1, n + 1)) * (C * mpp) ** 2 / 1e6
            white = ndimage.sum(fw, lab, range(1, n + 1)) * (C * mpp) ** 2 / 1e6
            # fill is outside the picture, so it reaches nodata or the raster's
            # edge; a big saturated roof is surrounded by picture and is kept
            edge = np.zeros_like(cells); edge[[0, -1], :] = True; edge[:, [0, -1]] = True
            touch = ndimage.maximum(((fz > 0.5) | edge).astype(np.uint8), lab, range(1, n + 1)) > 0
            keep = np.zeros(n + 1, bool); keep[1:] = (sizes >= MIN_KM2) & touch
            fill = keep[lab]
            print(f"  {name}: {fw.mean()*100:.1f}% near-white, {(white[keep[1:]]).sum():.1f} km2 of it in "
                  f"{int(keep.sum())} fill regions", flush=True)
        fill = ndimage.binary_dilation(fill, iterations=1)
        if not fill.any():
            print(f"  {name}: no fill", flush=True)
            return
        prof = ds.profile.copy()
        # .profile does not carry these; fetchlayer.py writes them
        prof.update(predictor=2, photometric='RGB', BIGTIFF='YES', num_threads='ALL_CPUS')
        tmp = src_p + '.tmp.tif'
        killed = 0
        st = np.ones((OPEN, OPEN), bool)
        with rasterio.open(tmp, 'w', **prof) as dst:
            for j in range(0, h, 32):
                jh = min(32, h - j)
                y0 = j * C; hh = min(jh * C, H - y0)
                r0 = max(0, y0 - M); r1 = min(H, y0 + hh + M)          # read with overlap
                a = ds.read([1, 2, 3], window=Window(0, r0, W, r1 - r0))
                c0, c1 = r0 // C, (r1 - 1) // C + 1
                m = np.repeat(np.repeat(fill[c0:c1], C, axis=0), C, axis=1)
                m = m[r0 - c0 * C:r0 - c0 * C + (r1 - r0), :W]
                k1 = m & whiteish(a)
                nz = (a.max(axis=0) > 0) & ~k1
                k2 = m & nz & ~ndimage.binary_opening(nz, structure=st)
                kill = (k1 | k2)[y0 - r0:y0 - r0 + hh]
                out = a[:, y0 - r0:y0 - r0 + hh]
                out[:, kill] = 0
                killed += int(kill.sum())
                dst.write(out, window=Window(0, y0, W, hh))
            print(f"  {name}: {killed/1e6:.1f} Mpx cleared; overviews", flush=True)
            dst.build_overviews([2, 4, 8, 16, 32, 64], Resampling.average)
    os.replace(tmp, src_p)
    print(f"  {name}: done in {time.time()-t0:.0f}s", flush=True)


if __name__ == '__main__':
    for n in sys.argv[1:]:
        one(n)
