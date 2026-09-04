#!/usr/bin/env python
"""Where on the map is this negative?

Expansion needs positions for frames no catalogue has ever placed, and the roll's
own numbering only carries you as far as the runs you have already solved. This
puts a single raw scan on the map by brute force: downsample it to the same coarse
grid as a placed mosaic, high-pass both so tone and haze drop out, and take the
normalised cross-correlation over the whole mosaic by FFT. The peak-to-runner-up
ratio says whether to believe it -- Detroit's grid repeats every 97.5 m, so a
weak peak means "somewhere along this street", not a position.

Usage:
  ./.venv/bin/python scripts/locate.py scans/fullres/861.jpg --against 1949:placedE
  ./.venv/bin/python scripts/locate.py --roll ha-16 --year 1949 --every 10 --against 1949:placedE
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import rasterio
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
from validate import P
from rewarp import SP
import dtmap

MPP = 6.0            # coarse enough to be fast, fine enough for arterials
GSD = 0.638          # metres per pixel on the raw scan (2300 m / 3602 px)


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def hp(a, k=9):
    from scipy.ndimage import uniform_filter
    a = a.astype(np.float32)
    m = a > 0
    b = a - uniform_filter(a, k)
    b[~m] = 0.0
    s = b[m].std() if m.any() else 1.0
    return b / (s + 1e-6), m.astype(np.float32)


def load_ref(ref):
    tag, suf = (ref.split(':', 1) + ['placedE'])[:2]
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{suf}.tif'))
    g = json.load(open(P('data', f'{tag}_{suf}_geo.json')))
    S, W, N, E = g['bbox']
    wm = (E - W) * dtmap.MLON; hm = (N - S) * dtmap.MLAT
    w = int(wm / MPP); h = int(hm / MPP)
    a = ds.read(1, out_shape=(h, w), resampling=rasterio.enums.Resampling.average)
    return a, (S, W, N, E), w, h


def locate(path, R, Rm, bbox, w, h):
    im = Image.open(path).convert('L')
    cw, ch = im.size
    im = im.crop((int(cw * .28), int(ch * .28), int(cw * .72), int(ch * .72)))   # the middle, least tilt
    tw = max(8, int(im.size[0] * GSD / MPP)); th = max(8, int(im.size[1] * GSD / MPP))
    t = np.asarray(im.resize((tw, th), Image.LANCZOS), dtype=np.float32)
    T, _ = hp(t)
    T -= T.mean()
    from numpy.fft import rfft2, irfft2
    sh = (h, w)
    F = rfft2(R, sh) * np.conj(rfft2(T, sh))
    c = irfft2(F, sh)
    # normalise by local reference energy so bright ground does not win
    E2 = irfft2(rfft2(R * R, sh) * np.conj(rfft2(np.ones_like(T), sh)), sh)
    cov = irfft2(rfft2(Rm, sh) * np.conj(rfft2(np.ones_like(T), sh)), sh)
    ok = cov > 0.85 * T.size
    score = np.where(ok, c / np.sqrt(np.maximum(E2, 1e-3)), -1e9)
    i = int(np.argmax(score)); y, x = divmod(i, w)
    best = score[y, x]
    m = score.copy()
    y0, y1 = max(0, y - th), min(h, y + th); x0, x1 = max(0, x - tw), min(w, x + tw)
    m[y0:y1, x0:x1] = -1e9
    second = float(m.max())
    S_, W_, N_, E_ = bbox
    lat = N_ - (y + th / 2) * MPP / dtmap.MLAT
    lon = W_ + (x + tw / 2) * MPP / dtmap.MLON
    return lat, lon, float(best), float(best / second if second > 0 else 99.0)


def main():
    ref = arg('--against', '1949:placedE')
    R0, bbox, w, h = load_ref(ref)
    R, Rm = hp(R0)
    print(f"reference {ref}: {w}x{h} at {MPP} m/px", flush=True)
    if '--roll' in sys.argv:
        year = arg('--year', '1949'); every = int(arg('--every', 10))
        cat = [i for i in json.load(open(P('data', 'dte_catalogue.json')))['items']
               if i['county'] == arg('--county', 'Wayne County') and i['year'] == year
               and i['section'] and i['section'].rsplit('-', 1)[0] == arg('--roll')]
        cat.sort(key=lambda i: int(i['section'].rsplit('-', 1)[1]))
        sel = [i for k, i in enumerate(cat) if k % every == 0]
        paths = [(i['section'], f"{SP}/fullres/{i['pointer']}.jpg") for i in sel]
    else:
        paths = [(os.path.basename(p), p) for p in sys.argv[1:] if p.endswith('.jpg')]
    for name, p in paths:
        if not os.path.exists(p):
            print(f"  {name}: no scan"); continue
        lat, lon, sc, ratio = locate(p, R, Rm, bbox, w, h)
        verdict = 'LOCKED  ' if ratio > 1.25 else 'weak    '
        print(f"  {name:14s} {verdict} {lat:.4f}, {lon:.4f}   peak/next {ratio:.2f}", flush=True)


if __name__ == '__main__':
    main()
