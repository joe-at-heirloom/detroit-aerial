"""Control stations: correlate the historical ridge response against the modern
ridge response on a grid of overlapping windows.

Road *vectors* covered only ~1.8% of a tile, which put a floor under station
accuracy. Modern imagery is dense everywhere, so the ridge-vs-ridge match is far
better conditioned."""
import numpy as np, math, json
from PIL import Image
import absorient
Image.MAX_IMAGE_PIXELS = None


def solve_grid(arr, mod, bbox, mpp, MLAT, MLON, win_m=2600.0, overlap=0.5,
               search_m=330.0, min_valid=0.55, log=print):
    """arr: historical raster; mod: modern raster resampled to the same grid."""
    S, W, N, E = bbox
    H, Wd = arr.shape
    minE = (W - (-83.0450)) * MLON
    maxN = (N - 42.3340) * MLAT
    wp = int(win_m / mpp)
    step = max(1, int(wp * overlap))
    out = []
    for yy in range(0, max(1, H - wp), step):
        for xx in range(0, max(1, Wd - wp), step):
            y1 = min(H, yy + wp); x1 = min(Wd, xx + wp)
            a = arr[yy:y1, xx:x1]
            if a.shape[0] < wp * 0.6 or a.shape[1] < wp * 0.6:
                continue
            if (a > 0).mean() < min_valid:
                continue
            b = mod[yy:y1, xx:x1]
            if (b > 0).mean() < min_valid:
                continue
            ra, va = absorient.ridge(a.astype(np.float32), mpp)
            rb, vb = absorient.ridge(b.astype(np.float32), mpp)
            A = ra * va; B = rb * vb
            if A.std() < 1e-6 or B.std() < 1e-6:
                continue
            A = A - A[va].mean(); A[~va] = 0
            B = B - B[vb].mean(); B[~vb] = 0
            C = np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(A) * np.conj(np.fft.rfft2(B)),
                                              s=A.shape))
            lim = int(search_m / mpp)
            cy, cx = A.shape[0] // 2, A.shape[1] // 2
            sub = C[max(0, cy - lim):cy + lim, max(0, cx - lim):cx + lim]
            if sub.size == 0:
                continue
            pk = np.unravel_index(np.argmax(sub), sub.shape)
            peak = sub[pk]
            mm = sub.copy()
            ex = max(1, int(140 / mpp))
            mm[max(0, pk[0] - ex):pk[0] + ex, max(0, pk[1] - ex):pk[1] + ex] = -9e9
            second = mm.max()
            cE = minE + (xx + a.shape[1] / 2) * mpp
            cN = maxN - (yy + a.shape[0] / 2) * mpp
            out.append(dict(cE=float(cE), cN=float(cN),
                            dE=float((pk[1] - lim) * mpp),
                            dN=float(-(pk[0] - lim) * mpp),
                            rot=0.0,
                            ratio=float(peak / max(second, 1e-9)),
                            sigma=float((peak - sub.mean()) / (sub.std() + 1e-9))))
    return out


def solve_grid_hier(arr, mod, bbox, mpp, MLAT, MLON, win_m=2600.0, overlap=0.5,
                    coarse_search=200.0, fine_search=None, min_valid=0.55,
                    log=print):
    """Control stations solved coarse-to-fine, so they cannot alias.

    `solve_grid` above searches +/-330 m in one stage against a ridge map that
    contains Detroit's side streets. Those repeat every 97.5 m, so that search
    holds six alias positions per axis and a station can come back confidently one
    whole block wrong -- which then gets baked into the warp. Here the coarse stage
    sees arterial widths only (mile grid, 1609 m, cannot alias) and the fine stage
    searches less than half a block around it."""
    import gridval
    if fine_search is None:
        fine_search = gridval.FINE_SEARCH
    S, W, N, E = bbox
    H, Wd = arr.shape
    minE = (W - (-83.0450)) * MLON
    maxN = (N - 42.3340) * MLAT
    wp = int(win_m / mpp)
    step = max(1, int(wp * overlap))
    out = []
    for yy in range(0, max(1, H - wp), step):
        for xx in range(0, max(1, Wd - wp), step):
            y1 = min(H, yy + wp); x1 = min(Wd, xx + wp)
            a = arr[yy:y1, xx:x1]
            if a.shape[0] < wp * 0.6 or a.shape[1] < wp * 0.6:
                continue
            if (a > 0).mean() < min_valid:
                continue
            b = mod[yy:y1, xx:x1]
            if (b > 0).mean() < min_valid:
                continue
            r = gridval.match_hier(a, b, mpp, coarse_search=coarse_search,
                                   fine_search=fine_search)
            if r is None or r['pegged']:
                continue
            out.append(dict(cE=float(minE + (xx + a.shape[1] / 2) * mpp),
                            cN=float(maxN - (yy + a.shape[0] / 2) * mpp),
                            dE=r['dE'], dN=r['dN'], rot=0.0,
                            ratio=r['ratio'], sigma=r['sigma'],
                            coarse_ratio=r['coarse_ratio']))
    return out


def reject_outliers(st, length=3000.0, k=2.5, min_neighbours=4):
    """Drop stations that disagree with their own neighbourhood.

    RBFWarp.robust_trim does this too, but doing it here as well means the linear
    trend -- which is fitted before the RBF and therefore sees every station -- is
    not dragged by a station that survived the confidence cut but is still wrong."""
    import numpy as _np, math as _m
    if len(st) < 12:
        return st, []
    P = _np.array([[s['cE'], s['cN']] for s in st], float)
    V = _np.array([[s['dE'], s['dN']] for s in st], float)
    keep, drop = [], []
    for i in range(len(st)):
        d2 = ((P - P[i]) ** 2).sum(1)
        w = _np.exp(-d2 / (2 * length ** 2)); w[i] = 0.0
        if (w > 0.1).sum() < min_neighbours:
            keep.append(st[i]); continue
        loc = (V * w[:, None]).sum(0) / w.sum()
        resid = _m.hypot(*(V[i] - loc))
        spread = _m.sqrt(max((w * ((V - loc) ** 2).sum(1)).sum() / w.sum(), 25.0))
        (keep if resid <= k * spread else drop).append(st[i])
    return keep, drop


def solve_grid_prep(t, bbox, MLAT, MLON, win_m=2600.0, overlap=0.5,
                    coarse_search=200.0, fine_search=None, min_valid=0.55,
                    prior=None):
    """solve_grid_hier on precomputed ridge maps -- same result, ~20x faster,
    which is what makes a dense station grid and cross-validation affordable."""
    import gridval
    mpp = t['mpp']
    if fine_search is None:
        fine_search = gridval.FINE_SEARCH
    S, W, N, E = bbox
    H, Wd = t['shape']
    minE = (W - (-83.0450)) * MLON
    maxN = (N - 42.3340) * MLAT
    wp = int(win_m / mpp)
    step = max(1, int(wp * overlap))
    out = []
    for yy in range(0, max(1, H - wp), step):
        for xx in range(0, max(1, Wd - wp), step):
            y1 = min(H, yy + wp); x1 = min(Wd, xx + wp)
            if y1 - yy < wp * 0.6 or x1 - xx < wp * 0.6:
                continue
            if t['hv'][yy:y1, xx:x1].mean() < min_valid:
                continue
            if t['mv'][yy:y1, xx:x1].mean() < min_valid:
                continue
            pr = prior.at((yy + y1) / 2, (xx + x1) / 2) if prior is not None else None
            r = gridval.match_hier_prep(t, yy, y1, xx, x1,
                                        coarse_search=coarse_search,
                                        fine_search=fine_search, prior=pr)
            if r is None or r['pegged']:
                continue
            out.append(dict(cE=float(minE + (xx + (x1 - xx) / 2) * mpp),
                            cN=float(maxN - (yy + (y1 - yy) / 2) * mpp),
                            y=float(yy + (y1 - yy) / 2), x=float(xx + (x1 - xx) / 2),
                            dE=r['dE'], dN=r['dN'], rot=0.0,
                            ratio=r['ratio'], sigma=r['sigma'],
                            fine_dE=r['fine_dE'], fine_dN=r['fine_dN'],
                            coarse_ratio=r['coarse_ratio']))
    return out
