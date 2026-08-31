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
