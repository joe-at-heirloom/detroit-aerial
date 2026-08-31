"""Grid validation that reports WHY a cell scored what it scored.

The old metric was a raw, unnormalised FFT cross-correlation of two ridge maps.
Where the historical mosaic only partly covers a cell, that correlation is biased
toward shifts which align the *coverage boundary* rather than the roads, so a cell
can report 200 m of error without the imagery being 200 m out -- and a cell can
report a small number for the same bad reason. Every number here therefore comes
with the diagnostics needed to decide whether it means anything:

  valid   fraction of the cell the historical mosaic actually covers
  ratio   peak / best peak outside a 140 m exclusion  (1.0 = no unique peak)
  edge    the peak landed on the search-window boundary, i.e. unbounded
  ovl     fraction of the cell overlapping at the winning shift

`method='mncc'` is Padfield masked normalised cross-correlation: it divides out
the overlap area and the local means/variances, so partial coverage stops biasing
the peak. `method='raw'` reproduces the old metric for comparison.
"""
import numpy as np, math
import absorient

try:
    from scipy.fft import next_fast_len
except Exception:
    def next_fast_len(n):
        while True:
            m = n
            for p in (2, 3, 5):
                while m % p == 0:
                    m //= p
            if m == 1:
                return n
            n += 1


def shift_image(I, dE_m, dN_m, mpp):
    """Move image content by (dE, dN) metres. +dE east, +dN north."""
    sx = int(round(dE_m / mpp)); sy = int(round(-dN_m / mpp))
    out = np.zeros_like(I)
    h, w = I.shape
    ys0, ys1 = max(0, -sy), min(h, h - sy)
    xs0, xs1 = max(0, -sx), min(w, w - sx)
    out[ys0 + sy:ys1 + sy, xs0 + sx:xs1 + sx] = I[ys0:ys1, xs0:xs1]
    return out


def _pad_fft(x, H, W):
    p = np.zeros((H, W), np.float32)
    p[:x.shape[0], :x.shape[1]] = x
    return np.fft.rfft2(p)


def _mncc(A, mA, B, mB, ly, lx, min_ovl=0.30):
    """Masked NCC surface, cropped to shifts within +/-(ly,lx), centre = zero shift.

    Convention matches the old code: corr = F(A) * conj(F(B)), so a positive x
    index past centre means A must move east to meet B."""
    h, w = A.shape
    H = next_fast_len(h + 2 * ly); W = next_fast_len(w + 2 * lx)
    FA, FmA, FAA = _pad_fft(A, H, W), _pad_fft(mA, H, W), _pad_fft(A * A, H, W)
    FB, FmB, FBB = _pad_fft(B, H, W), _pad_fft(mB, H, W), _pad_fft(B * B, H, W)
    def C(X, Y):
        return np.fft.irfft2(X * np.conj(Y), s=(H, W))
    N   = C(FmA, FmB)
    SA  = C(FA,  FmB)
    SB  = C(FmA, FB)
    SAA = C(FAA, FmB)
    SBB = C(FmA, FBB)
    SAB = C(FA,  FB)
    N = np.maximum(N, 1e-6)
    num = SAB - SA * SB / N
    va = np.maximum(SAA - SA * SA / N, 0.0)
    vb = np.maximum(SBB - SB * SB / N, 0.0)
    den = np.sqrt(va * vb)
    r = np.where(den > 1e-9, num / np.maximum(den, 1e-9), 0.0)
    r[N < min_ovl * h * w] = 0.0
    r = np.fft.fftshift(r)
    Nc = np.fft.fftshift(N)
    cy, cx = H // 2, W // 2
    return (r[cy - ly:cy + ly + 1, cx - lx:cx + lx + 1],
            Nc[cy - ly:cy + ly + 1, cx - lx:cx + lx + 1] / (h * w))


def _raw(A, mA, B, mB, ly, lx):
    """The old metric: mean-subtract inside the mask, zero outside, plain correlation."""
    a = A.copy(); b = B.copy()
    va = mA > 0.5; vb = mB > 0.5
    if va.sum() == 0 or vb.sum() == 0:
        return None, None
    a = a - a[va].mean(); a[~va] = 0
    b = b - b[vb].mean(); b[~vb] = 0
    C = np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(a) * np.conj(np.fft.rfft2(b)), s=a.shape))
    cy, cx = a.shape[0] // 2, a.shape[1] // 2
    return C[max(0, cy - ly):cy + ly + 1, max(0, cx - lx):cx + lx + 1], None


def match_cell(A, mA, B, mB, mpp, search_m=200.0, method='mncc', excl_m=140.0):
    ly = lx = int(search_m / mpp)
    if method == 'mncc':
        sub, ovl = _mncc(A, mA, B, mB, ly, lx)
    else:
        sub, ovl = _raw(A, mA, B, mB, ly, lx)
    if sub is None or sub.size == 0 or not np.isfinite(sub).any():
        return None
    pk = np.unravel_index(np.argmax(sub), sub.shape)
    peak = float(sub[pk])
    mm = sub.copy()
    ex = max(1, int(excl_m / mpp))
    mm[max(0, pk[0] - ex):pk[0] + ex, max(0, pk[1] - ex):pk[1] + ex] = -9e9
    second = float(mm.max())
    cy = sub.shape[0] // 2; cx = sub.shape[1] // 2
    dE = (pk[1] - cx) * mpp; dN = -(pk[0] - cy) * mpp
    edge = (pk[0] <= 1 or pk[1] <= 1 or pk[0] >= sub.shape[0] - 2 or pk[1] >= sub.shape[1] - 2)
    return dict(dE=float(dE), dN=float(dN), mag=float(math.hypot(dE, dN)),
                peak=peak, ratio=float(peak / max(abs(second), 1e-9)),
                sigma=float((peak - sub.mean()) / (sub.std() + 1e-9)),
                edge=bool(edge),
                ovl=float(ovl[pk]) if ovl is not None else float('nan'))


def grid(arr, mod, mpp, NY=16, NX=3, search_m=200.0, min_valid=0.45,
         method='mncc', prefilt=True):
    """arr: historical raster (0 = no data). mod: modern, same grid."""
    Hp, Wp = arr.shape
    ch = Hp // NY; cw = Wp // NX
    out = []
    for j in range(NY):
        for i in range(NX):
            a = arr[j * ch:(j + 1) * ch, i * cw:(i + 1) * cw]
            b = mod[j * ch:(j + 1) * ch, i * cw:(i + 1) * cw]
            vfrac = float((a > 0).mean())
            rec = dict(row=j, col=i, valid=vfrac)
            if vfrac < min_valid:
                rec['skip'] = 'coverage'; out.append(rec); continue
            ra, va = absorient.ridge(a.astype(np.float32), mpp)
            rb, vb = absorient.ridge(b.astype(np.float32), mpp)
            if ra.std() < 1e-6 or rb.std() < 1e-6:
                rec['skip'] = 'flat'; out.append(rec); continue
            A = (ra * va).astype(np.float32); B = (rb * vb).astype(np.float32)
            r = match_cell(A, va.astype(np.float32), B, vb.astype(np.float32),
                           mpp, search_m, method)
            if r is None:
                rec['skip'] = 'nopeak'; out.append(rec); continue
            rec.update(r); out.append(rec)
    return out


def summarise(recs, label='', usable=lambda r: 'skip' not in r):
    u = [r for r in recs if usable(r)]
    if not u:
        return f"{label}: no usable cells"
    m = np.array([r['mag'] for r in u])
    return (f"{label}: n={len(u)}/{len(recs)}  median {np.median(m):6.1f}  "
            f"p90 {np.percentile(m, 90):6.1f}  max {m.max():6.1f}  "
            f">25m {(m > 25).sum()}  >10m {(m > 10).sum()}")


def arterial_grid(arr, mpp, minE, maxN, NY=16, NX=3, search_m=200.0,
                  min_valid=0.45, min_div=0.25):
    """Independent per-cell check against the mile-grid arterials.

    Ridge-vs-modern-imagery can alias onto Detroit's ~100 m residential lattice
    when the historical and modern content differ enough that the true peak is not
    dominant. The arterials repeat every 1609 m, so this cannot alias inside a
    +/-200 m search. It is blunt -- roughly 5 m of noise floor, and it needs
    arterials running in two directions -- but it settles the only question that
    matters here: is a cell a few metres out, or a whole block out?"""
    S = absorient.arterial_segments(use_osm=False)
    Hp, Wp = arr.shape
    ch = Hp // NY; cw = Wp // NX
    out = []
    for j in range(NY):
        for i in range(NX):
            a = arr[j * ch:(j + 1) * ch, i * cw:(i + 1) * cw]
            vfrac = float((a > 0).mean())
            rec = dict(row=j, col=i, valid=vfrac)
            if vfrac < min_valid:
                rec['skip'] = 'coverage'; out.append(rec); continue
            cminE = minE + i * cw * mpp
            cmaxN = maxN - j * ch * mpp
            cE = cminE + a.shape[1] * mpp / 2
            cN = cmaxN - a.shape[0] * mpp / 2
            rad = max(a.shape) * mpp / 2
            div, nseg = absorient.bearing_diversity(S, cE, cN, rad)
            rec.update(div=float(div), nseg=int(nseg))
            if div < min_div:
                rec['skip'] = 'parallel-arterials-only'; out.append(rec); continue
            r, va = absorient.ridge(a.astype(np.float32), mpp)
            res = absorient.solve(r, va, S, cminE, cmaxN, mpp, mpp, search_m=search_m)
            if res is None:
                rec['skip'] = 'no-arterials'; out.append(rec); continue
            res['mag'] = float(math.hypot(res['dE'], res['dN']))
            rec.update(res); out.append(rec)
    return out
