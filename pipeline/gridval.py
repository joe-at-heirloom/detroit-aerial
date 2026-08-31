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
import numpy as np, math, json
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


def _parabolic(sub, pk):
    """Quadratic fit through the peak and its two neighbours on each axis."""
    out = []
    for ax in (0, 1):
        n = sub.shape[ax]; i = pk[ax]
        if i <= 0 or i >= n - 1:
            out.append(0.0); continue
        if ax == 0:
            a, b, c = sub[i - 1, pk[1]], sub[i, pk[1]], sub[i + 1, pk[1]]
        else:
            a, b, c = sub[pk[0], i - 1], sub[pk[0], i], sub[pk[0], i + 1]
        den = a - 2.0 * b + c
        out.append(0.0 if abs(den) < 1e-12 else
                   float(np.clip(0.5 * (a - c) / den, -0.5, 0.5)))
    return out[0], out[1]


def match_cell(A, mA, B, mB, mpp, search_m=200.0, method='mncc', excl_m=None):
    # The exclusion that defines "second peak" has to sit inside the search window,
    # otherwise it blanks the entire surface and every ratio comes back as zero.
    if excl_m is None:
        excl_m = min(140.0, search_m * 0.5)
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
    # Sub-pixel peak. The modern reference is 3.5 m/px resampled to 2.5, so an
    # integer-pixel peak quantises every measurement to +/-1.25 m for no reason.
    oy, ox = _parabolic(sub, pk)
    dE = (pk[1] + ox - cx) * mpp; dN = -(pk[0] + oy - cy) * mpp
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


# --- Coarse-to-fine matching -------------------------------------------------
#
# Detroit's residential streets repeat every 97.5 m here (measured off the modern
# imagery: autocorrelation r=0.46-0.67 on the E-W axis). Any single-stage search
# wider than half that period can lock onto the wrong block and return a confident
# number that is a whole block wrong. That is the trap the project already hit for
# absolute orientation -- it comes straight back through the validator if the
# validator searches +/-200 m over a ridge map that contains side streets.
#
# So: two stages.
#   coarse  arterial-width ridge only (34-80 m). Arterials lie on the 1609 m mile
#           grid, so a +/-200 m search cannot alias. Blunt, but unambiguous.
#   fine    all-road ridge, but searched only +/-FINE_SEARCH about the coarse
#           answer. With FINE_SEARCH < half a block, aliasing is impossible.
#
# A cell whose fine peak lands on the search boundary is reported as pegged rather
# than given a number -- an honest "cannot vouch for this" beats a fabricated one.

BLOCK_M = 97.5          # measured residential period, West Detroit
FINE_SEARCH = 40.0      # < BLOCK_M/2, so the fine stage cannot alias
COARSE_BAND = (34.0, 80.0)
COARSE_FLAT = 700.0


def match_hier(a, b, mpp, coarse_search=200.0, fine_search=FINE_SEARCH,
               coarse_band=COARSE_BAND, coarse_flat=COARSE_FLAT,
               fine_band=(18.0, 46.0), fine_flat=320.0, method='mncc'):
    """a: historical crop, b: modern crop (raw uint8/float, 0 = nodata).

    Returns the displacement of a's content relative to b's, plus the diagnostics
    needed to decide whether to believe it."""
    def R(I, band, flat):
        r, v = absorient.ridge(I.astype(np.float32), mpp, band, flat)
        return (r * v).astype(np.float32), v.astype(np.float32)

    Ac, mAc = R(a, coarse_band, coarse_flat)
    Bc, mBc = R(b, coarse_band, coarse_flat)
    if Ac.std() < 1e-9 or Bc.std() < 1e-9:
        return None
    c = match_cell(Ac, mAc, Bc, mBc, mpp, coarse_search, method)
    if c is None:
        return None

    # bring a into coarse register with b, then refine inside one block
    areg = shift_image(a, -c['dE'], -c['dN'], mpp)
    Af, mAf = R(areg, fine_band, fine_flat)
    Bf, mBf = R(b, fine_band, fine_flat)
    if Af.std() < 1e-9 or Bf.std() < 1e-9:
        return None
    f = match_cell(Af, mAf, Bf, mBf, mpp, fine_search, method)
    if f is None:
        return None

    dE = c['dE'] + f['dE']; dN = c['dN'] + f['dN']
    return dict(dE=float(dE), dN=float(dN), mag=float(math.hypot(dE, dN)),
                coarse_dE=c['dE'], coarse_dN=c['dN'], coarse_ratio=c['ratio'],
                fine_dE=f['dE'], fine_dN=f['dN'], ratio=f['ratio'],
                sigma=f['sigma'], pegged=bool(f['edge']), ovl=f.get('ovl', float('nan')))


def grid_hier(arr, mod, mpp, NY=16, NX=3, min_valid=0.45, **kw):
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
            r = match_hier(a, b, mpp, **kw)
            if r is None:
                rec['skip'] = 'nopeak'; out.append(rec); continue
            rec.update(r); out.append(rec)
    return out


# --- Precomputed ridge maps --------------------------------------------------
#
# Both the validator and the station solver cut hundreds of overlapping windows
# out of the same two rasters. Running the sato filter per window re-does the same
# work hundreds of times, and the coarse band (sigmas up to ~16 px) is the slow
# one. Compute each ridge map once over the whole raster instead, in overlapping
# strips so the background subtraction never sees a false edge, and cut windows
# out of that. The fine stage's "shift then re-ridge" becomes an offset crop from
# the full map, which is exact -- nothing is lost off the edge the way it is when
# a window is shifted in isolation.

from scipy.ndimage import uniform_filter as _uf
from skimage.filters import sato as _sato


def ridge_full(I, mpp, road_w=(18.0, 46.0), flat_m=320.0, strip=6144):
    valid = I > 0
    H, W = I.shape
    out = np.zeros((H, W), np.float32)
    marg = int(2 * flat_m / mpp)
    sig = np.linspace(max(1.0, road_w[0] / (2 * mpp)), road_w[1] / (2 * mpp), 5)
    for y0 in range(0, H, strip):
        y1 = min(H, y0 + strip)
        a0 = max(0, y0 - marg); a1 = min(H, y1 + marg)
        J = I[a0:a1].astype(np.float32)
        v = valid[a0:a1]
        J = J - _uf(J, size=int(max(5, flat_m / mpp)))
        J[~v] = 0
        r = _sato(J, sigmas=sig, black_ridges=False)
        r[~v] = 0
        out[y0:y1] = r[y0 - a0:y1 - a0]
    m = out.max()
    if m > 0:
        out /= m
    return out, valid


def _ds(I, k):
    """Block-mean downsample by k. Zeros are nodata and stay near zero."""
    if k == 1:
        return I
    H, W = I.shape
    H2, W2 = H // k * k, W // k * k
    return I[:H2, :W2].reshape(H2 // k, k, W2 // k, k).mean((1, 3))


def prepare(arr, mod, mpp, coarse_band=COARSE_BAND, coarse_flat=COARSE_FLAT,
            fine_band=(18.0, 46.0), fine_flat=320.0, coarse_ds=2, log=print):
    """Ridge maps for both rasters at both bands, computed once.

    The coarse band looks for 34-80 m roads, which are 14-32 px wide at 2.5 m/px.
    Nothing about that needs full resolution, and the sato filter's cost grows with
    both pixel count and sigma, so the coarse pair is built on a 2x downsample --
    about eight times cheaper, and its answer only has to be good enough to land
    inside the fine stage's +/-40 m window."""
    t = {'mpp': mpp, 'shape': arr.shape, 'ds': coarse_ds}
    for nm, I in (('h', arr), ('m', mod)):
        r, v = ridge_full(I.astype(np.float32), mpp, fine_band, fine_flat)
        t[nm + 'f'] = r
        t[nm + 'v'] = v
        Id = _ds(I.astype(np.float32), coarse_ds)
        rc, vc = ridge_full(Id, mpp * coarse_ds, coarse_band, coarse_flat)
        t[nm + 'c'] = rc
        t[nm + 'vc'] = vc
    return t


def _cut(M, y0, y1, x0, x1):
    """Crop with zero fill outside the raster, so an offset crop is always valid."""
    H, W = M.shape
    out = np.zeros((y1 - y0, x1 - x0), M.dtype)
    a0, a1 = max(0, y0), min(H, y1)
    b0, b1 = max(0, x0), min(W, x1)
    if a1 > a0 and b1 > b0:
        out[a0 - y0:a1 - y0, b0 - x0:b1 - x0] = M[a0:a1, b0:b1]
    return out


def match_hier_prep(t, y0, y1, x0, x1, coarse_search=200.0, fine_search=FINE_SEARCH,
                    method='mncc', hshift=(0.0, 0.0), prior=None, fine_iters=3):
    """hshift plants a known displacement on the historical side without redoing
    any filtering: moving the raster and then ridge-filtering it is the same as
    ridge-filtering and then moving it, except within a background-window of the
    data edge. That makes the planted-shift tests essentially free."""
    mpp = t['mpp']; k = t.get('ds', 1); cmpp = mpp * k
    sy = int(round(hshift[1] / mpp)); sx = int(round(-hshift[0] / mpp))
    y0 += sy; y1 += sy; x0 += sx; x1 += sx
    # coarse stage runs on the downsampled pair
    cy0, cy1, cx0, cx1 = y0 // k, y1 // k, x0 // k, x1 // k
    dy, dx = sy // k, sx // k
    Ac = _cut(t['hc'], cy0, cy1, cx0, cx1)
    mAc = _cut(t['hvc'], cy0, cy1, cx0, cx1).astype(np.float32)
    Bc = _cut(t['mc'], cy0 - dy, cy1 - dy, cx0 - dx, cx1 - dx)
    mBc = _cut(t['mvc'], cy0 - dy, cy1 - dy, cx0 - dx, cx1 - dx).astype(np.float32)
    if Ac.std() < 1e-9 or Bc.std() < 1e-9:
        return None
    if prior is not None:
        # a regional prior replaces the per-cell coarse solve, which has no
        # unique peak at this window size and is free to wander
        c = dict(dE=prior[0] + hshift[0], dN=prior[1] + hshift[1],
                 ratio=float('nan'), edge=False)
    else:
        c = match_cell(Ac, mAc, Bc, mBc, cmpp, coarse_search, method)
        if c is None:
            return None
    Bf = _cut(t['mf'], y0 - sy, y1 - sy, x0 - sx, x1 - sx)
    mBf = _cut(t['mv'], y0 - sy, y1 - sy, x0 - sx, x1 - sx).astype(np.float32)
    if Bf.std() < 1e-9:
        return None
    # Exact offset crop == shifting the historical window into register. Repeat it:
    # a single fine search is capped at +/-40 m so it cannot alias, but re-centring
    # on its answer and searching again extends the reach to ~120 m while every
    # individual search stays narrower than half a block. Stops as soon as a pass
    # moves less than one pixel.
    oy = int(round(c['dN'] / mpp)); ox = int(round(c['dE'] / mpp))
    f = None
    for _ in range(max(1, fine_iters)):
        Af = _cut(t['hf'], y0 - oy, y1 - oy, x0 + ox, x1 + ox)
        mAf = _cut(t['hv'], y0 - oy, y1 - oy, x0 + ox, x1 + ox).astype(np.float32)
        if Af.std() < 1e-9 or mAf.mean() < 0.30:
            return None
        f = match_cell(Af, mAf, Bf, mBf, mpp, fine_search, method)
        if f is None:
            return None
        sxp = int(round(f['dE'] / mpp)); syp = int(round(f['dN'] / mpp))
        if sxp == 0 and syp == 0:
            break
        ox += sxp; oy += syp
        if not f['edge']:
            f = dict(f, dE=f['dE'] - sxp * mpp, dN=f['dN'] - syp * mpp, edge=False)
            break
    dE = ox * mpp + f['dE']; dN = oy * mpp + f['dN']
    return dict(dE=float(dE), dN=float(dN), mag=float(math.hypot(dE, dN)),
                coarse_dE=c['dE'], coarse_dN=c['dN'], coarse_ratio=c['ratio'],
                # total movement away from the prior, not just the last step --
                # this is what tells the field iteration whether it has converged
                fine_dE=float(dE - c['dE']), fine_dN=float(dN - c['dN']),
                ratio=f['ratio'], sigma=f['sigma'], pegged=bool(f['edge']),
                ovl=f.get('ovl', float('nan')))


def grid_hier_prep(t, NY=16, NX=3, min_valid=0.45, prior=None, **kw):
    H, W = t['shape']
    ch = H // NY; cw = W // NX
    out = []
    for j in range(NY):
        for i in range(NX):
            y0, y1 = j * ch, (j + 1) * ch
            x0, x1 = i * cw, (i + 1) * cw
            vfrac = float(t['hv'][y0:y1, x0:x1].mean())
            rec = dict(row=j, col=i, valid=vfrac)
            if vfrac < min_valid:
                rec['skip'] = 'coverage'; out.append(rec); continue
            pr = prior.at((y0 + y1) / 2, (x0 + x1) / 2) if prior is not None else None
            r = match_hier_prep(t, y0, y1, x0, x1, prior=pr, **kw)
            if r is None:
                rec['skip'] = 'nopeak'; out.append(rec); continue
            rec.update(r); out.append(rec)
    return out


# --- Regional coarse field ---------------------------------------------------
#
# Running the coarse stage per cell does not work: one 3.6 x 2.0 km cell contains
# only about two arterials per axis, so the arterial correlation has no unique peak
# (measured coarse_ratio 1.03-1.25, i.e. nothing) and is free to jump 150 m. The
# fine stage then refines confidently around a wrong centre, which is worse than
# not having a coarse stage at all.
#
# The displacement being corrected is smooth at the kilometre scale, so the coarse
# estimate does not have to be local. Solve it on windows large enough to actually
# contain the mile grid, keep only the windows that locked, and interpolate. Each
# cell then gets a prior it can trust and only has to search +/-40 m around it.

def coarse_field(t, win_m=6000.0, overlap=0.5, search_m=250.0,
                 min_ratio=1.15, min_valid=0.40, log=print):
    k = t.get('ds', 1); cmpp = t['mpp'] * k
    H, W = t['hc'].shape
    wp = max(8, int(win_m / cmpp))
    step = max(1, int(wp * overlap))
    pts = []
    ys = list(range(0, max(1, H - wp + step), step))
    xs = list(range(0, max(1, W - wp + step), step))
    for yy in ys:
        for xx in xs:
            y1 = min(H, yy + wp); x1 = min(W, xx + wp)
            if y1 - yy < wp * 0.5 or x1 - xx < wp * 0.5:
                continue
            mA = t['hvc'][yy:y1, xx:x1].astype(np.float32)
            if mA.mean() < min_valid:
                continue
            A = t['hc'][yy:y1, xx:x1]
            B = t['mc'][yy:y1, xx:x1]
            mB = t['mvc'][yy:y1, xx:x1].astype(np.float32)
            if A.std() < 1e-9 or B.std() < 1e-9:
                continue
            r = match_cell(A, mA, B, mB, cmpp, search_m, 'mncc')
            if r is None or r['edge']:
                continue
            pts.append(dict(y=(yy + y1) / 2 * k, x=(xx + x1) / 2 * k,
                            dE=r['dE'], dN=r['dN'], ratio=r['ratio']))
    good = [p for p in pts if p['ratio'] >= min_ratio]
    # trim samples that disagree with the rest -- one window that locked onto the
    # wrong thing would otherwise pull the prior for a whole region
    if len(good) >= 5:
        V = np.array([[p['dE'], p['dN']] for p in good])
        c = np.median(V, axis=0)
        for _ in range(5):
            d = np.hypot(V[:, 0] - c[0], V[:, 1] - c[1])
            s = max(np.median(d) * 1.4826, 20.0)
            w = 1.0 / (1.0 + (d / (2 * s)) ** 2)
            c = (V * w[:, None]).sum(0) / w.sum()
        keep = np.hypot(V[:, 0] - c[0], V[:, 1] - c[1]) <= max(4 * s, 80.0)
        good = [p for p, kp in zip(good, keep) if kp]
    log(f"  coarse field: {len(pts)} windows, {len(good)} locked "
        f"(ratio >= {min_ratio})")
    return good


class Prior:
    """Gaussian-weighted interpolation of the locked coarse windows."""

    def __init__(self, pts, mpp, length_m=5000.0):
        self.mpp = mpp
        self.L = length_m / mpp
        if pts:
            self.P = np.array([[p['y'], p['x']] for p in pts], float)
            self.V = np.array([[p['dE'], p['dN']] for p in pts], float)
            self.W = np.array([min(p['ratio'], 3.0) for p in pts], float)
            self.fallback = np.median(self.V, axis=0)
        else:
            self.P = np.zeros((0, 2)); self.V = np.zeros((0, 2))
            self.W = np.zeros(0); self.fallback = np.zeros(2)

    def at(self, y, x):
        if len(self.P) == 0:
            return (0.0, 0.0)
        d2 = (self.P[:, 0] - y) ** 2 + (self.P[:, 1] - x) ** 2
        w = self.W * np.exp(-d2 / (2 * self.L ** 2))
        s = w.sum()
        if s < 1e-9:
            return tuple(self.fallback)
        return tuple((self.V * w[:, None]).sum(0) / s)


def prior_for(t, mpp, min_locked=20, log=print, **kw):
    """A regional prior, but only when enough of it locked.

    The coarse field is solved on 6 km windows of mile-grid arterials. Where those
    are thin -- 1956 runs west into country where the mile grid weakens, and its
    coverage is patchier -- few windows lock and the interpolated field is
    confidently wrong. That is worse than having no prior at all, because every
    cell then refines within +/-40 m of a wrong centre: measured on 1956's
    frame-corrected composite, a 12-window prior reported 101 m where a zero prior
    reported 11 m, for the same raster. 1961, whose field locks 28 windows, gives
    the same answer either way.

    A prior is only needed when the mosaic is further out than the fine search can
    walk. After the per-frame stage it never is."""
    pts = coarse_field(t, log=log, **kw)
    if len(pts) < min_locked:
        log(f"  coarse field locked only {len(pts)} windows (< {min_locked}); "
            f"using a zero prior instead -- a sparse one is worse than none")
        return Prior([], mpp)
    return Prior(pts, mpp)


def _sig(a):
    """Cheap content signature. A cached ridge map must never be reused after the
    raster underneath it changes -- and a re-warped mosaic has the same shape as
    the one it replaced, so shape alone is not enough to notice."""
    a = np.asarray(a)
    return [int(a.shape[0]), int(a.shape[1]),
            float(a[::97, ::89].astype(np.float64).sum())]


def _modsig(mod):
    return _sig(mod)


def prepare_cached(arr, mod, mpp, key, cachedir='/tmp/das_ridge', log=print, **kw):
    """prepare(), memoised on disk. The ridge maps take ~2 minutes to build and are
    identical every run, which otherwise dominates the cost of any experiment."""
    import os
    d = os.path.join(cachedir, key)
    names = ('hf', 'hv', 'hc', 'hvc', 'mf', 'mv', 'mc', 'mvc')
    meta = os.path.join(d, 'meta.json')
    if os.path.exists(meta):
        try:
            m = json.load(open(meta))
            if (m['mpp'] == mpp and tuple(m['shape']) == tuple(arr.shape)
                    and m.get('modsig') == _sig(mod)
                    and m.get('histsig') == _sig(arr)):
                t = {'mpp': mpp, 'shape': tuple(m['shape']), 'ds': m['ds']}
                for n in names:
                    t[n] = np.load(os.path.join(d, n + '.npy'), mmap_mode='r')
                log(f"  ridge maps from cache {d}")
                return t
        except Exception:
            pass
    t = prepare(arr, mod, mpp, log=log, **kw)
    os.makedirs(d, exist_ok=True)
    for n in names:
        np.save(os.path.join(d, n + '.npy'), np.ascontiguousarray(t[n]))
    json.dump(dict(mpp=mpp, shape=list(arr.shape), ds=t['ds'],
                   modsig=_sig(mod), histsig=_sig(arr)), open(meta, 'w'))
    return t
