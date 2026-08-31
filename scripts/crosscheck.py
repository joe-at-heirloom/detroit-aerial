#!/usr/bin/env python
"""Measure the same alignment with four independent implementations.

"Trust nothing" applied to the measurement itself. Everything reported so far
comes from one correlation routine written for this project. If that routine has a
systematic error, every number inherits it and no amount of internal
self-consistency would reveal it.

So the same check windows are measured four ways, by three third-party libraries
that share no code with this project and, in two cases, no algorithm family:

  ours      masked normalised cross-correlation of a ridge response (pipeline/gridval)
  skimage   scikit-image phase_cross_correlation, masked variant (Padfield)
  opencv    ECC maximisation -- gradient-based, not correlation-peak-based
  sitk      SimpleITK, Mattes mutual information -- an information-theoretic
            criterion that does not assume the two images have similar brightness
            at all, which is the right prior for 1949 film against modern satellite

Each backend is pinned by a planted-shift test before it is trusted, so its sign
and scale conventions are established empirically rather than assumed.

A caution that matters when reading the output. scikit-image's phase correlation
searches the *whole* image with no bound, so on a street grid that repeats every
97.5 m it is free to lock a block or several off -- and it does, returning 150-280 m
on windows where the bounded methods agree on ~20 m. That is the project's oldest
trap wearing a third-party badge. Results beyond `--bound` are therefore reported
as out-of-bound rather than counted as error. OpenCV's ECC and SimpleITK's mutual
information are both local, so they are the meaningful independent checks here.

Usage:  ./.venv/bin/python scripts/crosscheck.py --selftest
        ./.venv/bin/python scripts/crosscheck.py --block 1961 [--n 16]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, dtmap
import validate
from validate import P
import naipcheck
import rasterio
from rasterio.windows import Window


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def _ridge(a, mpp):
    """Road-like linear structure, via scikit-image's sato filter.

    Registering the raw film against modern imagery asks the two epochs to agree
    about trees, roof colour and land cover, all of which genuinely changed. The
    question this project actually cares about is whether the *streets* line up, so
    the backends can be handed a road response instead. The filter is
    scikit-image's, not ours, so the comparison stays independent of this project's
    correlation code."""
    r, v = gridval.ridge_full(a.astype(np.float32), mpp)
    out = (r * 255.0).astype(np.float32)
    out[~v] = 0
    return out


def _norm(a):
    a = a.astype(np.float32)
    v = a > 0
    if v.sum() < 50:
        return a, v
    m, s = a[v].mean(), a[v].std() + 1e-6
    o = (a - m) / s
    o[~v] = 0
    return o, v


# --- backends: each returns (dE, dN) in metres, historical relative to reference --

def m_ours(a, b, mpp, search_m=60.0):
    r = naipcheck.measure(a, b, mpp, search_m)
    if r is None or r['edge'] or r['ratio'] < 1.10:
        return None
    return r['dE'], r['dN']


def m_skimage(a, b, mpp, search_m=60.0):
    from skimage.registration import phase_cross_correlation
    A, va = _norm(a); B, vb = _norm(b)
    try:
        sh = phase_cross_correlation(B, A, reference_mask=vb, moving_mask=va,
                                     overlap_ratio=0.5)
    except TypeError:
        sh = phase_cross_correlation(B, A, reference_mask=vb, moving_mask=va)
    sh = sh[0] if isinstance(sh, tuple) else sh
    # phase_cross_correlation returns the shift that registers `moving` onto
    # `reference`; the displacement of the moving content is its negation. Pinned
    # by the self-test, not assumed.
    dy, dx = -float(sh[0]), -float(sh[1])
    return dx * mpp, -dy * mpp


def m_opencv(a, b, mpp, search_m=60.0):
    import cv2
    A, _ = _norm(a); B, _ = _norm(b)
    crit = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 500, 1e-7)
    # ECC is a local gradient method and will not find a 20 m shift from a cold
    # start, so it gets its own coarse-to-fine ladder rather than a hint from us
    warp = np.eye(2, 3, dtype=np.float32)
    for sc in (0.25, 0.5, 1.0):
        As = cv2.resize(A, None, fx=sc, fy=sc, interpolation=cv2.INTER_AREA)
        Bs = cv2.resize(B, None, fx=sc, fy=sc, interpolation=cv2.INTER_AREA)
        w = warp.copy(); w[0, 2] *= sc; w[1, 2] *= sc
        try:
            cv2.findTransformECC(Bs, As, w, cv2.MOTION_TRANSLATION, crit, None, 5)
        except cv2.error:
            return None
        warp = w.copy(); warp[0, 2] /= sc; warp[1, 2] /= sc
    dx, dy = float(warp[0, 2]), float(warp[1, 2])
    return dx * mpp, -dy * mpp


def m_sitk(a, b, mpp, search_m=60.0):
    import SimpleITK as sitk
    A, _ = _norm(a); B, _ = _norm(b)
    fx = sitk.GetImageFromArray(B.astype(np.float32))
    mv = sitk.GetImageFromArray(A.astype(np.float32))
    R = sitk.ImageRegistrationMethod()
    R.SetMetricAsMattesMutualInformation(numberOfHistogramBins=48)
    R.SetMetricSamplingStrategy(R.RANDOM)
    R.SetMetricSamplingPercentage(0.20, seed=7)
    R.SetInterpolator(sitk.sitkLinear)
    R.SetOptimizerAsRegularStepGradientDescent(8.0, 0.005, 600)
    R.SetOptimizerScalesFromPhysicalShift()
    R.SetInitialTransform(sitk.TranslationTransform(2), inPlace=False)
    R.SetShrinkFactorsPerLevel([8, 4, 2, 1])
    R.SetSmoothingSigmasPerLevel([4, 2, 1, 0])
    try:
        t = R.Execute(fx, mv)
    except RuntimeError:
        return None
    p = t.GetParameters()
    # sign pinned by the self-test
    return float(p[0]) * mpp, -float(p[1]) * mpp


BACKENDS = [('ours', m_ours), ('skimage', m_skimage),
            ('opencv', m_opencv), ('sitk', m_sitk)]


# SimpleITK is retained for the record but excluded by default: the real-data
# control measured it at 43.69 m median error recovering a planted shift, 0% within
# 5 m. Mutual information cannot register mid-century film against modern
# orthoimagery here, and it costs more runtime than the two that work. --with-sitk
# puts it back.
def available():
    out = []
    skip = set() if '--with-sitk' in sys.argv else {'sitk'}
    for name, fn in BACKENDS:
        if name in skip:
            continue
        try:
            if name == 'skimage':
                from skimage.registration import phase_cross_correlation  # noqa
            elif name == 'opencv':
                import cv2  # noqa
            elif name == 'sitk':
                import SimpleITK  # noqa
            out.append((name, fn))
        except Exception:
            pass
    return out


def selftest(mpp=1.0):
    """Pin every backend's sign convention on a planted shift before trusting it."""
    rng = np.random.default_rng(11)
    H = W = 700
    img = np.full((H, W), 60, np.float32)
    for y in range(0, H, 90): img[y:y + 5, :] = 200
    for x in range(0, W, 90): img[:, x:x + 5] = 200
    for _ in range(500):
        y, x = rng.integers(0, H - 16), rng.integers(0, W - 16)
        img[y:y + rng.integers(6, 15), x:x + rng.integers(6, 15)] = rng.integers(90, 245)
    img = np.clip(img + rng.normal(0, 5, (H, W)).astype(np.float32), 1, 255)
    print(f"Backend self-test ({mpp} m/px). Every backend must recover the plant.")
    ok_all = True
    for pE, pN in ((0.0, 0.0), (12.0, 0.0), (0.0, 12.0), (-18.0, 9.0)):
        hist = gridval.shift_image(img, pE, pN, mpp)
        row = f"  planted ({pE:+6.1f},{pN:+6.1f}) ->"
        for name, fn in available():
            try:
                r = fn(hist, img, mpp)
            except Exception as e:
                r = None
            if r is None:
                row += f"  {name}: FAIL"; ok_all = False
            else:
                good = abs(r[0] - pE) < 3 and abs(r[1] - pN) < 3
                ok_all &= good
                row += f"  {name}: ({r[0]:+5.1f},{r[1]:+5.1f}){'' if good else ' **'}"
        print(row, flush=True)
    print("  all backends agree with the plant" if ok_all
          else "  ** at least one backend disagrees -- do not trust it below")
    return ok_all


def realcontrol(tag, suffix, n, span, mpp, plants=(0.0, 25.0, -40.0)):
    """The control that actually matters: plant a known shift on the REAL historical
    window and see which backends report it back.

    A synthetic self-test only proves a backend's sign convention. It says nothing
    about whether the backend can match 1967 panchromatic film to modern
    orthoimagery at all -- and AROSICS, which passes every reasonable sanity check
    on satellite data, turns out to be unable to: it returns reliability 0 and does
    not recover a planted 30 m shift on this pair. A backend that cannot recover a
    planted shift here cannot be trusted to tell us the residual is small here
    either, however respectable its provenance."""
    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json'))); bbox = geo['bbox']
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{suffix}.tif'))
    S, W, N, E = bbox; px = int(round(span / mpp))

    def covered(bb):
        c0 = (bb[1]-W)/(E-W)*ds.width; r0 = (N-bb[2])/(N-S)*ds.height
        wp = (bb[3]-bb[1])/(E-W)*ds.width; hp = (bb[2]-bb[0])/(N-S)*ds.height
        a = ds.read(1, window=Window(c0, r0, wp, hp), out_shape=(24, 24),
                    boundless=True, fill_value=0)
        return (a > 0).mean() > 0.97

    wins = naipcheck.windows_over(bbox, covered, n, span)
    print(f"REAL-DATA CONTROL: {tag} ({suffix}) vs NAIP, {len(wins)} windows of "
          f"{span:.0f} m.\n  Each backend must report a planted shift back. "
          f"Scored as |measured change - planted|.\n", flush=True)
    err = {name: [] for name, _ in available()}
    got = {name: 0 for name, _ in available()}
    tried = 0
    for bb in wins:
        nap = naipcheck.fetch(bb, px, px)
        if nap is None or (nap > 0).mean() < 0.9:
            continue
        c0 = (bb[1]-W)/(E-W)*ds.width; r0 = (N-bb[2])/(N-S)*ds.height
        wp = (bb[3]-bb[1])/(E-W)*ds.width; hp = (bb[2]-bb[0])/(N-S)*ds.height
        hist = ds.read(1, window=Window(c0, r0, wp, hp), out_shape=(px, px),
                       boundless=True, fill_value=0)
        if (hist > 0).mean() < 0.9:
            continue
        tried += 1
        for name, fn in available():
            base = None
            try:
                base = fn(hist, nap, mpp)
            except Exception:
                pass
            if base is None:
                continue
            for pE in plants[1:]:
                try:
                    r = fn(gridval.shift_image(hist, pE, 0.0, mpp), nap, mpp)
                except Exception:
                    r = None
                if r is None:
                    continue
                got[name] += 1
                err[name].append(math.hypot(r[0] - base[0] - pE, r[1] - base[1]))
    print(f"  {tried} windows usable\n")
    trust = []
    for name in err:
        v = np.array(err[name]) if err[name] else None
        if v is None or len(v) == 0:
            print(f"  {name:8s}: never returned a comparable pair -- CANNOT VERIFY")
            continue
        good = (v <= 5.0).mean()
        # The median is what says whether a backend can measure this data at all;
        # a minority of windows will always defeat any method (water, land cleared
        # since the flight, a frame edge). Measured on 1967: ours 0.43 m / 71%,
        # scikit-image 1.00 m / 79%, OpenCV 0.27 m / 67% -- all fine. SimpleITK
        # 43.69 m / 0% -- mutual information cannot register this pair.
        ok = good >= 0.6 and np.median(v) <= 5.0
        print(f"  {name:8s}: n={len(v):3d}  median error recovering the plant "
              f"{np.median(v):6.2f} m  within 5 m: {good*100:4.0f}%   "
              f"{'TRUSTWORTHY on this data' if ok else 'NOT trustworthy on this data'}")
        if ok:
            trust.append(name)
    print(f"\n  backends that can measure this data: "
          f"{', '.join(trust) if trust else 'NONE'}")
    return trust


def main():
    if '--realcontrol' in sys.argv:
        realcontrol(arg('--block'), arg('--suffix', 'final'), int(arg('--n', 8)),
                    float(arg('--span', 1500)), float(arg('--mpp', 1.0)))
        return
    if '--selftest' in sys.argv:
        selftest(float(arg('--mpp', 1.0))); return
    tag = arg('--block'); suffix = arg('--suffix', 'final')
    n = int(arg('--n', 16)); span = float(arg('--span', 1500))
    mpp = float(arg('--mpp', 1.0)); px = int(round(span / mpp))
    bound = float(arg('--bound', 60.0))
    ridge = '--ridge' in sys.argv
    selftest(mpp)
    print()
    geo = json.load(open(P('data', f'{tag}_{suffix}_geo.json'))); bbox = geo['bbox']
    ds = rasterio.open(P('mosaics', f'detroit_{tag}_{suffix}.tif'))
    S, W, N, E = bbox

    def covered(bb):
        c0 = (bb[1] - W) / (E - W) * ds.width; r0 = (N - bb[2]) / (N - S) * ds.height
        wpx = (bb[3] - bb[1]) / (E - W) * ds.width
        hpx = (bb[2] - bb[0]) / (N - S) * ds.height
        a = ds.read(1, window=Window(c0, r0, wpx, hpx), out_shape=(24, 24),
                    boundless=True, fill_value=0)
        return (a > 0).mean() > 0.97

    wins = naipcheck.windows_over(bbox, covered, n, span)
    print(f"{tag} ({suffix}) vs USGS NAIP, {len(wins)} windows of {span:.0f} m, "
          f"on {'the road (ridge) response' if ridge else 'raw imagery'}:", flush=True)
    acc = {name: [] for name, _ in available()}
    for k, bb in enumerate(wins):
        nap = naipcheck.fetch(bb, px, px)
        if nap is None or (nap > 0).mean() < 0.9:
            continue
        c0 = (bb[1] - W) / (E - W) * ds.width; r0 = (N - bb[2]) / (N - S) * ds.height
        wpx = (bb[3] - bb[1]) / (E - W) * ds.width
        hpx = (bb[2] - bb[0]) / (N - S) * ds.height
        hist = ds.read(1, window=Window(c0, r0, wpx, hpx), out_shape=(px, px),
                       boundless=True, fill_value=0)
        if (hist > 0).mean() < 0.9:
            continue
        if ridge:
            hist_m, nap_m = _ridge(hist, mpp), _ridge(nap, mpp)
        else:
            hist_m, nap_m = hist, nap
        line = f"  {bb[0]+ (bb[2]-bb[0])/2:.4f},{bb[1]+(bb[3]-bb[1])/2:.4f}"
        for name, fn in available():
            try:
                r = fn(hist_m, nap_m, mpp)
            except Exception:
                r = None
            if r is None:
                line += f"  {name}: --"
            elif math.hypot(*r) > bound:
                line += f"  {name}: >{bound:.0f}"
            else:
                acc[name].append(r)
                line += f"  {name}: {math.hypot(*r):5.1f}"
        print(line, flush=True)
    print()
    for name in acc:
        v = acc[name]
        if not v:
            print(f"  {name:8s}: no result"); continue
        m = np.array([math.hypot(*r) for r in v])
        dE = np.array([r[0] for r in v]); dN = np.array([r[1] for r in v])
        print(f"  {name:8s}: n={len(v):3d}  median {np.median(m):5.2f}  "
              f"p90 {np.percentile(m,90):5.2f}  max {m.max():5.2f} m   "
              f"bias dE {np.median(dE):+5.2f} dN {np.median(dN):+5.2f}")
    json.dump({k: v for k, v in acc.items()}, open(f'/tmp/crosscheck_{tag}.json', 'w'))


if __name__ == '__main__':
    main()
