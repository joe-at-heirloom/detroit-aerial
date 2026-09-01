"""Per-frame corrections that keep the block in one piece.

The per-frame stage as first written moved every negative independently to
best-match modern imagery. That fixes each frame's absolute position and destroys
the block's internal geometry: on 1949 the corrections had a median of 90 m but a
p90 of 196 m and a max of 312 m, and that ~200 m of *spread* is relative error
injected into a block whose bundle adjustment had made it self-consistent to
0.5-1.7 m. A road crossing the seam between two such frames visibly kinks, which
no amount of "the area is in the right place on average" makes acceptable -- and
which the grid metric cannot see, because a 3.6 x 2.0 km cell contains several
frames and averages their disagreement away.

The two things being asked for are in tension and have very different precisions:

  relative  where each frame sits with respect to its neighbours. The bundle
            already knows this to about a metre.
  absolute  where the block as a whole sits on the map. Known only from matching
            film against modern imagery, which is noisy, and sometimes wrong.

So solve them together instead of letting the absolute observations overwrite the
relative ones. For corrections c_i, minimise

    sum_i  w_i |c_i - d_i|^2   +   sum_(i,j) v_ij |c_i - c_j|^2

the first term pulling each frame toward its own measurement, the second holding
neighbours together. With v >> w the result is a correction that varies smoothly
along the strip -- which is what bundle drift actually looks like over a 30 km
flight line -- while refusing the frame-to-frame jumps that tear the mosaic.

That is a graph-Laplacian-regularised least squares, and it is small: a few dozen
frames, one linear solve per axis.
"""
import math
import numpy as np


def neighbours(sol, recs, frac=0.75):
    """Pairs of frames that actually overlap, with a weight that falls off as the
    overlap shrinks -- a pair sharing a sliver constrains less than one sharing
    half its area."""
    out = []
    for i, a in enumerate(recs):
        for b in recs[i + 1:]:
            ra = math.hypot(sol[a]['gw'], sol[a]['gh']) / 2
            rb = math.hypot(sol[b]['gw'], sol[b]['gh']) / 2
            d = math.hypot(sol[a]['e'] - sol[b]['e'], sol[a]['n'] - sol[b]['n'])
            lim = frac * (ra + rb)
            if 1.0 < d < lim:
                out.append((a, b, float(max(0.05, 1.0 - d / lim))))
    return out


def _solve_axis(recs, d, w, pairs, ridge=1e-6):
    n = len(recs)
    idx = {r: i for i, r in enumerate(recs)}
    A = np.zeros((n, n)); b = np.zeros(n)
    for i in range(n):
        A[i, i] += w[i]; b[i] += w[i] * d[i]
    for a, bb, v in pairs:
        ia, ib = idx[a], idx[bb]
        A[ia, ia] += v; A[ib, ib] += v
        A[ia, ib] -= v; A[ib, ia] -= v
    A += np.eye(n) * ridge
    return np.linalg.solve(A, b)


def solve(sol, obs, sigma_abs=20.0, sigma_rel=2.0, sigma_rot_abs=0.30,
          sigma_rot_rel=0.05, log=print):
    """obs: {rec: dict(dE, dN, rot, ratio)} -- what each frame measured on its own.

    Returns {rec: dict(dE, dN, rot)} for EVERY frame in sol, including frames that
    never matched: with the block held together they inherit a position from their
    neighbours automatically, rather than needing a separate borrowing rule.

    sigma_rel is the precision of the bundle's relative geometry; sigma_abs the
    precision of a single frame's match against modern imagery. The ratio between
    them is the whole design -- it decides how much the block is allowed to deform
    to chase the absolute measurements."""
    recs = sorted([r for r in sol if sol[r].get('ok')])
    if not recs:
        return {}
    pr = neighbours(sol, recs)
    wa = 1.0 / sigma_abs ** 2
    vr = 1.0 / sigma_rel ** 2
    pairs = [(a, b, vr * f) for a, b, f in pr]

    dE = np.zeros(len(recs)); dN = np.zeros(len(recs)); dR = np.zeros(len(recs))
    w = np.zeros(len(recs)); wr = np.zeros(len(recs))
    for i, r in enumerate(recs):
        o = obs.get(r)
        if not o:
            continue
        # a confident match is worth more, but never more than a few times a weak one
        conf = min(max(o.get('ratio', 1.0) - 1.0, 0.02) / 0.5, 3.0)
        dE[i] = o['dE']; dN[i] = o['dN']; dR[i] = o.get('rot', 0.0)
        w[i] = wa * conf
        wr[i] = (1.0 / sigma_rot_abs ** 2) * conf
    cE = _solve_axis(recs, dE, w, pairs)
    cN = _solve_axis(recs, dN, w, pairs)
    rpairs = [(a, b, (1.0 / sigma_rot_rel ** 2) * f) for a, b, f in pr]
    cR = _solve_axis(recs, dR, wr, rpairs)

    meas = np.array([1.0 if obs.get(r) else 0.0 for r in recs])
    sh = np.hypot(dE - cE, dN - cN)[meas > 0]
    spread_in = np.hypot(dE, dN)[meas > 0]
    log(f"    joint adjustment: {len(recs)} frames, {len(pairs)} overlapping pairs, "
        f"{int(meas.sum())} measured")
    if len(sh):
        log(f"      per-frame measurements had spread {np.percentile(spread_in,90)-np.median(spread_in):5.1f} m "
            f"(median {np.median(spread_in):5.1f}, p90 {np.percentile(spread_in,90):5.1f}); "
            f"the adjustment moved them {np.median(sh):5.1f} m to hold the block together")
    return {r: dict(dE=float(cE[i]), dN=float(cN[i]), rot=float(cR[i]))
            for i, r in enumerate(recs)}


# --- measured relative geometry ----------------------------------------------
#
# The formulation above assumes the bundle's relative geometry is right and only
# the block's absolute placement is in doubt. That holds for 1961 and 1967, whose
# bundles close to 0.5-0.8 m. It does NOT hold for 1949 and 1956: measured pair by
# pair, the 1949 bundle's own frames disagree with each other by a median of
# 22.6 m, p90 60 m, before anything has been corrected. Anchoring to that is
# anchoring to the problem.
#
# So measure the relative geometry too, from the imagery, in map space, right now
# -- the same quantity the original tie-point stage estimated, but re-observed on
# the frames as they are currently placed. Then the adjustment has two sets of real
# observations and can satisfy both: neighbours agree with each other AND the block
# sits on the map.
#
# Each frame is rendered once and reused for every pair it belongs to; rendering
# per pair would repeat the expensive part hundreds of times.

def render_all(sol, imgdir, minE, maxN, W, H, mpp, crop=0.92, log=print):
    import frameadjust
    out = {}
    recs = sorted([r for r in sol if sol[r].get('ok')])
    for i, r in enumerate(recs):
        v = frameadjust.render_frame(r, sol[r], imgdir, minE, maxN, W, H, mpp, crop=crop)
        if v is None:
            continue
        img, x0, y0 = v
        if (img > 0).mean() < 0.2:
            continue
        rg, vd = None, None
        out[r] = (img, x0, y0)
        if (i + 1) % 20 == 0:
            log(f"    rendered {i+1}/{len(recs)} frames")
    return out


def relative_observations(rend, mpp, search_m=90.0, min_ratio=1.12, log=print):
    """Measured offset between every overlapping pair, from the imagery itself."""
    import gridval
    keys = sorted(rend)
    obs = []
    for i, a in enumerate(keys):
        ia, xa, ya = rend[a]
        for b in keys[i + 1:]:
            ib, xb, yb = rend[b]
            x0 = max(xa, xb); y0 = max(ya, yb)
            x1 = min(xa + ia.shape[1], xb + ib.shape[1])
            y1 = min(ya + ia.shape[0], yb + ib.shape[0])
            if x1 - x0 < 180 or y1 - y0 < 180:
                continue
            A = ia[y0 - ya:y1 - ya, x0 - xa:x1 - xa]
            B = ib[y0 - yb:y1 - yb, x0 - xb:x1 - xb]
            if (A > 0).mean() < 0.5 or (B > 0).mean() < 0.5:
                continue
            ra, va = gridval.ridge_full(A.astype(np.float32), mpp)
            rb, vb = gridval.ridge_full(B.astype(np.float32), mpp)
            if ra.std() < 1e-9 or rb.std() < 1e-9:
                continue
            m = gridval.match_cell((ra * va).astype(np.float32), va.astype(np.float32),
                                   (rb * vb).astype(np.float32), vb.astype(np.float32),
                                   mpp, search_m, 'mncc')
            if m is None or m['edge'] or m['ratio'] < min_ratio:
                continue
            # A sits (dE, dN) from B: correction_a - correction_b should equal that
            obs.append((a, b, m['dE'], m['dN'], float(m['ratio'])))
    if obs:
        v = np.array([math.hypot(o[2], o[3]) for o in obs])
        log(f"    relative geometry: {len(obs)} pairs measured, neighbours currently "
            f"disagree by median {np.median(v):.1f} m (p90 {np.percentile(v,90):.1f}, "
            f"max {v.max():.1f})")
    return obs


def solve_full(sol, abs_obs, rel_obs, sigma_abs=18.0, sigma_rel=3.0, log=print):
    """Corrections satisfying measured relative AND absolute observations.

    This is a bundle adjustment with two observation types, which is what the
    problem actually is. The relative observations say how the frames must sit with
    respect to each other; the absolute ones say where the whole thing belongs. The
    sigmas set which wins where they conflict -- relative tighter, because a
    frame-to-frame match on identical film is a far better measurement than film
    against modern satellite."""
    recs = sorted([r for r in sol if sol[r].get('ok')])
    idx = {r: i for i, r in enumerate(recs)}
    n = len(recs)
    if n == 0:
        return {}
    wa = 1.0 / sigma_abs ** 2
    vr = 1.0 / sigma_rel ** 2

    def axis(dabs, wabs, rel):
        A = np.zeros((n, n)); b = np.zeros(n)
        for i in range(n):
            A[i, i] += wabs[i]; b[i] += wabs[i] * dabs[i]
        for a_, b_, r_, w_ in rel:
            ia, ib = idx[a_], idx[b_]
            A[ia, ia] += w_; A[ib, ib] += w_
            A[ia, ib] -= w_; A[ib, ia] -= w_
            b[ia] += w_ * r_; b[ib] -= w_ * r_
        A += np.eye(n) * 1e-6
        return np.linalg.solve(A, b)

    dE = np.zeros(n); dN = np.zeros(n); wab = np.zeros(n)
    for r, o in abs_obs.items():
        if r not in idx:
            continue
        i = idx[r]
        conf = min(max(o.get('ratio', 1.0) - 1.0, 0.02) / 0.5, 3.0)
        dE[i] = o['dE']; dN[i] = o['dN']; wab[i] = wa * conf
    relE = [(a, b, de, vr * min(max(rt - 1.0, 0.02) / 0.5, 3.0))
            for a, b, de, dn, rt in rel_obs]
    relN = [(a, b, dn, vr * min(max(rt - 1.0, 0.02) / 0.5, 3.0))
            for a, b, de, dn, rt in rel_obs]
    cE = axis(dE, wab, relE); cN = axis(dN, wab, relN)
    # residuals, so the caller can see whether the two observation sets could be
    # reconciled or were fighting
    if rel_obs:
        rr = [math.hypot((cE[idx[a]] - cE[idx[b]]) - de,
                         (cN[idx[a]] - cN[idx[b]]) - dn)
              for a, b, de, dn, rt in rel_obs]
        log(f"    after adjustment: relative residual median {np.median(rr):.1f} m "
            f"(p90 {np.percentile(rr,90):.1f})")
    if abs_obs:
        ar = [math.hypot(cE[idx[r]] - o['dE'], cN[idx[r]] - o['dN'])
              for r, o in abs_obs.items() if r in idx]
        log(f"    after adjustment: absolute residual median {np.median(ar):.1f} m "
            f"(p90 {np.percentile(ar,90):.1f})")
    return {r: dict(dE=float(cE[i]), dN=float(cN[i]), rot=0.0)
            for i, r in enumerate(recs)}


def observe(sol, imgdir, t, minE, maxN, mpp, prior=None, crop=0.92,
            abs_search=None, rel_search=90.0, log=print):
    """Render every frame once, then take both kinds of observation off it.

    Rendering is the expensive step and both observation sets need the same
    rendered frames, so doing them together costs little more than doing either
    alone."""
    import gridval, frameadjust
    H, W = t['shape']
    if abs_search is None:
        abs_search = gridval.FINE_SEARCH
    rend = render_all(sol, imgdir, minE, maxN, W, H, mpp, crop=crop, log=log)
    log(f"    rendered {len(rend)} frames")

    abs_obs = {}
    for r, (img, x0, y0) in rend.items():
        y1 = y0 + img.shape[0]; x1 = x0 + img.shape[1]
        rf, vf = gridval.ridge_full(img.astype(np.float32), mpp)
        if rf.std() < 1e-9:
            continue
        pr = prior.at((y0 + y1) / 2, (x0 + x1) / 2) if prior is not None else (0.0, 0.0)
        m = frameadjust._match_local(rf, vf, t, y0, y1, x0, x1, mpp, pr)
        if m is None or m['pegged'] or m['ratio'] < 1.10:
            continue
        abs_obs[r] = m
    if abs_obs:
        v = np.array([math.hypot(o['dE'], o['dN']) for o in abs_obs.values()])
        log(f"    absolute: {len(abs_obs)}/{len(rend)} frames matched modern, "
            f"median {np.median(v):.1f} m (p90 {np.percentile(v,90):.1f})")
    rel_obs = relative_observations(rend, mpp, search_m=rel_search, log=log)
    del rend
    return abs_obs, rel_obs
