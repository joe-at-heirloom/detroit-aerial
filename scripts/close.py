#!/usr/bin/env python
"""Close the block first, place it second. Nothing in between.

This replaces the per-frame correction and the RBF warp, both of which are deleted.
They produced excellent absolute numbers by tearing the mosaic apart: measured pair
by pair, 1949's frames disagreed with each other by 22.6 m out of the bundle and
56.4 m after the per-frame stage "improved" it.

The architecture the evidence supports:

  STAGE A -- close the block.  Solve where each frame sits RELATIVE to its
    neighbours, using only frame-to-frame observations measured now, in map space,
    on the frames themselves. Modern imagery is not allowed to touch an individual
    frame at any point in this stage. Same-epoch frames look alike, often taken
    seconds apart, so this is an easy and trustworthy measurement -- unlike matching
    film to satellite across seventy years, which is where every bad control point
    in this project has come from.

  STAGE B -- place the closed block.  Move the whole thing as ONE rigid body onto
    the map. A handful of parameters against many well-spread observations, so no
    single bad match can bend anything. A rigid transform cannot tear a seam, which
    is the entire point.

Measured on 1961, the disagreement between overlapping frames is about two-thirds a
constant per-pair offset and one-third variation across each overlap. Stage A can
remove the first. The second is per-frame tilt that a shift cannot represent, and it
sets the floor at roughly 8-9 m until there is a camera model. That floor is
published, not hidden.

Usage:  ./.venv/bin/python scripts/close.py 1961 [--dry] [--rigid similarity|shift]
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust, frameadjust, dtmap
import validate
from validate import P, MPP, report
from rebuild import placements, composite, SP
import rasterio

MPP_TIE = 2.0          # frame-to-frame work; same epoch, so resolution is cheap here


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def canvas(sol, mpp):
    recs = [r for r in sol if sol[r].get('ok')]
    E = [sol[r]['e'] - sol[r].get('dE', 0.0) for r in recs]
    N = [sol[r]['n'] - sol[r].get('dN', 0.0) for r in recs]
    g = max(math.hypot(sol[r]['gw'], sol[r]['gh']) / 2 for r in recs)
    minE = min(E) - g; maxN = max(N) + g
    return minE, maxN, int((max(E) + g - minE) / mpp), int((maxN - (min(N) - g)) / mpp)


def seam_stats(rend, mpp):
    obs = blockadjust.relative_observations(rend, mpp, search_m=90.0,
                                            min_ratio=1.10, log=lambda *_: None)
    if not obs:
        return None, []
    v = np.array([math.hypot(o[2], o[3]) for o in obs])
    return v, obs


def close_block(sol, imgdir, mpp=MPP_TIE, rounds=3, log=print):
    """Stage A: minimise frame-to-frame disagreement. No modern imagery involved."""
    for it in range(rounds):
        minE, maxN, W, H = canvas(sol, mpp)
        rend = blockadjust.render_all(sol, imgdir, minE, maxN, W, H, mpp,
                                      crop=0.95, log=lambda *_: None)
        v, obs = seam_stats(rend, mpp)
        del rend
        if v is None:
            log("    no overlap observations; stopping"); return sol
        log(f"  round {it}: {len(obs)} overlapping pairs, frames disagree by "
            f"median {np.median(v):5.1f}  p90 {np.percentile(v,90):5.1f}  "
            f"max {v.max():5.1f} m", flush=True)
        if np.median(v) < 1.0:
            break
        recs = sorted([r for r in sol if sol[r].get('ok')])
        idx = {r: i for i, r in enumerate(recs)}
        n = len(recs)

        use = [(a_, b_, de, dn, rt) for a_, b_, de, dn, rt in obs
               if a_ in idx and b_ in idx]

        def solve(extra):
            """extra: per-observation robustness multiplier."""
            A = np.zeros((n, n)); bE = np.zeros(n); bN = np.zeros(n)
            for k, (a_, b_, de, dn, rt) in enumerate(use):
                w = min(max(rt - 1.0, 0.02) / 0.5, 3.0) * extra[k]
                ia, ib = idx[a_], idx[b_]
                A[ia, ia] += w; A[ib, ib] += w
                A[ia, ib] -= w; A[ib, ia] -= w
                bE[ia] += w * de; bE[ib] -= w * de
                bN[ia] += w * dn; bN[ib] -= w * dn
            # gauge: relative observations fix only the SHAPE of the block, so pin
            # the mean correction to zero and let stage B decide where it all goes
            A += np.ones((n, n)) * 1e-3
            A += np.eye(n) * 1e-9
            return np.linalg.solve(A, bE), np.linalg.solve(A, bN)

        # Iteratively reweighted: a pair that disagrees with the emerging consensus
        # is downweighted rather than allowed to drag its two frames. Without this,
        # six one-block blunders among sixty pairs put 14 m of error into the
        # solution -- and this scene aliases at 97.5 m, so blunders are the norm.
        extra = np.ones(len(use))
        for _ in range(6):
            cE, cN = solve(extra)
            res = np.array([math.hypot((cE[idx[a_]] - cE[idx[b_]]) - de,
                                       (cN[idx[a_]] - cN[idx[b_]]) - dn)
                            for a_, b_, de, dn, rt in use])
            sg = max(np.median(res) * 1.4826, 2.0)
            extra = 1.0 / (1.0 + (res / (2.0 * sg)) ** 2)
        rej = int((extra < 0.25).sum())
        log(f"           {len(use)} observations, {rej} downweighted as blunders; "
            f"consensus residual median {np.median(res):.1f} m", flush=True)
        cE -= cE.mean(); cN -= cN.mean()
        for r in recs:
            sol[r]['dE'] += float(cE[idx[r]]); sol[r]['dN'] += float(cN[idx[r]])
        mv = np.hypot(cE, cN)
        log(f"           corrections applied: median {np.median(mv):5.1f}  "
            f"max {mv.max():5.1f} m", flush=True)
    # final measurement
    minE, maxN, W, H = canvas(sol, mpp)
    rend = blockadjust.render_all(sol, imgdir, minE, maxN, W, H, mpp, crop=0.95,
                                  log=lambda *_: None)
    v, obs = seam_stats(rend, mpp)
    del rend
    if v is not None:
        log(f"  CLOSED: {len(obs)} pairs, disagreement median {np.median(v):5.1f}  "
            f"p90 {np.percentile(v,90):5.1f}  max {v.max():5.1f} m", flush=True)
    return sol


def place_block(sol, imgdir, t, mpp, grid, kind='similarity', log=print):
    """Stage B: one rigid move of the whole block onto the map.

    `grid` is the (minE, maxN) of the raster the modern ridge maps in `t` were
    built on. The frames MUST be rendered into that same frame of reference --
    rendering them into a canvas derived from the current solution instead compares
    them against the wrong part of the modern imagery, which is a silent and
    total failure."""
    minE, maxN = grid
    H, W = t['shape']
    # A prior is NOT optional here. The bundle's absolute error is tens of metres and
    # Detroit's street grid repeats every 97.5 m, so a bare +/-40 m fine search locks
    # onto the nearest alias and reports a confident, tiny, wrong shift -- which is
    # exactly what happened: 3.6 m reported where ~66 m was needed. The regional
    # coarse field is solved on 6 km windows of mile-grid arterials, which cannot
    # alias, and prior_for refuses it if too few windows lock.
    pri = gridval.prior_for(t, mpp, log=log)
    abs_obs, _ = blockadjust.observe(sol, imgdir, t, minE, maxN, mpp,
                                     prior=pri, log=lambda *_: None)
    if len(abs_obs) < 6:
        log(f"    only {len(abs_obs)} frames matched modern; not placing"); return sol
    P_ = np.array([[sol[r]['e'] - sol[r].get('dE', 0.0),
                    sol[r]['n'] - sol[r].get('dN', 0.0)] for r in abs_obs])
    D = np.array([[o['dE'], o['dN']] for o in abs_obs.values()])
    w = np.array([min(max(o['ratio'] - 1.0, 0.02) / 0.5, 3.0) for o in abs_obs.values()])
    # robust: trim observations that disagree with the consensus, iteratively
    keep = np.ones(len(D), bool)
    for _ in range(6):
        c = np.average(D[keep], axis=0, weights=w[keep])
        r = np.hypot(D[:, 0] - c[0], D[:, 1] - c[1])
        s = max(np.median(r[keep]) * 1.4826, 5.0)
        keep = r < 3.0 * s
        if keep.sum() < 6:
            keep = r < np.percentile(r, 75); break
    log(f"    {len(abs_obs)} frames matched modern, {int(keep.sum())} kept after "
        f"outlier trim; shift {c[0]:+.1f}, {c[1]:+.1f} m", flush=True)
    if kind == 'shift':
        for r in sol:
            if sol[r].get('ok'):
                sol[r]['dE'] += float(c[0]); sol[r]['dN'] += float(c[1])
        return sol
    # similarity: shift plus a small rotation and scale about the block centre
    Q = P_[keep]; V = D[keep]; ww = w[keep]
    ctr = np.average(Q, axis=0, weights=ww)
    X = Q - ctr
    M = np.zeros((4, 4)); rhs = np.zeros(4)
    for (x, y), (de, dn), wi in zip(X, V, ww):
        for row, val in ((np.array([1, 0, x, -y]), de), (np.array([0, 1, y, x]), dn)):
            M += wi * np.outer(row, row); rhs += wi * row * val
    try:
        p = np.linalg.solve(M + np.eye(4) * 1e-9, rhs)
    except np.linalg.LinAlgError:
        p = np.array([c[0], c[1], 0.0, 0.0])
    tx, ty, a, b = p
    log(f"    similarity: shift {tx:+.1f}, {ty:+.1f} m; scale {1+a:.6f}; "
        f"rotation {math.degrees(math.atan2(b, 1+a)):+.4f} deg", flush=True)
    # A similarity moves the frames' POSITIONS and also rotates and scales the
    # frames themselves. Applying it to centres alone shears the block apart --
    # measured, that took internal consistency from 0.2 m back to 4.4 m.
    scale = math.hypot(1.0 + a, b)
    rot_deg = math.degrees(math.atan2(b, 1.0 + a))
    # D(p) is where the content SITS relative to truth: rotated by +theta, scaled
    # by s. The correction is the inverse. Positions move by -D (which is what
    # `dE += D` does, since position = e - dE), and each frame must therefore be
    # un-rotated by theta and its footprint divided by s. The previous version
    # added theta and multiplied by s -- doubling the error instead of removing it,
    # which at 0.77% over a 3.3 km frame is ~25 m in every overlap. That, not
    # "scale absorbing error", is why consistency went 0.2 -> 21.2 m.
    for r in sol:
        if not sol[r].get('ok'):
            continue
        x = (sol[r]['e'] - sol[r].get('dE', 0.0)) - ctr[0]
        y = (sol[r]['n'] - sol[r].get('dN', 0.0)) - ctr[1]
        sol[r]['dE'] += float(tx + a * x - b * y)
        sol[r]['dN'] += float(ty + b * x + a * y)
        sol[r]['rot'] = float(sol[r]['rot'] - rot_deg)
        sol[r]['gw'] = float(sol[r]['gw'] / scale)
        sol[r]['gh'] = float(sol[r]['gh'] / scale)
    return sol


def run(tag, dry=False, kind='similarity'):
    t0 = time.time()
    print(f"\n===== close {tag} =====", flush=True)
    sol, ppm = placements(tag)
    imgdir = f"{SP}/fullres"
    print("STAGE A -- close the block (frame-to-frame only, no modern imagery)", flush=True)
    sol = close_block(sol, imgdir)
    # save Stage A on its own: it is the expensive part and stage B is cheap to redo
    json.dump({r: dict(dE=sol[r]['dE'], dN=sol[r]['dN'], rot=sol[r]['rot'],
                       gw=sol[r]['gw'], gh=sol[r]['gh'])
               for r in sol if sol[r].get('ok')},
              open(P('data', f'stageA_{tag}.json'), 'w'))
    print("STAGE B -- place the closed block as one rigid body", flush=True)
    arr, mod, bbox = validate.load(tag, 'pre')
    t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_pre_modern', log=lambda *_: None)
    grid = ((bbox[1] - dtmap.LON0) * dtmap.MLON, (bbox[2] - dtmap.LAT0) * dtmap.MLAT)
    sol = place_block(sol, imgdir, t, MPP, grid, kind=kind)
    json.dump({r: dict(dE=sol[r]['dE'], dN=sol[r]['dN'], rot=sol[r]['rot'],
                       gw=sol[r]['gw'], gh=sol[r]['gh'])
               for r in sol if sol[r].get('ok')},
              open(P('data', f'closed_{tag}.json'), 'w'))
    if dry:
        print(f"  [dry] {time.time()-t0:.0f}s"); return
    geo = composite(tag, sol, ppm, 'closed')
    arr2, mod2, _ = validate.load(tag, 'closed')
    t2 = gridval.prepare_cached(arr2, mod2, MPP, f'{tag}_closed_{int(t0)}')
    pz = gridval.Prior([], MPP)
    print("  absolute accuracy of the closed, placed block:", flush=True)
    report(gridval.grid_hier_prep(t2, prior=pz), f'  {tag} 16x3')
    report(gridval.grid_hier_prep(t2, NY=32, NX=6, prior=pz), f'  {tag} 32x6')
    print(f"  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    dry = '--dry' in sys.argv
    kind = arg('--rigid', 'similarity')
    for tg in [a for a in sys.argv[1:] if not a.startswith('--') and a != kind]:
        run(tg, dry, kind)
