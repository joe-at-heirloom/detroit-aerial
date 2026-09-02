#!/usr/bin/env python
"""Stage A, second form: close the BLOCK, not just each line.

close.py solves a translation per frame from frame-to-frame observations. That
closes each flight line along its length -- 0.1 m -- and leaves the lines
disagreeing with each other by up to 109 m, because a translation per frame
cannot express what is actually wrong: each line was one pass at one height with
one heading, so each line carries its own scale and rotation, and along-track ties
are blind to a line's scale (every frame in the line shares it). Only sidelap ties
see it, and close.py let 225 of 241 sidelap pairs fall out before the solve.

Unknowns here: a translation per frame, plus a scale and a rotation per line about
that line's centre. Observations: the measured offset between every overlapping
pair, along-track AND cross-line, the latter searched wider because that is where
the disagreement is. For an along-track pair the per-line terms cancel exactly, so
this reduces to close.py; for a cross-line pair the difference in the two lines'
scale and rotation at the overlap point is what the observation constrains.

Relative observations cannot see the block's overall translation, scale or
rotation. Those are pinned to zero here and left for the absolute step -- which
then has a four-parameter problem against a block that is consistent throughout,
instead of a four-parameter problem against four blocks that disagree.

Measured before this existed (crossdiag.py, 1961): 151 of 241 sidelap pairs hit
the search edge at +/-150 m and the 41 that locked disagreed by a median of 114 m,
max 195 m. The lines sit 100-200 m apart. So the cross-line search defaults to
+/-250 m; same-epoch film matches itself strongly enough that the 97.5 m lattice
is not the threat there that it is against satellite imagery, and the IRLS still
downweights any blunder that gets through.

Usage:  ./.venv/bin/python scripts/close2.py 1961 [--rounds 4] [--from stageA] [--cross 250]
"""
import sys, os, json, math, itertools, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust
from validate import P
from rebuild import placements, SP
import close as C
import seamclass as SC


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def pair_observations(rend, sol, lab, mpp, search_along=90.0, search_cross=150.0,
                      min_ratio=1.08, min_ovl_px=150):
    """Measured offset for every overlapping pair, with the overlap midpoint."""
    keys = sorted(rend)
    out = []
    for a, b in itertools.combinations(keys, 2):
        ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
        x0 = max(xa, xb); y0 = max(ya, yb)
        x1 = min(xa + ia.shape[1], xb + ib.shape[1]); y1 = min(ya + ia.shape[0], yb + ib.shape[0])
        if x1 - x0 < min_ovl_px or y1 - y0 < min_ovl_px:
            continue
        A = ia[y0 - ya:y1 - ya, x0 - xa:x1 - xa]; B = ib[y0 - yb:y1 - yb, x0 - xb:x1 - xb]
        if (A > 0).mean() < 0.4 or (B > 0).mean() < 0.4:
            continue
        cross = lab[a] != lab[b]
        ra, ma = gridval.ridge_full(A.astype(np.float32), mpp)
        rb, mb = gridval.ridge_full(B.astype(np.float32), mpp)
        if ra.std() < 1e-9 or rb.std() < 1e-9:
            continue
        m = gridval.match_cell((ra * ma).astype(np.float32), ma.astype(np.float32),
                               (rb * mb).astype(np.float32), mb.astype(np.float32),
                               mpp, search_cross if cross else search_along, 'mncc')
        if m is None or m['edge'] or m['ratio'] < min_ratio:
            continue
        out.append(dict(a=a, b=b, dE=m['dE'], dN=m['dN'], ratio=m['ratio'], cross=cross,
                        qE=(x0 + x1) / 2, qN=(y0 + y1) / 2))   # canvas px; converted by caller
    return out


def _pc(A, B, pad=2):
    """Hann-windowed phase correlation, zero-padded so a large shift is not
    circular. Returns (dx, dy, sharpness, at_boundary) with the same sign convention
    as gridval.match_cell: F(A) * conj(F(B)), positive dx = A sits east of B."""
    h, w = A.shape
    wn = np.outer(np.hanning(h), np.hanning(w)).astype(np.float32)
    a = (A - A.mean()) * wn; b = (B - B.mean()) * wn
    H, W = h * pad, w * pad
    F = np.fft.rfft2(a, s=(H, W)) * np.conj(np.fft.rfft2(b, s=(H, W)))
    F /= (np.abs(F) + 1e-9)
    c = np.fft.fftshift(np.fft.irfft2(F, s=(H, W)))
    pk = np.unravel_index(np.argmax(c), c.shape)
    peak = float(c[pk]); sharp = peak / (float(c.std()) + 1e-12)
    dy = pk[0] - H // 2; dx = pk[1] - W // 2
    at_edge = abs(dx) > w * 0.45 or abs(dy) > h * 0.45
    return dx, dy, sharp, at_edge


def cross_observations(rend, sol, lab, mpp, win_m=384.0, min_valid=0.75, min_sharp=7.0,
                       min_windows=3, agree_m=15.0, max_shift_m=300.0, use_ridge=False,
                       log=print):
    """Sidelap offsets by multi-window phase correlation.

    One big normalised correlation over a sidelap fails at large shifts: with a
    +/-250 m search, 1967's cross-line 'matches' came back at a median of 239 m --
    the search boundary -- because the normalised correlation of the few pixels
    still overlapping at a big shift is high by chance. The solver then fitted
    those blunders (per-line scales of 5%, frames moved 700 m).

    Same-epoch film against itself gives phase correlation a delta peak, which is
    why the original tie-point stage closed to 0.5 m with it. So: tile the sidelap
    into small windows, phase-correlate each, and accept the pair only if several
    windows independently agree on the offset.

    Window size is set by coverage, not by peak strength. Measured on 30 real
    sidelaps (pcdiag.py): every window that reached the peak test passed it with
    median sharpness 15-20 against a threshold of 7, but 78% of windows never got
    there because one frame or the other had under 60% content in them. The frame
    bounding box is sized by its diagonal to fit any rotation, so two sidelapping
    boxes intersect ~2.2 km wide where the film actually overlaps ~0.9 km. 800 m
    windows accepted 3 of 30 sidelaps; 512 m accepted 10; 384 m with a 75% floor
    is the choice. Raw film and the ridge response performed identically. Window-to-window agreement is the
    confidence measure, and a pair whose windows disagree is dropped rather than
    guessed."""
    import itertools
    import gridval as _g
    win_px = max(64, int(round(win_m / mpp)))
    keys = sorted(rend)
    out = []; tried = 0
    feat = {}
    def F(r):
        # the feature phase-correlated: raw film, or the road-ridge response, which
        # survives the shadow change between lines flown hours apart
        if r not in feat:
            img = rend[r][0]
            feat[r] = _g.ridge_full(img.astype(np.float32), mpp)[0] if use_ridge else img
        return feat[r]
    for a, b in itertools.combinations(keys, 2):
        if lab[a] == lab[b]:
            continue
        ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
        fa, fb = F(a), F(b)
        x0 = max(xa, xb); y0 = max(ya, yb)
        x1 = min(xa + ia.shape[1], xb + ib.shape[1]); y1 = min(ya + ia.shape[0], yb + ib.shape[0])
        if x1 - x0 < win_px or y1 - y0 < win_px:
            continue
        tried += 1
        offs = []
        step = win_px // 2
        for wy in range(y0, y1 - win_px + 1, step):
            for wx in range(x0, x1 - win_px + 1, step):
                va = (ia[wy - ya:wy - ya + win_px, wx - xa:wx - xa + win_px] > 0).mean()
                vb = (ib[wy - yb:wy - yb + win_px, wx - xb:wx - xb + win_px] > 0).mean()
                if va < min_valid or vb < min_valid:
                    continue
                A = fa[wy - ya:wy - ya + win_px, wx - xa:wx - xa + win_px].astype(np.float32)
                B = fb[wy - yb:wy - yb + win_px, wx - xb:wx - xb + win_px].astype(np.float32)
                if A.std() < (1e-4 if use_ridge else 5) or B.std() < (1e-4 if use_ridge else 5):
                    continue
                dx, dy, sharp, edge = _pc(A, B)
                if edge or sharp < min_sharp or math.hypot(dx, dy) * mpp > max_shift_m:
                    continue
                offs.append((dx * mpp, -dy * mpp, sharp))
        if len(offs) < min_windows:
            continue
        O = np.array(offs)
        c = np.median(O[:, :2], axis=0)
        for _ in range(3):
            d = np.hypot(O[:, 0] - c[0], O[:, 1] - c[1])
            m = d <= agree_m
            if m.sum() < min_windows:
                break
            c = np.median(O[m, :2], axis=0)
        d = np.hypot(O[:, 0] - c[0], O[:, 1] - c[1]); m = d <= agree_m
        if m.sum() < min_windows:
            continue
        frac = m.sum() / len(O)
        out.append(dict(a=a, b=b, dE=float(c[0]), dN=float(c[1]),
                        ratio=1.0 + 0.5 * frac * min(m.sum(), 8) / 8.0,   # 1.0..1.5
                        cross=True, nwin=int(m.sum()), nwin_all=int(len(O)),
                        qE=(x0 + x1) / 2, qN=(y0 + y1) / 2))
    log(f"           cross-line: {len(out)} pairs from {tried} sidelaps, "
        f"windows agreeing median {np.median([o['nwin'] for o in out]) if out else 0:.0f}")
    return out


def solve_block(sol, obs, recs, lab, minE, maxN, mpp, log=print):
    """Per-frame translation + per-line (scale, rotation), linear least squares,
    iteratively reweighted. Returns corrections in dE-units (position moves by -c)
    and per-line (sigma, rho)."""
    idx = {r: i for i, r in enumerate(recs)}
    nf = len(recs); nl = max(lab.values()) + 1
    # unknown vector: [cE_0..cE_nf-1, cN_0.., sigma_0..sigma_nl-1, rho_0..]
    iE = lambda i: i; iN = lambda i: nf + i; iS = lambda k: 2 * nf + k; iR = lambda k: 2 * nf + nl + k
    nu = 2 * nf + 2 * nl
    pos = {r: np.array([sol[r]['e'] - sol[r]['dE'], sol[r]['n'] - sol[r]['dN']]) for r in recs}
    ctr = {k: np.mean([pos[r] for r in recs if lab[r] == k], axis=0) for k in range(nl)}

    rows = []   # (coeff dict for E-eq, coeff dict for N-eq, rhsE, rhsN, weight, is_cross)
    for o in obs:
        a, b = o['a'], o['b']
        if a not in idx or b not in idx:
            continue
        q = np.array([minE + o['qE'] * mpp, maxN - o['qN'] * mpp])
        ka, kb = lab[a], lab[b]
        eE = {iE(idx[a]): 1.0, iE(idx[b]): -1.0}
        eN = {iN(idx[a]): 1.0, iN(idx[b]): -1.0}
        # content moves by -c + sigma (q-c_k) + rho R90(q-c_k); in dE-units the
        # per-line part enters with a minus sign. R90(v) = (-vN, vE).
        for k, sgn in ((ka, -1.0), (kb, +1.0)):
            v = q - ctr[k]
            eE[iS(k)] = eE.get(iS(k), 0.0) + sgn * v[0]
            eE[iR(k)] = eE.get(iR(k), 0.0) + sgn * (-v[1])
            eN[iS(k)] = eN.get(iS(k), 0.0) + sgn * v[1]
            eN[iR(k)] = eN.get(iR(k), 0.0) + sgn * v[0]
        w = min(max(o['ratio'] - 1.0, 0.02) / 0.5, 3.0)
        rows.append((eE, eN, o['dE'], o['dN'], w, o['cross']))
    if not rows:
        return None

    def build(extra):
        A = np.zeros((nu, nu)); rhs = np.zeros(nu)
        for (eE, eN, rE, rN, w, _), x in zip(rows, extra):
            ww = w * x
            for coeffs, r in ((eE, rE), (eN, rN)):
                ks = list(coeffs.items())
                for i1, v1 in ks:
                    rhs[i1] += ww * v1 * r
                    for i2, v2 in ks:
                        A[i1, i2] += ww * v1 * v2
        # gauge: mean translation, mean scale, mean rotation pinned to zero
        big = 1e3
        for grp in ([iE(i) for i in range(nf)], [iN(i) for i in range(nf)],
                    [iS(k) for k in range(nl)], [iR(k) for k in range(nl)]):
            for i1 in grp:
                for i2 in grp:
                    A[i1, i2] += big / len(grp)
        A += np.eye(nu) * 1e-9
        return np.linalg.solve(A, rhs)

    extra = np.ones(len(rows))
    for it in range(8):
        x = build(extra)
        res = []
        for eE, eN, rE, rN, w, _ in rows:
            pe = sum(v * x[i] for i, v in eE.items()); pn = sum(v * x[i] for i, v in eN.items())
            res.append(math.hypot(pe - rE, pn - rN))
        res = np.array(res)
        sg = max(np.median(res) * 1.4826, 2.0)
        extra = 1.0 / (1.0 + (res / (2.0 * sg)) ** 2)
    cross = np.array([r[5] for r in rows])
    for nm, m in (('along', ~cross), ('cross', cross)):
        if m.any():
            log(f"           fit residual {nm:5s}: n={int(m.sum())}  median {np.median(res[m]):5.1f}  "
                f"p90 {np.percentile(res[m],90):5.1f}  max {res[m].max():5.1f} m   downweighted {(extra[m]<0.25).sum()}")
    cE = {r: float(x[iE(idx[r])]) for r in recs}; cN = {r: float(x[iN(idx[r])]) for r in recs}
    sig = {k: float(x[iS(k)]) for k in range(nl)}; rho = {k: float(x[iR(k)]) for k in range(nl)}
    return cE, cN, sig, rho, ctr


def apply_corrections(sol, recs, lab, cE, cN, sig, rho, ctr):
    for r in recs:
        k = lab[r]
        p = np.array([sol[r]['e'] - sol[r]['dE'], sol[r]['n'] - sol[r]['dN']])
        v = p - ctr[k]
        # centre moves by -c + sigma v + rho R90(v)   (content model above)
        mvE = -cE[r] + sig[k] * v[0] + rho[k] * (-v[1])
        mvN = -cN[r] + sig[k] * v[1] + rho[k] * (v[0])
        sol[r]['dE'] -= mvE; sol[r]['dN'] -= mvN            # position = e - dE
        sol[r]['gw'] *= (1.0 + sig[k]); sol[r]['gh'] *= (1.0 + sig[k])
        sol[r]['rot'] += math.degrees(rho[k])
    return sol


def seam_report(rend, sol, lab, mpp, label, log=print, search_cross=250.0):
    along = [o for o in pair_observations(rend, sol, lab, mpp, search_cross=1.0)
             if not o['cross']]
    cross = cross_observations(rend, sol, lab, mpp, max_shift_m=search_cross, log=log)
    obs = along + cross
    for nm, want in (('along', False), ('cross', True)):
        v = np.array([math.hypot(o['dE'], o['dN']) for o in obs if o['cross'] == want])
        if len(v):
            log(f"  {label} {nm:5s}: {len(v):3d} pairs  median {np.median(v):5.1f}  "
                f"p90 {np.percentile(v,90):5.1f}  max {v.max():5.1f} m  over 10 m {(v>10).sum()}")
        else:
            log(f"  {label} {nm:5s}:   0 pairs")
    return obs


def run(tag, rounds=3, start='stageA', search_cross=250.0):
    t0 = time.time()
    sol, ppm = placements(tag)
    if start and os.path.exists(P('data', f'{start}_{tag}.json')):
        saved = json.load(open(P('data', f'{start}_{tag}.json')))
        for r in list(sol):
            v = saved.get(str(r)) or saved.get(r)
            if v:
                sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
                if 'rot' in v: sol[r]['rot'] = v['rot']
                if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
        print(f"\n===== close2 {tag}: starting from {start} =====", flush=True)
    else:
        print(f"\n===== close2 {tag}: starting from the bundle =====", flush=True)
    recs = sorted([r for r in sol if sol[r].get('ok')])
    lab = SC.lines_of(sol, recs)
    nl = max(lab.values()) + 1
    print(f"  {len(recs)} frames, {nl} flight lines", flush=True)
    mpp = C.MPP_TIE
    imgdir = f"{SP}/fullres"
    for it in range(rounds):
        minE, maxN, W, H = C.canvas(sol, mpp)
        rend = blockadjust.render_all(sol, imgdir, minE, maxN, W, H, mpp, crop=0.95,
                                      log=lambda *_: None)
        obs = seam_report(rend, sol, lab, mpp, f"round {it}", log=print, search_cross=search_cross)
        del rend
        if not obs:
            print("  no observations; stopping"); break
        r = solve_block(sol, obs, recs, lab, minE, maxN, mpp)
        if r is None:
            break
        cE, cN, sig, rho, ctr = r
        for k in range(nl):
            print(f"           line {k}: scale {sig[k]*1e2:+.3f}%  rotation {math.degrees(rho[k]):+.4f} deg", flush=True)
        mv = np.array([math.hypot(cE[r_], cN[r_]) for r_ in recs])
        print(f"           frame translations: median {np.median(mv):.1f}  max {mv.max():.1f} m", flush=True)
        # Physical sanity before anything is applied. A flight line's scale error is
        # a flying-height error; 1.5% is 75 m of altitude on a 5000 m sortie, already
        # generous. A solution asking for more, or moving frames further than twice
        # the disagreement it is trying to remove, is fitting blunders, and applying
        # it would tear the block worse than it found it. Refuse and say so.
        # the disagreement being removed is the CROSS-line one; along-track pairs are
        # already at 0.1 m and would make any move look implausible
        dis = np.array([math.hypot(o['dE'], o['dN']) for o in obs if o['cross']] or [0.0])
        if max(abs(v) for v in sig.values()) > 0.015 or np.median(mv) > 2.0 * np.median(dis) + 20:
            print(f"           REFUSED: implausible solution (max |scale| {max(abs(v) for v in sig.values())*1e2:.2f}%, "
                  f"median move {np.median(mv):.0f} m vs disagreement {np.median(dis):.0f} m); not applied", flush=True)
            break
        sol = apply_corrections(sol, recs, lab, cE, cN, sig, rho, ctr)
        if np.median(mv) < 0.5 and max(abs(v) for v in sig.values()) < 2e-5:
            break
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, imgdir, minE, maxN, W, H, mpp, crop=0.95, log=lambda *_: None)
    seam_report(rend, sol, lab, mpp, "CLOSED", log=print, search_cross=search_cross)
    del rend
    json.dump({r: dict(dE=sol[r]['dE'], dN=sol[r]['dN'], rot=sol[r]['rot'],
                       gw=sol[r]['gw'], gh=sol[r]['gh']) for r in recs},
              open(P('data', f'stageA2_{tag}.json'), 'w'))
    print(f"  saved data/stageA2_{tag}.json  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    run(sys.argv[1], int(arg('--rounds', 3)), arg('--from', 'stageA'), float(arg('--cross', 250.0)))
