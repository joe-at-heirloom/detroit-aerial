#!/usr/bin/env python
"""Stage A, third form: a tie-point bundle adjustment with a similarity per frame.

Why the second form failed, measured: a scale and rotation per flight LINE, fitted
to one offset per overlapping pair, took 1961's cross-line disagreement from 63 to
21 m by pushing along-track from 0.2 to 9 m. Least squares sacrificed 50 good
observations to serve 121 bad-fitting ones because the model could not satisfy
both -- the lines are not rigidly scaled and rotated relative to each other.

The structural mistake was upstream of the model: one offset per pair. A single
offset determines a translation and nothing else; scale and rotation are visible
only as the offset VARYING across the overlap, and collapsing the windows to a
median throws exactly that away. So here every phase-correlation window is its
own tie point, with its own location, for along-track and cross-line pairs alike
-- thousands of them -- and the unknowns are a translation, a scale and a rotation
for every frame. Along-track windows pin neighbouring frames' scales together
(0.2 m agreement over a 2 km overlap allows 0.03% of difference); cross-line
windows pin the lines to each other. Nothing is averaged before the solve.

The block's overall translation, scale and rotation are unobservable from ties and
are pinned to zero for the absolute stage. A mild prior holds each frame's scale
and rotation near zero so a frame with few windows is not free to wander.

Usage:  ./.venv/bin/python scripts/close3.py 1961 [--rounds 3] [--from stageA]
        [--model similarity|affine|translate]
"""
import sys, os, json, math, itertools, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import blockadjust
from validate import P
from rebuild import placements, SP
import close as C
import close2
import seamclass as SC


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def content_box(img, x0, y0):
    rows = np.where((img > 0).mean(1) > 0.5)[0]; cols = np.where((img > 0).mean(0) > 0.5)[0]
    if len(rows) == 0 or len(cols) == 0:
        return None
    return (x0 + cols[0], y0 + rows[0], x0 + cols[-1] + 1, y0 + rows[-1] + 1)


def tie_windows(rend, lab, mpp, win_m=384.0, step_frac=0.5, min_valid=0.75, min_sharp=7.0,
                max_shift_m=300.0, agree_m=20.0, min_windows=3, feature='raw', log=print):
    """Every window of every overlapping pair, as its own tie point.

    feature='ridge' phase-correlates the road-ridge response instead of the raw
    film. On 1949 -- grainy film, three lines, open land -- raw-texture windows over
    fields returned random peaks that scattered by +/-45 m within a pair, and the
    consensus filter then left the worst sidelaps with 2-3 windows, i.e. no
    constraint. Roads persist where fields do not; a ridge window with no road in
    it is flat and skipped rather than wrong."""
    import gridval
    win = max(64, int(round(win_m / mpp))); step = max(16, int(win * step_frac))
    keys = sorted(rend)
    boxes = {r: content_box(rend[r][0], rend[r][1], rend[r][2]) for r in keys}
    feat = {}
    if feature == 'ridge':
        for r in keys:
            feat[r] = gridval.ridge_full(rend[r][0].astype(np.float32), mpp)[0]
    flat_thr = 1e-4 if feature == 'ridge' else 5.0
    ties = []; npairs = {'along': 0, 'cross': 0}; rejected = 0
    for a, b in itertools.combinations(keys, 2):
        ta, tb = boxes[a], boxes[b]
        if ta is None or tb is None:
            continue
        x0 = max(ta[0], tb[0]); y0 = max(ta[1], tb[1]); x1 = min(ta[2], tb[2]); y1 = min(ta[3], tb[3])
        if x1 - x0 < win or y1 - y0 < win:
            continue
        ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
        cross = lab[a] != lab[b]
        cand = []
        for wy in range(y0, y1 - win + 1, step):
            for wx in range(x0, x1 - win + 1, step):
                A = ia[wy - ya:wy - ya + win, wx - xa:wx - xa + win]
                B = ib[wy - yb:wy - yb + win, wx - xb:wx - xb + win]
                if (A > 0).mean() < min_valid or (B > 0).mean() < min_valid:
                    continue
                if feature == 'ridge':
                    A = feat[a][wy - ya:wy - ya + win, wx - xa:wx - xa + win]
                    B = feat[b][wy - yb:wy - yb + win, wx - xb:wx - xb + win]
                A = A.astype(np.float32); B = B.astype(np.float32)
                if A.std() < flat_thr or B.std() < flat_thr:
                    continue
                dx, dy, sharp, edge = close2._pc(A, B)
                if edge or sharp < min_sharp or math.hypot(dx, dy) * mpp > max_shift_m:
                    continue
                cand.append(dict(a=a, b=b, qx=wx + win / 2, qy=wy + win / 2,
                                 dE=dx * mpp, dN=-dy * mpp, sharp=float(sharp), cross=cross))
        # Coarse-to-fine, per pair. A 384 m window on a 97.5 m street lattice can
        # lock a period away and still be sharp; measured, along-track windows
        # disagreed by 15 m median (p90 71) inside overlaps that agree to 0.2 m as
        # a whole. So the pair's robust offset is the prior, and a window is a tie
        # only if it agrees with it within `agree_m` -- which keeps the genuine
        # +/-10 m variation across the overlap that scale and rotation live in.
        if len(cand) < min_windows:
            rejected += len(cand); continue
        O = np.array([[c['dE'], c['dN']] for c in cand])
        ctr = np.median(O, axis=0)
        for _ in range(4):
            d = np.hypot(O[:, 0] - ctr[0], O[:, 1] - ctr[1]); m = d <= agree_m
            if m.sum() < min_windows:
                break
            ctr = np.median(O[m], axis=0)
        d = np.hypot(O[:, 0] - ctr[0], O[:, 1] - ctr[1]); m = d <= agree_m
        if m.sum() < min_windows:
            rejected += len(cand); continue
        rejected += int((~m).sum())
        ties.extend(c for c, keep in zip(cand, m) if keep)
        npairs['cross' if cross else 'along'] += 1
    log(f"           ties: {sum(1 for t in ties if not t['cross'])} along-track windows on "
        f"{npairs['along']} pairs, {sum(1 for t in ties if t['cross'])} cross-line windows on "
        f"{npairs['cross']} pairs; {rejected} windows rejected as off the pair's consensus")
    return ties


def seam_stats(ties, label, log=print):
    out = {}
    for nm, want in (('along', False), ('cross', True)):
        v = np.array([math.hypot(t['dE'], t['dN']) for t in ties if t['cross'] == want])
        if len(v):
            log(f"  {label} {nm:5s}: {len(v):5d} windows  median {np.median(v):5.1f}  "
                f"p90 {np.percentile(v,90):5.1f}  max {v.max():5.1f} m")
            out[nm] = float(np.median(v))
        else:
            out[nm] = None
    return out


def solve(sol, ties, recs, minE, maxN, mpp, model='similarity', prior_m=0.25, log=print):
    """prior_m sets how strongly each frame's scale and rotation are held near zero:
    a scale of 1% costs the same as a (10/prior_m)^2 * 0.25 m^2 residual. At 0.25, 1%
    costs like a 20 m residual on one tie -- enough that a frame with few ties cannot
    wander to 4.5% (which happened at 1.0), not enough to suppress the ~0.5% that is
    physically there and supported by hundreds of ties."""
    """Per-frame [cE, cN, sigma, rho] (affine adds a shear pair). Returns dict of
    corrections. Sign convention as close2: content of frame i at q moves by
    -c_i + sigma_i (q-p_i) + rho_i R90(q-p_i); a tie 'a sits d from b' requires
    c_a - c_b - S_a(q) + S_b(q) = d  where S is the per-frame linear part."""
    idx = {r: i for i, r in enumerate(recs)}
    nf = len(recs)
    npar = {'translate': 2, 'similarity': 4, 'affine': 6}[model]
    nu = nf * npar
    pos = {r: np.array([sol[r]['e'] - sol[r]['dE'], sol[r]['n'] - sol[r]['dN']]) for r in recs}
    L = 1000.0    # lever-arm scale so all unknowns are O(1): sigma*L, rho*L are metres/km

    def lin_coeffs(i, v):
        """coefficients of frame i's unknowns in the E and N equations, for lever v."""
        cE = {i * npar + 0: 1.0}; cN = {i * npar + 1: 1.0}
        if npar >= 4:
            # -sigma*v - rho*R90(v)   with unknowns scaled by L
            cE[i * npar + 2] = -v[0] / L; cE[i * npar + 3] = -(-v[1]) / L
            cN[i * npar + 2] = -v[1] / L; cN[i * npar + 3] = -(v[0]) / L
        if npar == 6:
            # extra affine terms: dE += s1*vE + s2*vN ; dN += s3*vE - s1*vN (traceless shear)
            cE[i * npar + 4] = -v[0] / L; cE[i * npar + 5] = -v[1] / L
            cN[i * npar + 4] = +v[1] / L; cN[i * npar + 5] = -v[0] / L
        return cE, cN

    rows = []
    for t in ties:
        a, b = t['a'], t['b']
        if a not in idx or b not in idx:
            continue
        q = np.array([minE + t['qx'] * mpp, maxN - t['qy'] * mpp])
        eE = {}; eN = {}
        for r, sgn in ((a, 1.0), (b, -1.0)):
            cE, cN = lin_coeffs(idx[r], q - pos[r])
            for k, v in cE.items(): eE[k] = eE.get(k, 0.0) + sgn * v
            for k, v in cN.items(): eN[k] = eN.get(k, 0.0) + sgn * v
        w = min(t['sharp'] / 10.0, 3.0)
        rows.append((eE, eN, t['dE'], t['dN'], w, t['cross']))
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
        big = 1e3
        for j in range(npar):                       # gauge: mean of every parameter = 0
            grp = [i * npar + j for i in range(nf)]
            for i1 in grp:
                for i2 in grp:
                    A[i1, i2] += big / len(grp)
        for i in range(nf):                         # mild prior on the linear terms
            for j in range(2, npar):
                A[i * npar + j, i * npar + j] += prior_m ** -2 * 0.25
        A += np.eye(nu) * 1e-9
        return np.linalg.solve(A, rhs)

    extra = np.ones(len(rows))
    for it in range(8):
        x = build(extra)
        res = np.array([math.hypot(sum(v * x[i] for i, v in eE.items()) - rE,
                                   sum(v * x[i] for i, v in eN.items()) - rN)
                        for eE, eN, rE, rN, w, _ in rows])
        sg = max(np.median(res) * 1.4826, 1.5)
        extra = 1.0 / (1.0 + (res / (2.5 * sg)) ** 2)
    cross = np.array([r[5] for r in rows])
    for nm, m in (('along', ~cross), ('cross', cross)):
        if m.any():
            log(f"           fit residual {nm:5s}: n={int(m.sum()):5d}  median {np.median(res[m]):5.1f}  "
                f"p90 {np.percentile(res[m],90):5.1f}  max {res[m].max():6.1f} m   downweighted {(extra[m]<0.25).sum()}")
    out = {}
    for r in recs:
        i = idx[r]
        d = dict(cE=float(x[i * npar]), cN=float(x[i * npar + 1]), sig=0.0, rho=0.0, s1=0.0, s2=0.0)
        if npar >= 4:
            d['sig'] = float(x[i * npar + 2] / L); d['rho'] = float(x[i * npar + 3] / L)
        if npar == 6:
            d['s1'] = float(x[i * npar + 4] / L); d['s2'] = float(x[i * npar + 5] / L)
        out[r] = d
    return out


def apply(sol, recs, corr):
    for r in recs:
        d = corr[r]
        sol[r]['dE'] += d['cE']; sol[r]['dN'] += d['cN']          # position = e - dE
        sol[r]['gw'] *= (1.0 + d['sig']); sol[r]['gh'] *= (1.0 + d['sig'])
        sol[r]['rot'] += math.degrees(d['rho'])
    return sol


def run(tag, rounds=3, start='stageA', model='similarity', feature='raw', out='stageA3'):
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
    recs = sorted([r for r in sol if sol[r].get('ok')])
    lab = SC.lines_of(sol, recs)
    print(f"\n===== close3 {tag}: from {start} -> {out}, model {model}, feature {feature}, {len(recs)} frames, "
          f"{max(lab.values())+1} lines =====", flush=True)
    if model == 'affine':
        print("  NOTE: affine shear terms are solved but not applied (render_frame is a "
              "similarity); use them only to read the residual", flush=True)
    mpp = C.MPP_TIE; imgdir = f"{SP}/fullres"
    for it in range(rounds):
        minE, maxN, W, H = C.canvas(sol, mpp)
        rend = blockadjust.render_all(sol, imgdir, minE, maxN, W, H, mpp, crop=0.95,
                                      log=lambda *_: None)
        ties = tie_windows(rend, lab, mpp, feature=feature)
        del rend
        st = seam_stats(ties, f"round {it}")
        if not ties:
            print("  no ties; stopping"); break
        corr = solve(sol, ties, recs, minE, maxN, mpp, model=model)
        if corr is None:
            break
        sg = np.array([abs(corr[r]['sig']) for r in recs]); rh = np.array([abs(corr[r]['rho']) for r in recs])
        mv = np.array([math.hypot(corr[r]['cE'], corr[r]['cN']) for r in recs])
        print(f"           per-frame |scale| p50 {np.median(sg)*1e2:.3f}% p90 {np.percentile(sg,90)*1e2:.3f}% max {sg.max()*1e2:.3f}%   "
              f"|rotation| p50 {math.degrees(np.median(rh)):.3f} p90 {math.degrees(np.percentile(rh,90)):.3f} max {math.degrees(rh.max()):.3f} deg", flush=True)
        print(f"           translations median {np.median(mv):.1f} max {mv.max():.1f} m", flush=True)
        dis = st.get('cross') or 0.0
        # Scale bound 3%: with the prior at 0.25 the median is ~0.4% and the tail
        # belongs to poorly-tied edge frames. Rotation is deliberately unbounded --
        # the bundle estimated each frame's crab from one central window swept in
        # 0.5 deg steps, where rotation is invisible, so ~1 deg errors are credible
        # and round 1's re-measurement is the test, not a threshold.
        if sg.max() > 0.03 or np.median(mv) > 2.0 * dis + 20:
            print(f"           REFUSED: implausible (scale {sg.max()*1e2:.2f}%, move {np.median(mv):.0f} m "
                  f"vs cross disagreement {dis:.0f} m)", flush=True); break
        sol = apply(sol, recs, corr)
        json.dump({r: dict(dE=sol[r]['dE'], dN=sol[r]['dN'], rot=sol[r]['rot'],
                           gw=sol[r]['gw'], gh=sol[r]['gh']) for r in recs},
                  open(P('data', f'{out}_{tag}.json'), 'w'))
        if np.median(mv) < 0.3 and sg.max() < 3e-5:
            break
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, imgdir, minE, maxN, W, H, mpp, crop=0.95, log=lambda *_: None)
    ties = tie_windows(rend, lab, mpp, feature=feature, log=lambda *_: None); del rend
    seam_stats(ties, "CLOSED")
    json.dump({r: dict(dE=sol[r]['dE'], dN=sol[r]['dN'], rot=sol[r]['rot'],
                       gw=sol[r]['gw'], gh=sol[r]['gh']) for r in recs},
              open(P('data', f'{out}_{tag}.json'), 'w'))
    print(f"  saved data/{out}_{tag}.json  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    run(sys.argv[1], int(arg('--rounds', 3)), arg('--from', 'stageA'), arg('--model', 'similarity'),
        arg('--feature', 'raw'), arg('--out', 'stageA3'))
