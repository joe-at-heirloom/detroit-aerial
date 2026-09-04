#!/usr/bin/env python
"""Look at the closed block's absolute error field, then fit the simplest model
that explains it -- and say what each model leaves behind.

Every failed placement so far fitted a model without checking it was the right
shape. This prints the raw field first, so the shape is visible before any fit,
then compares models by held-out residual, not in-sample residual.

Models, in order of freedom:
  shift        2 params   (what placeshift.py did: left 72 m)
  similarity   4 params   rotation + isotropic scale
  affine       6 params   adds shear and anisotropic scale
  quadratic   12 params   one gentle bend per axis
  per-line affine        one affine per flight line (frames only)

Usage:  ./.venv/bin/python scripts/fieldfit.py 1961 closedA [--use cells|frames]
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from validate import P


def design(E, N, kind, E0, N0, sc):
    x = (np.asarray(E) - E0) / sc; y = (np.asarray(N) - N0) / sc
    one = np.ones_like(x); z = np.zeros_like(x)
    if kind == 'shift':
        return np.stack([one, z], 1), np.stack([z, one], 1)
    if kind == 'similarity':
        # dE = tx + a x - b y ; dN = ty + b x + a y
        return np.stack([one, z, x, -y], 1), np.stack([z, one, y, x], 1)
    if kind == 'affine':
        return (np.stack([one, z, x, y, z, z], 1), np.stack([z, one, z, z, x, y], 1))
    if kind == 'quadratic':
        q = [one, x, y, x * x, x * y, y * y]
        zq = [z] * 6
        return np.stack(q + zq, 1), np.stack(zq + q, 1)
    if kind == 'strip':
        # A flight block is a strip: 11 km wide and 49 km long once 1949 roll ha-17
        # is solved whole. Its error bends along the strip and barely across it, so
        # spending parameters equally in both directions (the full 2D cubic) buys
        # nothing across and not enough along, and leave-one-out rejects it. This is
        # quartic in N, linear in E, with two cross terms: 8 per axis against the
        # cubic's 10, aimed where the bend actually is.
        q = [one, x, y, y * y, y ** 3, y ** 4, x * y, x * y * y]
        zq = [z] * 8
        return np.stack(q + zq, 1), np.stack(zq + q, 1)
    if kind == 'cubic':
        # A 33 km block (1956) leaves a cubic residual after the quadratic: the
        # along-track error is a parabola whose curvature itself changes along
        # the block. 20 parameters; only leave-one-out may choose it.
        q = [one, x, y, x * x, x * y, y * y, x ** 3, x * x * y, x * y * y, y ** 3]
        zq = [z] * 10
        return np.stack(q + zq, 1), np.stack(zq + q, 1)
    raise ValueError(kind)


def fit(E, N, dE, dN, w, kind, E0, N0, sc, iters=8):
    AE, AN = design(E, N, kind, E0, N0, sc)
    A = np.vstack([AE, AN]); b = np.concatenate([dE, dN])
    ww = np.concatenate([w, w]).astype(float)
    for _ in range(iters):
        sw = np.sqrt(ww)[:, None]
        p = np.linalg.lstsq(A * sw, b * np.sqrt(ww), rcond=None)[0]
        r = A @ p - b
        rr = np.hypot(r[:len(E)], r[len(E):])
        s = max(np.median(rr) * 1.4826, 2.0)
        rob = 1.0 / (1.0 + (rr / (2.5 * s)) ** 2)
        ww = np.concatenate([w * rob, w * rob])
    return p


def predict(p, E, N, kind, E0, N0, sc):
    AE, AN = design(E, N, kind, E0, N0, sc)
    return AE @ p, AN @ p


def lines_by_E(E, gap=800.0):
    """Flight-line index per observation, by clustering easting."""
    E = np.asarray(E); order = np.argsort(E); Es = E[order]
    lab = np.zeros(len(E), int); k = 0
    for i, oi in enumerate(order):
        if i > 0 and Es[i] - Es[i - 1] > gap:
            k += 1
        lab[oi] = k
    return lab


def fit_perline(E, N, dE, dN, w, sc):
    """One similarity per flight line. Physically motivated: each line is one
    pass at one height with one heading, so it has one scale and one rotation."""
    lab = lines_by_E(E)
    lines = []
    for k in range(lab.max() + 1):
        m = lab == k
        if m.sum() < 5:
            continue
        E0 = float(np.average(E[m], weights=w[m])); N0 = float(np.average(N[m], weights=w[m]))
        p = fit(E[m], N[m], dE[m], dN[m], w[m], 'similarity', E0, N0, sc)
        lines.append(dict(cE=E0, cN=N0, p=p.tolist(), n=int(m.sum())))
    return dict(kind='perline', lines=lines, sc=sc)


def predict_model(model, E, N):
    """Displacement (dE, dN) at (E, N) for any saved model kind."""
    E = np.atleast_1d(np.asarray(E, float)); N = np.atleast_1d(np.asarray(N, float))
    if model['kind'] != 'perline':
        return predict(np.array(model['p']), E, N, model['kind'], model['E0'], model['N0'], model['sc'])
    L = sorted(model['lines'], key=lambda l: l['cE'])
    sc = model['sc']
    cE = np.array([l['cE'] for l in L])
    # evaluate every line's similarity, then blend linearly in E between the two
    # nearest line centres -- continuous across the sidelap, exact at line centres
    preds = [predict(np.array(l['p']), E, N, 'similarity', l['cE'], l['cN'], sc) for l in L]
    # Each line's own similarity applies over the line's body; the transition to
    # the neighbouring line happens only across the sidelap, a band of width
    # `blend` centred on the midpoint between line centres. Blending over the whole
    # inter-line gap would smear each line's correction with its neighbour's across
    # most of the line -- and the correction is per line for a physical reason.
    blend = float(model.get('blend', 1200.0))
    pe = np.zeros_like(E); pn = np.zeros_like(N)
    for i in range(len(E)):
        x = E[i]
        j = int(np.argmin(np.abs(cE - x)))          # nearest line
        pe[i], pn[i] = preds[j][0][i], preds[j][1][i]
        # is x inside a sidelap band with a neighbour?
        for jn in (j - 1, j + 1):
            if 0 <= jn < len(L):
                mid = 0.5 * (cE[j] + cE[jn]); d = (x - mid) * (1 if jn > j else -1)
                # d in (-blend/2, +blend/2): 0 = midpoint, + toward neighbour
                if abs(x - mid) < blend / 2:
                    t = 0.5 + d / blend               # 0.5 at midpoint
                    pe[i] = (1 - t) * preds[j][0][i] + t * preds[jn][0][i]
                    pn[i] = (1 - t) * preds[j][1][i] + t * preds[jn][1][i]
    return pe, pn


def loo_perline(E, N, dE, dN, w, sc):
    out = []
    n = len(E)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        md = fit_perline(E[m], N[m], dE[m], dN[m], w[m], sc)
        pe, pn = predict_model(md, E[i:i+1], N[i:i+1])
        out.append(math.hypot(pe[0] - dE[i], pn[0] - dN[i]))
    return np.array(out)


def loo(E, N, dE, dN, w, kind, E0, N0, sc):
    """Leave-one-out residual: the honest number."""
    out = []
    n = len(E)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        if m.sum() < 8:
            return None
        p = fit(E[m], N[m], dE[m], dN[m], w[m], kind, E0, N0, sc)
        pe, pn = predict(p, E[i:i+1], N[i:i+1], kind, E0, N0, sc)
        out.append(math.hypot(pe[0] - dE[i], pn[0] - dN[i]))
    return np.array(out)


def main():
    tag, label = sys.argv[1], sys.argv[2]
    use = sys.argv[sys.argv.index('--use') + 1] if '--use' in sys.argv else 'cells'
    F = json.load(open(P('data', f'field_{tag}_{label}.json')))

    # --- show the raw field on the coarse grid ---
    cells = [c for c in F['cells_16x3'] if 'skip' not in c and not c.get('pegged')]
    print(f"\n{tag} {label}: raw 16x3 field (dE,dN in m; . = no lock)")
    rows = sorted(set(c['row'] for c in F['cells_16x3']))
    for r in rows:
        line = f"  row {r:2d} "
        for col in range(3):
            c = next((x for x in cells if x['row'] == r and x['col'] == col and x['ratio'] >= 1.15), None)
            line += f"  {c['dE']:+6.1f},{c['dN']:+6.1f}" if c else "        .      "
        print(line)

    if use == 'frames' and F.get('frames'):
        obs = [(f['e'], f['n'], f['dE'], f['dN'], f['ratio']) for f in F['frames']]
        src = f"{len(obs)} frames"
    else:
        g = [c for c in F['cells_32x6'] if 'skip' not in c and not c.get('pegged') and c['ratio'] >= 1.15]
        obs = [(c['cE'], c['cN'], c['dE'], c['dN'], c['ratio']) for c in g]
        src = f"{len(obs)} cells (32x6, ratio>=1.15)"
    E = np.array([o[0] for o in obs]); N = np.array([o[1] for o in obs])
    dE = np.array([o[2] for o in obs]); dN = np.array([o[3] for o in obs])
    w = np.array([min(max(o[4] - 1.0, 0.02) / 0.5, 3.0) for o in obs])
    E0, N0 = float(np.average(E, weights=w)), float(np.average(N, weights=w))
    sc = 10000.0
    mag = np.hypot(dE, dN)
    print(f"\nfitting to {src}: raw error median {np.median(mag):.1f}  p90 {np.percentile(mag,90):.1f} m")
    print(f"  {'model':12s} {'params':>6s}  {'in-sample med':>13s} {'p90':>6s}   {'leave-one-out med':>17s} {'p90':>6s} {'max':>6s}")
    best = None
    models = {}
    for kind in ('shift', 'similarity', 'affine', 'quadratic', 'cubic', 'strip'):
        p = fit(E, N, dE, dN, w, kind, E0, N0, sc)
        pe, pn = predict(p, E, N, kind, E0, N0, sc)
        r = np.hypot(pe - dE, pn - dN)
        lo = loo(E, N, dE, dN, w, kind, E0, N0, sc)
        lmed = float(np.median(lo)) if lo is not None else float('nan')
        print(f"  {kind:12s} {len(p):6d}  {np.median(r):13.1f} {np.percentile(r,90):6.1f}   "
              f"{lmed:17.1f} {np.percentile(lo,90) if lo is not None else float('nan'):6.1f} "
              f"{lo.max() if lo is not None else float('nan'):6.1f}")
        models[kind] = dict(kind=kind, p=p.tolist(), E0=E0, N0=N0, sc=sc,
                            loo_median=lmed, insample_median=float(np.median(r)))
        if kind == 'similarity':
            tx, ty, a, b = p          # x, y were scaled by sc, so slopes are a/sc, b/sc
            s_ = math.hypot(1 + a / sc, b / sc); th = math.degrees(math.atan2(b / sc, 1 + a / sc))
            print(f"               -> shift {tx:+.1f},{ty:+.1f} m   scale {s_:.5f} ({(s_-1)*1e2:+.3f}%)"
                  f"   rotation {th:+.4f} deg")
        if kind == 'affine':
            tx, ty, axx, axy, ayx, ayy = p
            print(f"               -> shift {tx:+.1f},{ty:+.1f} m   d(dE)/dE {axx/sc*1e2:+.3f}%  d(dE)/dN {axy/sc*1e2:+.3f}%"
                  f"   d(dN)/dE {ayx/sc*1e2:+.3f}%  d(dN)/dN {ayy/sc*1e2:+.3f}%")
        if lo is not None and (best is None or lmed < best[1]):
            best = (kind, lmed)

    # per-line similarity: needs enough observations per line; prefer frames
    src_pl = None
    if F.get('frames') and len(F['frames']) >= 20:
        fr = F['frames']
        Ep = np.array([f['e'] for f in fr]); Np = np.array([f['n'] for f in fr])
        dEp = np.array([f['dE'] for f in fr]); dNp = np.array([f['dN'] for f in fr])
        wp = np.array([min(max(f['ratio'] - 1.0, 0.02) / 0.5, 3.0) for f in fr]); src_pl = 'frames'
    else:
        Ep, Np, dEp, dNp, wp = E, N, dE, dN, w; src_pl = 'cells'
    md = fit_perline(Ep, Np, dEp, dNp, wp, sc)
    if md['lines']:
        pe, pn = predict_model(md, Ep, Np); r = np.hypot(pe - dEp, pn - dNp)
        lo = loo_perline(Ep, Np, dEp, dNp, wp, sc)
        npar = 4 * len(md['lines'])
        print(f"  {'perline':12s} {npar:6d}  {np.median(r):13.1f} {np.percentile(r,90):6.1f}   "
              f"{np.median(lo):17.1f} {np.percentile(lo,90):6.1f} {lo.max():6.1f}   (on {src_pl})")
        for l in md['lines']:
            tx, ty, a, b = l['p']
            s_ = math.hypot(1 + a / sc, b / sc); th = math.degrees(math.atan2(b / sc, 1 + a / sc))
            print(f"               line E~{l['cE']:.0f}: n={l['n']}  shift {tx:+.1f},{ty:+.1f}  "
                  f"scale {(s_-1)*1e2:+.3f}%  rot {th:+.4f} deg")
        models['perline'] = dict(md, loo_median=float(np.median(lo)), insample_median=float(np.median(r)))
        if best is None or float(np.median(lo)) < best[1]:
            best = ('perline', float(np.median(lo)))

    # per-flight-line structure, if frames are available
    if F.get('frames'):
        fr = F['frames']
        Ef = np.array([f['e'] for f in fr]); Nf = np.array([f['n'] for f in fr])
        order = np.argsort(Ef); Es = Ef[order]
        lines = []; cur = [order[0]]
        for i in range(1, len(order)):
            if Es[i] - Es[i-1] > 800: lines.append(cur); cur = []
            cur.append(order[i])
        lines.append(cur)
        print(f"\nper-flight-line view ({len(lines)} lines, from {len(fr)} matched frames):")
        for k, idx in enumerate(lines):
            idx = np.array(idx)
            if len(idx) < 4: continue
            n_ = Nf[idx]; de = np.array([fr[i]['dE'] for i in idx]); dn = np.array([fr[i]['dN'] for i in idx])
            sE = np.polyfit(n_, de, 1); sN = np.polyfit(n_, dn, 1)
            print(f"  line {k} (E~{Ef[idx].mean():.0f}): {len(idx)} frames, mean err ({de.mean():+.1f},{dn.mean():+.1f}),"
                  f"  along-track slope dE {sE[0]*1e2:+.3f}%  dN {sN[0]*1e2:+.3f}%  (N-scale error)")

    choice = best[0] if best else 'similarity'
    print(f"\nbest by leave-one-out: {choice} ({best[1]:.1f} m)" if best else "")
    json.dump(dict(models=models, choice=choice, source=src),
              open(P('data', f'fieldfit_{tag}_{label}.json'), 'w'))
    print(f"saved data/fieldfit_{tag}_{label}.json")


if __name__ == '__main__':
    main()
