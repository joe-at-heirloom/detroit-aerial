#!/usr/bin/env python
"""Solve the film years' corrections together, on measurements alone (plan step 5).

Each film block was placed on its own -- 1961 against modern imagery, 1949 direct
to modern, 1956 chained through 1961 -- so each carries its own residual, and a
wipe between two of them shows both. The pair matrix (crossmatrix.py) measures
every film year against the reference AND against every other film year. Here a
correction field per film year is fitted to all of those at once:

    film i vs reference   m = d_i(x)
    film i vs film j      m = d_i(x) - d_j(x)

where d_i is how far layer i's content sits from the ground; the correction a
raster would take is -d_i. The reference (NAIP 2016) is held fixed. Esri hi is
never fitted and is the independent check. Other modern years are not used.

Model order is earned, not chosen: shift, similarity, affine, quadratic and cubic
are each scored on SPATIALLY withheld cells (4 km blocks, 5 folds), on both kinds
of pair, and against Esri hi. Cells over --gate metres are treated as unverifiable
(changed ground: 1956's north-west locks 100 m off on suburbs that were farmland)
and listed, not fitted.

Nothing here moves a raster.

Usage:  ./.venv/bin/python scripts/jointfit.py runs/crossmatrix_west_adj.json
            [--ref l2016] [--check esrihi] [--gate 30] [--block 4000] [--out runs/jointfit_west.json]
"""
import sys, os, json, math, argparse
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import dtmap

FILM = ['b1949', 'b1956', 'b1961', 'b1967']
ORDERS = ['none', 'shift', 'similarity', 'affine', 'quadratic', 'cubic']
SC = 10000.0          # coordinate scale, metres


def basis(order, x, y):
    """Rows for dE and dN of one layer's field at (x, y), normalised coords.
    Returns (2, k)."""
    if order == 'none':
        return np.zeros((2, 0))
    if order == 'shift':
        return np.array([[1, 0], [0, 1]], float)
    if order == 'similarity':
        return np.array([[1, 0, x, -y], [0, 1, y, x]], float)
    deg = {'affine': 1, 'quadratic': 2, 'cubic': 3}[order]
    t = [x ** i * y ** j for i in range(deg + 1) for j in range(deg + 1 - i)]
    k = len(t); B = np.zeros((2, 2 * k))
    B[0, :k] = t; B[1, k:] = t
    return B


def cells(run, ref, check, min_ratio, gate, want_unmeasured=False):
    """Usable observations, the ones gated out as unverifiable, and (optionally)
    where a film layer has imagery against the reference but nothing measured."""
    obs, gated, unmeasured = [], [], []
    for p in run['pairs']:
        a, b = p['a'], p['b']
        kinds = {a, b}
        if not (kinds <= set(FILM) | {ref, check}) or not (kinds & set(FILM)):
            continue
        if kinds & {ref, check} == kinds:
            continue
        for c in p.get('cells', []):
            E = (c['lon'] - dtmap.LON0) * dtmap.MLON; N = (c['lat'] - dtmap.LAT0) * dtmap.MLAT
            bad = ('skip' in c and c['skip'] != 'coverage') or ('skip' not in c and (
                c['pegged'] or c['ratio'] < min_ratio or c.get('track', 0.0) > 5.0))
            if 'skip' in c or bad:
                film_side = 'va' if a in FILM else 'vb'
                if bad or (c['skip'] == 'coverage' and c.get(film_side, 0) >= 0.45):
                    if ref in kinds:
                        unmeasured.append(dict(a=a, b=b, E=E, N=N, dE=0.0, dN=0.0, lat=c['lat'], lon=c['lon']))
                continue
            o = dict(a=a, b=b, E=E, N=N, dE=c['dE'], dN=c['dN'], lat=c['lat'], lon=c['lon'])
            if math.hypot(c['dE'], c['dN']) > gate:
                gated.append(o)
                if ref in kinds:
                    unmeasured.append(dict(o, dE=0.0, dN=0.0))
            else:
                obs.append(o)
    return (obs, gated, unmeasured) if want_unmeasured else (obs, gated)


class Model:
    def __init__(self, order, obs):
        self.order = order if isinstance(order, dict) else {f: order for f in FILM}
        cen = {}
        for f in FILM:
            P = [(o['E'], o['N']) for o in obs if f in (o['a'], o['b'])]
            cen[f] = np.mean(P, axis=0) if P else np.zeros(2)
        self.cen = cen
        self.kf = {f: basis(self.order[f], 0.0, 0.0).shape[1] for f in FILM}
        self.col, c = {}, 0
        for f in FILM:
            self.col[f] = c; c += self.kf[f]
        self.k = c

    def B(self, f, E, N):
        return basis(self.order[f], (E - self.cen[f][0]) / SC, (N - self.cen[f][1]) / SC)

    def rows(self, o):
        R = np.zeros((2, self.k))
        for side, s in (('a', 1.0), ('b', -1.0)):
            f = o[side]
            if f in self.col and self.kf[f]:
                R[:, self.col[f]:self.col[f] + self.kf[f]] += s * self.B(f, o['E'], o['N'])
        return R

    def fit(self, obs, ridge=1e-3):
        self.p = np.zeros(self.k)
        if not self.k or not obs:
            return self
        X = np.vstack([self.rows(o) for o in obs]); y = np.concatenate([[o['dE'], o['dN']] for o in obs])
        # a prior observation ('no change here') has weight w0 and is exempt from the
        # robust down-weighting that would otherwise discard it exactly where it matters
        w0 = np.repeat([o.get('w', 1.0) for o in obs], 2)
        prior = np.repeat([o.get('prior', False) for o in obs], 2)
        w = w0.copy()
        for _ in range(8):
            A = X * w[:, None]
            self.p = np.linalg.solve(A.T @ A + ridge * np.eye(X.shape[1]), A.T @ (y * w))
            r = y - X @ self.p; s = max(np.median(np.abs(r[~prior])) * 1.4826, 0.5)
            w = w0 / np.sqrt(1.0 + (r / (2 * s)) ** 2)           # Cauchy
            w[prior] = w0[prior]
        return self

    def d(self, f, E, N):
        if f not in self.col or not self.kf[f]:
            return np.zeros(2)
        return self.B(f, E, N) @ self.p[self.col[f]:self.col[f] + self.kf[f]]

    def resid(self, o):
        pred = self.d(o['a'], o['E'], o['N']) - self.d(o['b'], o['E'], o['N'])
        return math.hypot(o['dE'] - pred[0], o['dN'] - pred[1])


def stats(v):
    v = np.asarray(v)
    return dict(n=int(len(v)), median=float(np.median(v)) if len(v) else None,
                p90=float(np.percentile(v, 90)) if len(v) else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('--ref', default='l2016'); ap.add_argument('--check', default='esrihi')
    ap.add_argument('--min-ratio', type=float, default=1.15)
    ap.add_argument('--gate', type=float, default=30.0)
    ap.add_argument('--block', type=float, default=4000.0)
    ap.add_argument('--folds', type=int, default=5)
    ap.add_argument('--out', default=os.path.join(ROOT, 'runs', 'jointfit_west.json'))
    ap.add_argument('--damp', type=float, default=0.0,
                    help="weight of 'no change' pseudo-observations where a film has imagery but nothing measured")
    ap.add_argument('--force', default='', help='fix orders, e.g. b1956=quadratic,b1949=none')
    a = ap.parse_args()
    run = json.load(open(a.run))
    allobs, gated, unmeas = cells(run, a.ref, a.check, a.min_ratio, a.gate, want_unmeasured=True)
    fitobs = [o for o in allobs if a.check not in (o['a'], o['b'])]
    priors = [dict(o, w=a.damp, prior=True) for o in unmeas] if a.damp > 0 else []
    forced = dict(kv.split('=') for kv in a.force.split(',') if kv)
    if priors:
        print(f"prior: {len(priors)} 'no change' pseudo-observations at weight {a.damp} where film has "
              f"imagery but nothing was measured (fitting only, never scored)")
    chk = [o for o in allobs if a.check in (o['a'], o['b'])]
    kind = lambda o: 'ref' if a.ref in (o['a'], o['b']) else 'film'
    print(f"observations: {sum(kind(o)=='ref' for o in fitobs)} film-vs-{a.ref}, "
          f"{sum(kind(o)=='film' for o in fitobs)} film-vs-film, {len(chk)} film-vs-{a.check} (never fitted); "
          f"{len(gated)} over {a.gate:.0f} m gated as unverifiable", flush=True)
    for f in FILM:
        g = [o for o in gated if f in (o['a'], o['b'])]
        if g:
            print(f"    gated {f}: {len(g)} cells, e.g. {g[0]['lat']:.3f},{g[0]['lon']:.3f}")

    # spatial folds: 4 km blocks, shuffled
    rng = np.random.default_rng(0)
    blk = {}
    for o in fitobs:
        key = (int(o['E'] // a.block), int(o['N'] // a.block))
        blk.setdefault(key, rng.integers(0, a.folds))
    fold = np.array([blk[(int(o['E'] // a.block), int(o['N'] // a.block))] for o in fitobs])

    report = dict(run=a.run, ref=a.ref, check=a.check, gate=a.gate, block=a.block, damp=a.damp,
                  forced=forced, models={})
    print(f"\n  held-out residual by 4 km block, median (p90) m; '{a.check}' uses the full fit and is never fitted")
    hdr = ''.join(f"{f[1:]:>14s}" for f in FILM)
    print(f"  {'model':12s} {'film/ref':>12s} {'film/film':>12s} {'vs '+a.check:>12s} |{hdr}")
    for order in ORDERS:
        held = {'ref': [], 'film': []}; per = {f: [] for f in FILM}
        for k in range(a.folds):
            tr = [o for o, fo in zip(fitobs, fold) if fo != k]
            te = [o for o, fo in zip(fitobs, fold) if fo == k]
            m = Model(order, tr).fit(tr + priors)
            for o in te:
                r = m.resid(o); held[kind(o)].append(r)
                if kind(o) == 'ref':
                    per[o['a'] if o['a'] in FILM else o['b']].append(r)
        full = Model(order, fitobs).fit(fitobs + priors)
        ck = [full.resid(o) for o in chk]
        rec = dict(ref=stats(held['ref']), film=stats(held['film']), check=stats(ck),
                   per_layer={f: stats(v) for f, v in per.items()})
        report['models'][order] = rec
        fmt = lambda s: f"{s['median']:5.2f} ({s['p90']:4.1f})" if s['n'] else '      -     '
        print(f"  {order:12s} {fmt(rec['ref']):>12s} {fmt(rec['film']):>12s} {fmt(rec['check']):>12s} |"
              + ''.join(f"{fmt(rec['per_layer'][f]):>14s}" for f in FILM), flush=True)
    # Per-layer orders. A narrow block (1949) cannot carry a cubic -- withheld 4 km
    # blocks are extrapolated wildly -- while 1956 needs one. Coordinate descent over
    # each layer's order, scored on the same withheld blocks; Esri hi stays out of
    # the choice and is reported only for the winner.
    def cv(order):
        held = {'ref': [], 'film': []}; per = {f: [] for f in FILM}
        for k in range(a.folds):
            tr = [o for o, fo in zip(fitobs, fold) if fo != k]
            te = [o for o, fo in zip(fitobs, fold) if fo == k]
            m = Model(order, tr).fit(tr + priors)
            for o in te:
                r = m.resid(o); held[kind(o)].append(r)
                if kind(o) == 'ref':
                    per[o['a'] if o['a'] in FILM else o['b']].append(r)
        return held, per
    score = lambda h: float(np.median(h['ref']) + np.median(h['film']))
    best = {f: forced.get(f, 'none') for f in FILM}
    cur = score(cv(best)[0])
    for _ in range(3):
        changed = False
        for f in FILM:
            if f in forced:
                continue
            for order in ORDERS:
                trial = dict(best, **{f: order})
                sc_ = score(cv(trial)[0])
                if sc_ < cur - 1e-6:
                    cur, best, changed = sc_, trial, True
        if not changed:
            break
    held, per = cv(best)
    full = Model(best, fitobs).fit(fitobs + priors)
    ck = [full.resid(o) for o in chk]
    ckper = {f: stats([full.resid(o) for o in chk if f in (o['a'], o['b'])]) for f in FILM}
    ck0 = {f: stats([math.hypot(o['dE'], o['dN']) for o in chk if f in (o['a'], o['b'])]) for f in FILM}
    report['best'] = best
    report['best_cv'] = dict(ref=stats(held['ref']), film=stats(held['film']), check=stats(ck),
                             per_layer={f: stats(v) for f, v in per.items()},
                             check_per_layer=ckper, check_per_layer_uncorrected=ck0)
    fmt = lambda s: f"{s['median']:5.2f} ({s['p90']:4.1f})" if s['n'] else '      -     '
    print(f"\n  per-layer orders (coordinate descent on held-out blocks): {best}")
    print(f"  {'':12s} {fmt(stats(held['ref'])):>12s} {fmt(stats(held['film'])):>12s} {fmt(stats(ck)):>12s} |"
          + ''.join(f"{fmt(stats(per[f])):>14s}" for f in FILM))
    print("  against Esri hi, never fitted, per layer: uncorrected -> corrected")
    for f in FILM:
        print(f"    {f}: {fmt(ck0[f])} -> {fmt(ckper[f])}")
    # the separate (not joint) solve, for comparison: film-vs-film observations dropped
    sep = [o for o in fitobs if kind(o) == 'ref']
    for label, obs in (('joint', fitobs), ('separate', sep)):
        m = Model(best, obs).fit(obs + priors)
        ff = [m.resid(o) for o in fitobs if kind(o) == 'film']
        report[f'{label}_film_film_in_sample'] = stats(ff)
    print(f"  film/film in-sample with {best}: joint {report['joint_film_film_in_sample']['median']:.2f} m, "
          f"separate (reference only) {report['separate_film_film_in_sample']['median']:.2f} m")
    m = Model(best, fitobs).fit(fitobs + priors)
    report['params'] = {f: m.p[m.col[f]:m.col[f] + m.kf[f]].tolist() for f in FILM}
    report['centres'] = {f: m.cen[f].tolist() for f in FILM}
    report['gated'] = gated
    json.dump(report, open(a.out, 'w'), indent=1)
    print(f"\nwrote {a.out}")


if __name__ == '__main__':
    main()
