"""Locally-weighted (RBF) warp.

A global polynomial is fine in the middle of a block and diverges past the last
control station, which is why the 1967 end chunks were 80-106 m. A Gaussian-kernel
local regression instead decays toward the nearest stations' consensus, so
extrapolation degrades gracefully rather than exploding."""
import numpy as np, json, math


def _prep(stations, min_ratio=1.10):
    st = [s for s in stations if s.get('ratio', 0) >= min_ratio]
    P = np.array([[s['cE'], s['cN']] for s in st], float)
    V = np.array([[s['dE'], s['dN'], s.get('rot', 0.0)] for s in st], float)
    W = np.array([min(s.get('ratio', 1.0), 4.0) for s in st], float)
    return P, V, W


def robust_trim(P, V, W, length=3500.0, k=3.0):
    """Drop stations that disagree with their own neighbourhood."""
    keep = np.ones(len(P), bool)
    for i in range(len(P)):
        d2 = ((P - P[i]) ** 2).sum(1)
        w = W * np.exp(-d2 / (2 * length ** 2))
        w[i] = 0.0
        if w.sum() < 1e-6:
            continue
        loc = (V[:, :2] * w[:, None]).sum(0) / w.sum()
        resid = math.hypot(*(V[i, :2] - loc))
        spread = math.sqrt(max((w * ((V[:, :2] - loc) ** 2).sum(1)).sum() / w.sum(), 1.0))
        if resid > k * spread:
            keep[i] = False
    return keep


class RBFWarp:
    """Gaussian-kernel local-linear regression with ridge damping."""

    def __init__(self, stations, length=4000.0, ridge=1e-3, min_ratio=1.10,
                 trim=True, local_linear=True):
        P, V, W = _prep(stations, min_ratio)
        if trim and len(P) > 12:
            keep = robust_trim(P, V, W, length)
            P, V, W = P[keep], V[keep], W[keep]
        self.P, self.V, self.W = P, V, W
        self.length = float(length)
        self.ridge = float(ridge)
        self.local_linear = bool(local_linear)
        self.n = len(P)

    def at(self, cE, cN):
        p = np.array([cE, cN], float)
        d2 = ((self.P - p) ** 2).sum(1)
        w = self.W * np.exp(-d2 / (2 * self.length ** 2))
        s = w.sum()
        if s < 1e-9:
            j = int(np.argmin(d2))
            return tuple(self.V[j])
        # effective sample size decides whether a local plane is supportable
        ess = s ** 2 / max((w ** 2).sum(), 1e-12)
        if self.local_linear and ess >= 6.0:
            X = np.column_stack([np.ones(len(self.P)),
                                 (self.P[:, 0] - p[0]) / self.length,
                                 (self.P[:, 1] - p[1]) / self.length])
            sw = np.sqrt(w)[:, None]
            A = X * sw
            G = A.T @ A + self.ridge * np.eye(3) * max(s, 1.0)
            out = []
            for k in range(3):
                b = A.T @ (self.V[:, k] * np.sqrt(w))
                out.append(float(np.linalg.solve(G, b)[0]))
            return tuple(out)
        return tuple((self.V * w[:, None]).sum(0) / s)

    def dump(self, path):
        json.dump(dict(kind='rbf', length=self.length, ridge=self.ridge,
                       local_linear=self.local_linear,
                       P=self.P.tolist(), V=self.V.tolist(), W=self.W.tolist()),
                  open(path, 'w'))

    @staticmethod
    def load(path):
        d = json.load(open(path))
        o = RBFWarp.__new__(RBFWarp)
        o.P = np.array(d['P']); o.V = np.array(d['V']); o.W = np.array(d['W'])
        o.length = d['length']; o.ridge = d['ridge']
        o.local_linear = d.get('local_linear', True); o.n = len(o.P)
        return o


# --- trend + RBF as one fitted object ----------------------------------------
# The trend is fitted first and the RBF models what the trend leaves behind, so
# outside station coverage the RBF term decays to zero and the trend carries.

def trend_fit(st, deg=1):
    import numpy as _np
    P = _np.array([[s['cE'], s['cN']] for s in st], float)
    V = _np.array([[s['dE'], s['dN']] for s in st], float)
    W = _np.array([min(s.get('ratio', 1.0), 4.0) for s in st], float)
    mx, my = P[:, 0].mean(), P[:, 1].mean()
    sx, sy = max(P[:, 0].std(), 1.0), max(P[:, 1].std(), 1.0)
    import warpapply
    A = warpapply.design((P[:, 0] - mx) / sx, (P[:, 1] - my) / sy, deg)
    sw = _np.sqrt(W)[:, None]
    C = [_np.linalg.lstsq(A * sw, V[:, k] * _np.sqrt(W), rcond=None)[0] for k in range(2)]
    return dict(deg=deg, mx=float(mx), my=float(my), sx=float(sx), sy=float(sy),
                C=[c.tolist() for c in C])


def trend_at(t, cE, cN):
    import numpy as _np, warpapply
    A = warpapply.design(_np.array([(cE - t['mx']) / t['sx']]),
                         _np.array([(cN - t['my']) / t['sy']]), t['deg'])
    return (float((A @ _np.array(t['C'][0]))[0]), float((A @ _np.array(t['C'][1]))[0]))


def fit(st, length=800.0, deg=1, **kw):
    """Fit trend then RBF-on-residual. Returns (trend, RBFWarp, warp_at)."""
    t = trend_fit(st, deg)
    res = []
    for s in st:
        a = trend_at(t, s['cE'], s['cN'])
        res.append(dict(cE=s['cE'], cN=s['cN'], dE=s['dE'] - a[0], dN=s['dN'] - a[1],
                        rot=0.0, ratio=s.get('ratio', 1.0)))
    R = RBFWarp(res, length=length, min_ratio=0.0, **kw)
    def warp_at(cE, cN, t=t, R=R):
        a = trend_at(t, cE, cN); d = R.at(cE, cN)
        return (a[0] + d[0], a[1] + d[1], 0.0)
    return t, R, warp_at


def _spatial_folds(st, folds, seed):
    """Fold assignment by location, not at random.

    Control windows overlap by design, so a randomly held-out station almost always
    has a near-duplicate left in the training set. Random-fold cross-validation
    therefore rewards ever-shorter length scales -- it is scoring interpolation
    between neighbours, not prediction. Holding out contiguous patches instead
    makes the held-out error mean what it is supposed to mean."""
    import numpy as _np
    P = _np.array([[s['cE'], s['cN']] for s in st], float)
    cell = 2500.0
    gx = _np.floor((P[:, 0] - P[:, 0].min()) / cell).astype(int)
    gy = _np.floor((P[:, 1] - P[:, 1].min()) / cell).astype(int)
    keys = gx + 1000 * gy
    uniq = _np.unique(keys)
    rng = _np.random.default_rng(seed)
    assign = {k: int(v) for k, v in zip(uniq, rng.integers(0, folds, len(uniq)))}
    return _np.array([assign[k] for k in keys])


def cv_length(st, lengths=(250., 400., 600., 800., 1200., 1800., 2400.),
              folds=6, seeds=3, deg=1, spatial=True):
    """Pick the RBF length scale by held-out station error.

    Reported as median held-out residual, not mean: one aliased station left in the
    set would otherwise choose the length scale for all the others."""
    import numpy as _np, math as _m
    best = None; table = []
    for L in lengths:
        errs = []
        for seed in range(seeds):
            rng = _np.random.default_rng(seed)
            idx = rng.permutation(len(st))
            fold_of = _spatial_folds(st, folds, seed) if spatial else None
            for f in range(folds):
                te = (set(_np.nonzero(fold_of == f)[0].tolist()) if spatial
                      else set(idx[f::folds].tolist()))
                tr = [s for k, s in enumerate(st) if k not in te]
                ho = [s for k, s in enumerate(st) if k in te]
                if len(tr) < 12 or not ho:
                    continue
                _, _, w = fit(tr, length=L, deg=deg)
                for s in ho:
                    p = w(s['cE'], s['cN'])
                    errs.append(_m.hypot(s['dE'] - p[0], s['dN'] - p[1]))
        if not errs:
            continue
        med = float(_np.median(errs)); p90 = float(_np.percentile(errs, 90))
        table.append((L, med, p90, len(errs)))
        if best is None or med < best[1]:
            best = (L, med, p90)
    return best, table
