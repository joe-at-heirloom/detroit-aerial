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
