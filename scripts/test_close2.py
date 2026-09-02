#!/usr/bin/env python
"""Sign and yield tests for the sidelap matcher and the block solver.

The expected offset is DERIVED from how the test moves the content, not typed by
hand -- three times running the matcher was right and the hand-written
expectation was wrong about north. The convention throughout the project:
the reported (dE, dN) is where A's content sits relative to B's, east and north
positive. Moving B's content by +px columns and +py rows puts it px*mpp east and
py*mpp SOUTH of A, so A sits (-px*mpp, +py*mpp) from B.

Run:  ./.venv/bin/python scripts/test_close2.py
"""
import sys, os, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import close2
from scipy.ndimage import gaussian_filter


def expected(px, py, mpp):
    return (-px * mpp, +py * mpp)


def synthetic_city(seed=1, n=1400, period=49):
    rng = np.random.default_rng(seed)
    big = gaussian_filter(rng.normal(0, 1, (n, n)).astype(np.float32), 3) * 40 + 128
    for y in range(0, n, period): big[y:y + 3, :] += 60
    for x in range(0, n, period): big[:, x:x + 3] += 60
    return big


def test_sidelap_matcher():
    mpp = 2.0; big = synthetic_city()
    ok = True
    for px, py in ((90, 30), (105, -40), (20, 10), (-70, 55)):
        A = np.clip(big, 1, 255).astype(np.uint8)
        B = np.clip(np.roll(np.roll(big, py, 0), px, 1), 1, 255).astype(np.uint8)
        A[:, :300] = 0; B[:, 900:] = 0
        obs = close2.cross_observations({'a': (A, 0, 0), 'b': (B, 0, 0)}, {'a': {}, 'b': {}},
                                        {'a': 0, 'b': 1}, mpp, log=lambda *_: None)
        eE, eN = expected(px, py, mpp)
        got = obs[0] if obs else None
        good = got is not None and abs(got['dE'] - eE) < 3 and abs(got['dN'] - eN) < 3 and got['nwin'] >= 3
        print(f"  sidelap planted ({eE:+.0f},{eN:+.0f}) m -> "
              f"{'(%+.0f,%+.0f) %d/%d windows' % (got['dE'], got['dN'], got['nwin'], got['nwin_all']) if got else 'nothing'}"
              f"  {'OK' if good else 'FAIL'}")
        ok &= good
    return ok


def test_block_solver():
    """Two lines, line 1 carries +0.5% scale: the solver must return the correction
    that cancels it, and applying it must leave zero disagreement."""
    sol = {}; recs = []; lab = {}
    for k, E in enumerate((-2300., 0.)):
        for j in range(6):
            r = f"{k}_{j}"; sol[r] = dict(ok=True, e=E, n=j * 1300., dE=0., dN=0., rot=0., gw=3400., gh=3200.)
            recs.append(r); lab[r] = k
    ctr = {k: np.mean([[sol[r]['e'], sol[r]['n']] for r in recs if lab[r] == k], axis=0) for k in (0, 1)}
    true_sig = {0: 0.0, 1: 0.005}
    mpp = 2.0; minE = -5000.; maxN = 10000.
    def shift(r, q):
        v = q - ctr[lab[r]]; return true_sig[lab[r]] * v
    obs = []
    for i, a in enumerate(recs):
        for b in recs[i + 1:]:
            pa = np.array([sol[a]['e'], sol[a]['n']]); pb = np.array([sol[b]['e'], sol[b]['n']])
            if np.linalg.norm(pa - pb) > 2700: continue
            q = (pa + pb) / 2; d = shift(a, q) - shift(b, q)
            obs.append(dict(a=a, b=b, dE=d[0], dN=d[1], ratio=2.0, cross=lab[a] != lab[b],
                            qE=(q[0] - minE) / mpp, qN=(maxN - q[1]) / mpp))
    cE, cN, sig, rho, c = close2.solve_block(sol, obs, recs, lab, minE, maxN, mpp, log=lambda *_: None)
    diff = (sig[1] - sig[0])
    good = abs(diff + 0.005) < 1e-4       # correction cancels the +0.5% error
    print(f"  solver: line1-line0 scale correction {diff*1e2:+.3f}% (expect -0.500)  {'OK' if good else 'FAIL'}")
    return good


if __name__ == '__main__':
    a = test_sidelap_matcher(); b = test_block_solver()
    print("ALL OK" if a and b else "FAILURES")
    sys.exit(0 if a and b else 1)
