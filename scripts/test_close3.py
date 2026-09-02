#!/usr/bin/env python
"""close3 solver test: per-frame scale and rotation planted on a synthetic two-line
block, observed through dense tie windows, must be recovered and must leave zero
disagreement once applied. Expected signs are derived, not typed."""
import sys, os, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
import close3

def main():
    mpp = 2.0; minE = -6000.; maxN = 10000.
    sol = {}; recs = []; lab = {}
    for k, E in enumerate((-2300., 0.)):
        for j in range(6):
            r = f"{k}_{j}"; sol[r] = dict(ok=True, e=E, n=j * 1300., dE=0., dN=0., rot=0., gw=3400., gh=3200.)
            recs.append(r); lab[r] = k
    # truth: line 1 frames carry +0.5% scale, frame 0_3 carries +0.3 deg rotation,
    # frame 1_4 an extra +0.2% -- per FRAME, not per line
    true = {r: dict(sig=(0.005 if lab[r] == 1 else 0.0) + (0.002 if r == '1_4' else 0.0),
                    rho=(math.radians(0.3) if r == '0_3' else 0.0)) for r in recs}
    pos = {r: np.array([sol[r]['e'], sol[r]['n']]) for r in recs}
    def content(r, q):
        v = q - pos[r]; t = true[r]
        return t['sig'] * v + t['rho'] * np.array([-v[1], v[0]])
    ties = []
    for i, a in enumerate(recs):
        for b in recs[i + 1:]:
            if np.linalg.norm(pos[a] - pos[b]) > 2700: continue
            # dense windows across the overlap of the two footprints
            lo = np.maximum(pos[a] - [1700, 1600], pos[b] - [1700, 1600]); hi = np.minimum(pos[a] + [1700, 1600], pos[b] + [1700, 1600])
            for qx in np.arange(lo[0] + 192, hi[0] - 192, 192):
                for qy in np.arange(lo[1] + 192, hi[1] - 192, 192):
                    q = np.array([qx, qy]); d = content(a, q) - content(b, q)
                    ties.append(dict(a=a, b=b, qx=(q[0] - minE) / mpp, qy=(maxN - q[1]) / mpp,
                                     dE=d[0], dN=d[1], sharp=20.0, cross=lab[a] != lab[b]))
    corr = close3.solve(sol, ties, recs, minE, maxN, mpp, model='similarity', log=lambda *_: None)
    ms = np.mean([true[r]['sig'] for r in recs]); mr = np.mean([true[r]['rho'] for r in recs])
    ok = True
    for r in ('0_0', '1_0', '1_4', '0_3'):
        es = -(true[r]['sig'] - ms); er = -(true[r]['rho'] - mr)     # correction cancels the error, mean pinned
        gs, gr = corr[r]['sig'], corr[r]['rho']
        good = abs(gs - es) < 2e-4 and abs(gr - er) < 2e-4
        print(f"  {r}: scale corr {gs*1e2:+.3f}% (expect {es*1e2:+.3f})   rot corr {math.degrees(gr):+.3f} deg (expect {math.degrees(er):+.3f})  {'OK' if good else 'FAIL'}")
        ok &= good
    # apply and re-evaluate disagreement with the corrected geometry
    sol2 = close3.apply({r: dict(sol[r]) for r in recs}, recs, corr)
    res = []
    for t in ties:
        q = np.array([minE + t['qx'] * mpp, maxN - t['qy'] * mpp])
        def now(r):
            p0 = pos[r]; p = p0 - np.array([sol2[r]['dE'], sol2[r]['dN']])
            s = sol2[r]['gw'] / 3400. - 1.; th = math.radians(sol2[r]['rot']); v = q - p0
            return content(r, q) + (p - p0) + s * v + th * np.array([-v[1], v[0]])
        res.append(np.linalg.norm(now(t['a']) - now(t['b'])))
    print(f"  residual after applying: median {np.median(res):.3f}  max {np.max(res):.3f} m   "
          f"(before: median {np.median([math.hypot(t['dE'], t['dN']) for t in ties]):.1f})")
    ok &= np.max(res) < 0.05
    print("ALL OK" if ok else "FAILURES"); sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
