#!/usr/bin/env python
"""Why does multi-window phase correlation accept so few sidelap pairs?

On 1961 it accepted 21 of 224 sidelaps, each with the minimum three agreeing
windows. This looks at every window of a sample of sidelaps and reports what
killed it -- coverage, flat texture, weak peak, boundary -- and tries the feature
and window-size variants that might rescue it: raw film versus the road-ridge
response (lines flown hours apart have different shadows; roads do not move), and
800 m versus 512 m windows (a 1.1 km sidelap fits only one column of 800 m ones).

Usage:  ./.venv/bin/python scripts/pcdiag.py 1961 stageA [--n 24]
"""
import sys, os, json, math, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust
from validate import P
from rebuild import placements, SP
import close as C, seamclass as SC, close2


def main():
    tag, stage = sys.argv[1], sys.argv[2]
    n = int(sys.argv[sys.argv.index('--n') + 1]) if '--n' in sys.argv else 24
    sol = close2.__dict__['placements'](tag)[0] if False else None
    sol, ppm = placements(tag)
    saved = json.load(open(P('data', f'{stage}_{tag}.json')))
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r].update(dE=v['dE'], dN=v['dN'], rot=v.get('rot', sol[r]['rot']),
                          gw=v.get('gw', sol[r]['gw']), gh=v.get('gh', sol[r]['gh']))
    recs = sorted(r for r in sol if sol[r].get('ok')); lab = SC.lines_of(sol, recs)
    mpp = C.MPP_TIE
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp, crop=0.95,
                                  log=lambda *_: None)
    # ridge maps per frame, for the ridge variant
    ridge = {}
    pairs = []
    for a, b in itertools.combinations(sorted(rend), 2):
        if lab[a] == lab[b]:
            continue
        ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
        x0 = max(xa, xb); y0 = max(ya, yb)
        x1 = min(xa + ia.shape[1], xb + ib.shape[1]); y1 = min(ya + ia.shape[0], yb + ib.shape[0])
        if x1 - x0 >= 256 and y1 - y0 >= 256:
            pairs.append((a, b, x0, y0, x1, y1))
    rng = np.random.default_rng(0)
    sample = [pairs[i] for i in rng.choice(len(pairs), min(n, len(pairs)), replace=False)]
    print(f"{tag} {stage}: {len(pairs)} sidelaps, sampling {len(sample)}")
    variants = [('raw  800m', False, 400), ('raw  512m', False, 256), ('ridge 800m', True, 400), ('ridge 512m', True, 256)]
    from collections import Counter
    summary = {v[0]: Counter() for v in variants}
    accepted = {v[0]: 0 for v in variants}
    for a, b, x0, y0, x1, y1 in sample:
        ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
        if a not in ridge:
            ridge[a] = gridval.ridge_full(ia.astype(np.float32), mpp)[0]
        if b not in ridge:
            ridge[b] = gridval.ridge_full(ib.astype(np.float32), mpp)[0]
        for name, use_ridge, win in variants:
            src_a = ridge[a] if use_ridge else ia; src_b = ridge[b] if use_ridge else ib
            offs = []
            step = win // 2
            for wy in range(y0, y1 - win + 1, step):
                for wx in range(x0, x1 - win + 1, step):
                    A = src_a[wy - ya:wy - ya + win, wx - xa:wx - xa + win].astype(np.float32)
                    B = src_b[wy - yb:wy - yb + win, wx - xb:wx - xb + win].astype(np.float32)
                    va = (ia[wy - ya:wy - ya + win, wx - xa:wx - xa + win] > 0).mean()
                    vb = (ib[wy - yb:wy - yb + win, wx - xb:wx - xb + win] > 0).mean()
                    if va < 0.85 or vb < 0.85:
                        summary[name]['coverage<0.85' if min(va, vb) >= 0.6 else 'coverage<0.6'] += 1; continue
                    if A.std() < (1e-4 if use_ridge else 5) or B.std() < (1e-4 if use_ridge else 5):
                        summary[name]['flat'] += 1; continue
                    dx, dy, sharp, edge = close2._pc(A, B)
                    if edge:
                        summary[name]['boundary'] += 1; continue
                    if sharp < 7.0:
                        summary[name]['weak(sharp<7)'] += 1; continue
                    if math.hypot(dx, dy) * mpp > 300:
                        summary[name]['>300m'] += 1; continue
                    summary[name]['ok'] += 1
                    offs.append((dx * mpp, -dy * mpp, sharp))
            if len(offs) >= 3:
                O = np.array(offs); c = np.median(O[:, :2], axis=0)
                d = np.hypot(O[:, 0] - c[0], O[:, 1] - c[1])
                if (d <= 12).sum() >= 3:
                    accepted[name] += 1
    print(f"\n  {'variant':11s} {'pairs accepted':>14s}   window outcomes")
    for name, _, _ in variants:
        print(f"  {name:11s} {accepted[name]:4d}/{len(sample):<9d}   {dict(summary[name])}")
    # sharpness distribution of the non-failed windows, raw vs ridge, to see if the
    # threshold is the issue
    print("\n  sharpness of windows that reached the peak test (raw 512m vs ridge 512m):")
    for name, use_ridge, win in (variants[1], variants[3]):
        sh = []
        for a, b, x0, y0, x1, y1 in sample[:8]:
            ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
            src_a = ridge[a] if use_ridge else ia; src_b = ridge[b] if use_ridge else ib
            for wy in range(y0, y1 - win + 1, win // 2):
                for wx in range(x0, x1 - win + 1, win // 2):
                    A = src_a[wy - ya:wy - ya + win, wx - xa:wx - xa + win].astype(np.float32)
                    B = src_b[wy - yb:wy - yb + win, wx - xb:wx - xb + win].astype(np.float32)
                    va = (ia[wy - ya:wy - ya + win, wx - xa:wx - xa + win] > 0).mean()
                    vb = (ib[wy - yb:wy - yb + win, wx - xb:wx - xb + win] > 0).mean()
                    if va < 0.85 or vb < 0.85 or A.std() < 1e-6 or B.std() < 1e-6:
                        continue
                    sh.append(close2._pc(A, B)[2])
        if sh:
            sh = np.array(sh)
            print(f"    {name}: n={len(sh)}  median {np.median(sh):.1f}  p75 {np.percentile(sh,75):.1f}  p90 {np.percentile(sh,90):.1f}  >=7: {(sh>=7).sum()}  >=5: {(sh>=5).sum()}")


if __name__ == '__main__':
    main()
