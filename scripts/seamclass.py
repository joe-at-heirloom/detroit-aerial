#!/usr/bin/env python
"""Split the seam metric into along-track and cross-line pairs.

A block of parallel flight lines can be perfectly closed along each line and still
have the lines disagree with each other -- along-track ties cannot see a per-line
scale error, only sidelap ties can. A single "median 0.19 m" over all pairs hides
that if along-track pairs dominate the count. This reports the two classes
separately, and also how many cross-line pairs exist geometrically versus how many
survived the confidence filter, so a class being silently dropped is visible.

No external reference is involved anywhere here.

Usage:  ./.venv/bin/python scripts/seamclass.py 1961 stageA
"""
import sys, os, json, math, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust
from validate import P
from rebuild import placements, SP
import close as C


def lines_of(sol, recs, gap=800.0):
    E = np.array([sol[r]['e'] for r in recs])
    order = np.argsort(E); Es = E[order]
    lab = {}; k = 0
    for i, oi in enumerate(order):
        if i > 0 and Es[i] - Es[i - 1] > gap:
            k += 1
        lab[recs[oi]] = k
    return lab


def main():
    tag, stage = sys.argv[1], sys.argv[2]
    sol, ppm = placements(tag)
    saved = json.load(open(P('data', f'{stage}_{tag}.json')))
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
            if 'rot' in v: sol[r]['rot'] = v['rot']
            if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
    recs = sorted([r for r in sol if sol[r].get('ok')])
    lab = lines_of(sol, recs)
    mpp = C.MPP_TIE
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp,
                                  crop=0.95, log=lambda *_: None)
    print(f"{tag} {stage}: {len(rend)} frames rendered, {max(lab.values())+1} flight lines")

    # geometric candidates: pairs whose rendered extents overlap by > 360 m both ways
    cand = {'along': 0, 'cross': 0}
    ext = {r: (x0, y0, x0 + im.shape[1], y0 + im.shape[0]) for r, (im, x0, y0) in rend.items()}
    for a, b in itertools.combinations(sorted(rend), 2):
        ax0, ay0, ax1, ay1 = ext[a]; bx0, by0, bx1, by1 = ext[b]
        ow = min(ax1, bx1) - max(ax0, bx0); oh = min(ay1, by1) - max(ay0, by0)
        if ow >= 180 and oh >= 180:
            cand['along' if lab[a] == lab[b] else 'cross'] += 1

    for thr in (1.10, 1.05):
        obs = blockadjust.relative_observations(rend, mpp, search_m=90.0, min_ratio=thr,
                                                log=lambda *_: None)
        cls = {'along': [], 'cross': []}
        for a, b, de, dn, rt in obs:
            cls['along' if lab[a] == lab[b] else 'cross'].append(math.hypot(de, dn))
        print(f"\n  min_ratio {thr}:")
        for k in ('along', 'cross'):
            v = np.array(cls[k])
            if len(v):
                print(f"    {k:5s}-line pairs: {len(v):3d} measured of {cand[k]:3d} geometric   "
                      f"median {np.median(v):5.1f}  p90 {np.percentile(v,90):5.1f}  "
                      f"max {v.max():5.1f} m   over 10 m {(v>10).sum()}")
            else:
                print(f"    {k:5s}-line pairs:   0 measured of {cand[k]:3d} geometric")
        if thr == 1.10:
            # where along the strip are the cross-line disagreements? scale error
            # shows up as disagreement growing toward the block ends
            xs = [(0.5 * (sol[a]['n'] + sol[b]['n']), math.hypot(de, dn), de, dn)
                  for a, b, de, dn, rt in obs if lab[a] != lab[b]]
            if xs:
                xs.sort()
                print("    cross-line disagreement by position along the strip (N, |d|, dE, dN):")
                for n_, m_, de, dn in xs:
                    print(f"      N {n_:+8.0f}   {m_:5.1f}   ({de:+6.1f},{dn:+6.1f})")


if __name__ == '__main__':
    main()
