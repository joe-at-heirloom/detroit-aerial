#!/usr/bin/env python
"""Composite a block from a saved close.py solution, then measure it both ways.

Stage A is slow -- it renders every frame several times -- so the solution is saved
and applied separately. Both numbers get reported, because either one alone is
misleading: internal consistency without absolute placement is a beautiful mosaic in
the wrong place, and absolute placement without internal consistency is what this
project shipped for two sessions.
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust, dtmap
import validate
from validate import P, MPP, report
from rebuild import placements, composite, SP
import close as C


def run(tag):
    t0 = time.time()
    sol, ppm = placements(tag)
    which = 'closed'
    for a in sys.argv[1:]:
        if a in ('--stageA', '--closed'):
            which = a[2:]
    saved = json.load(open(P('data', f'{which}_{tag}.json')))
    n = 0
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
            if 'rot' in v: sol[r]['rot'] = v['rot']
            if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
            n += 1
    print(f"\n===== {tag}: applying the {which} solution ({n} frames) =====", flush=True)
    composite(tag, sol, ppm, 'closed')

    print("  internal consistency (no external reference of any kind):", flush=True)
    mpp = C.MPP_TIE
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp,
                                  crop=0.95, log=lambda *_: None)
    v, obs = C.seam_stats(rend, mpp)
    del rend
    if v is not None:
        print(f"    {len(obs)} overlapping pairs: median {np.median(v):5.2f}  "
              f"p90 {np.percentile(v,90):5.1f}  max {v.max():5.1f} m   "
              f"over 10 m: {(v>10).sum()}", flush=True)

    print("  absolute placement against modern imagery:", flush=True)
    arr, mod, bbox = validate.load(tag, 'closed')
    t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_closed_{int(t0)}',
                               log=lambda *_: None)
    pz = gridval.Prior([], MPP)
    report(gridval.grid_hier_prep(t, prior=pz), f'    {tag} 3.6 x 2.0 km cells')
    report(gridval.grid_hier_prep(t, NY=32, NX=6, prior=pz), f'    {tag} 1.8 x 1.0 km cells')
    print(f"  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    for tg in [a for a in sys.argv[1:] if not a.startswith('--')]:
        run(tg)
