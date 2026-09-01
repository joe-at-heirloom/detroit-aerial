#!/usr/bin/env python
"""Place a closed block: translation, plus a gentle along-strip bend if one is measured.

Stage B was written three times and got it wrong three times -- wrong coordinate
frame, similarity applied to centres but not frames, and a fine search with no
alias-proof prior. Each produced a plausible small number rather than a visible
failure. The last version fitted a 0.78% scale that was not real scale at all but
the four-parameter fit absorbing residual shape error, and applying it took internal
consistency from 8.3 m to 21.2 m -- worse than the bundle it started from.

So: no more bespoke placement. Composite the closed block, measure where it sits
with `gridval`'s alias-proof grid measurement -- the one instrument here that has
been calibrated against planted shifts on every block -- and move the whole thing by
one translation.

A translation is the only transform that CANNOT change internal consistency. Every
frame moves identically, so no seam can open. Scale and rotation are deliberately
not fitted: on this data they absorb error rather than remove it, and the evidence
that they were real was never there.

Usage:  ./.venv/bin/python scripts/placeshift.py 1961 [--stage stageA]
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


def load_stage(tag, stage):
    sol, ppm = placements(tag)
    saved = json.load(open(P('data', f'{stage}_{tag}.json')))
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
            if 'rot' in v: sol[r]['rot'] = v['rot']
            if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
    return sol, ppm


def seam(sol, tag, label=''):
    mpp = C.MPP_TIE
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp,
                                  crop=0.95, log=lambda *_: None)
    v, obs = C.seam_stats(rend, mpp)
    del rend
    if v is None:
        print(f"    {label}: no overlap observations"); return
    print(f"    {label}: {len(obs)} pairs, median {np.median(v):5.2f}  "
          f"p90 {np.percentile(v,90):5.1f}  max {v.max():5.1f} m  over 10 m {(v>10).sum()}",
          flush=True)


def main():
    tag = [a for a in sys.argv[1:] if not a.startswith('--')][0]
    stage = sys.argv[sys.argv.index('--stage') + 1] if '--stage' in sys.argv else 'stageA'
    t0 = time.time()
    sol, ppm = load_stage(tag, stage)
    print(f"\n===== place {tag} by translation (from {stage}) =====", flush=True)
    print("  internal consistency of the closed block, before placing:", flush=True)
    seam(sol, tag, 'closed')

    composite(tag, sol, ppm, 'closedA')
    arr, mod, bbox = validate.load(tag, 'closedA')
    t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_closedA_{int(t0)}',
                               log=lambda *_: None)
    pri = gridval.prior_for(t, MPP)
    cells = gridval.grid_hier_prep(t, NY=16, NX=3, prior=pri)
    good = [c for c in cells if 'skip' not in c and not c.get('pegged')
            and c['ratio'] >= 1.15]
    if len(good) < 6:
        print(f"  only {len(good)} cells measurable; cannot place"); return
    dE = np.array([c['dE'] for c in good]); dN = np.array([c['dN'] for c in good])
    w = np.array([min(max(c['ratio'] - 1.0, 0.02) / 0.5, 3.0) for c in good])
    keep = np.ones(len(dE), bool)
    for _ in range(6):
        cx = np.average(dE[keep], weights=w[keep]); cy = np.average(dN[keep], weights=w[keep])
        r = np.hypot(dE - cx, dN - cy)
        s = max(np.median(r[keep]) * 1.4826, 4.0)
        keep = r < 3.0 * s
        if keep.sum() < 6:
            break
    print(f"  {len(good)} cells measured, {int(keep.sum())} kept; "
          f"block sits {cx:+.1f}, {cy:+.1f} m from where it belongs", flush=True)
    # is the offset uniform, or does it vary along the strip? if it varies, a
    # translation is not enough and that must be said rather than fitted away
    yy = np.array([c['row'] for c in good])[keep]
    if len(set(yy.tolist())) > 4:
        aE = np.polyfit(yy, dE[keep], 1)[0]; aN = np.polyfit(yy, dN[keep], 1)[0]
        span = yy.max() - yy.min()
        print(f"  along-strip trend: {aE*span:+.1f} m east, {aN*span:+.1f} m north "
              f"end to end (a translation cannot remove this)", flush=True)

    # A closed block is locally perfect and globally bent. Local pairwise constraints
    # fix local shape exactly but accumulate drift down a 50-link chain, so the strip
    # ends up a chain of perfectly-fitted links curving away from reality. Measured on
    # 1961 after closing: 94 m east and 180 m north end to end.
    #
    # Removing that needs the relative geometry to change slightly -- but only by the
    # per-step amount, roughly the total drift divided by the number of frames, a few
    # metres. That is a good trade for 180 m of absolute error, and unlike a free
    # per-frame correction it is a SMOOTH function of position along the strip, so
    # neighbours still move almost identically. Degree stays low deliberately: this
    # is here to remove measured drift, not to soak up residual error, which is what
    # the deleted warp and the 0.78% scale fit were doing.
    deg = int(sys.argv[sys.argv.index('--bend') + 1]) if '--bend' in sys.argv else 1
    yy_all = np.array([c['row'] for c in good], float)
    fitE = np.polyfit(yy_all[keep], dE[keep], deg)
    fitN = np.polyfit(yy_all[keep], dN[keep], deg)
    # map a frame's north coordinate onto the same row index the cells used
    S_, W_, N_, E_ = bbox
    maxN_m = (N_ - dtmap.LAT0) * dtmap.MLAT
    rows = 16
    cellH = (N_ - S_) * dtmap.MLAT / rows
    if deg <= 0:
        print("  applying translation only", flush=True)
    else:
        print(f"  applying translation plus a degree-{deg} along-strip bend", flush=True)
    for r in sol:
        if not sol[r].get('ok'):
            continue
        n_m = sol[r]['n'] - sol[r].get('dN', 0.0)
        row = (maxN_m - n_m) / cellH
        row = float(np.clip(row, yy_all.min(), yy_all.max()))
        sol[r]['dE'] += float(np.polyval(fitE, row))
        sol[r]['dN'] += float(np.polyval(fitN, row))
    json.dump({r: dict(dE=sol[r]['dE'], dN=sol[r]['dN'], rot=sol[r]['rot'],
                       gw=sol[r]['gw'], gh=sol[r]['gh'])
               for r in sol if sol[r].get('ok')},
              open(P('data', f'placed_{tag}.json'), 'w'))
    composite(tag, sol, ppm, 'placed')
    print("  RESULT", flush=True)
    seam(sol, tag, 'internal consistency')
    arr2, mod2, _ = validate.load(tag, 'placed')
    t2 = gridval.prepare_cached(arr2, mod2, MPP, f'{tag}_placed_{int(t0)}',
                                log=lambda *_: None)
    pri2 = gridval.prior_for(t2, MPP, log=lambda *_: None)
    report(gridval.grid_hier_prep(t2, prior=pri2), '    absolute 3.6 x 2.0 km')
    report(gridval.grid_hier_prep(t2, NY=32, NX=6, prior=pri2), '    absolute 1.8 x 1.0 km')
    print(f"  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
