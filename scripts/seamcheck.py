#!/usr/bin/env python
"""Is the mosaic internally consistent, or just correct on average?

Every accuracy figure in this project so far answers one question: does this AREA
sit in the right place against modern imagery. None of them answers the question a
viewer actually notices first: does a road stay straight across the join between
two negatives.

Those are different, and the per-frame stage is exactly what pulls them apart. It
moves each negative independently to best-match modern imagery, deliberately
overriding the relative tie constraints that made the block internally consistent.
Each frame can then be individually right on average while two neighbours disagree
along their shared edge -- and a road crossing that edge kinks. A 3.6 x 2.0 km
grid cell contains several frames, so four of them disagreeing inside one cell
still reports a small median.

This measures the disagreement directly: render each pair of overlapping frames
into map space with their final placements, and correlate the two in the region
they share. Perfectly consistent frames give zero. What comes back is the size of
the kink a viewer sees at that seam.

Usage:  ./.venv/bin/python scripts/seamcheck.py 1949 [--suffix final] [--max 120]
"""
import sys, os, json, math, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, frameadjust, dtmap
import validate
from validate import P, MPP
from rebuild import placements, SP


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def load_placements(tag, use_adj=True):
    """Frame placements as actually composited: bundle solution plus whatever the
    per-frame stage decided."""
    sol, ppm = placements(tag)
    fp = P('data', f'frameadj_{tag}.json')
    if use_adj and os.path.exists(fp):
        fx = json.load(open(fp))
        for r in list(sol):
            v = fx.get(str(r)) or fx.get(r)
            if v:
                sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
                if 'rot' in v:
                    sol[r]['rot'] = v['rot']       # stored absolute
    return sol, ppm


def main():
    tag = [a for a in sys.argv[1:] if not a.startswith('--')][0]
    maxpairs = int(arg('--max', 120))
    mpp = float(arg('--mpp', 2.0))
    for use_adj, label in ((False, 'BEFORE per-frame (bundle only)'),
                           (True, 'AFTER per-frame correction')):
        sol, ppm = load_placements(tag, use_adj)
        recs = [r for r in sol if sol[r].get('ok')]
        # frame centres and radius, to find genuinely overlapping pairs
        C = {r: (sol[r]['e'] - sol[r].get('dE', 0.0),
                 sol[r]['n'] - sol[r].get('dN', 0.0)) for r in recs}
        rad = {r: math.hypot(sol[r]['gw'], sol[r]['gh']) / 2 for r in recs}
        pairs = []
        for a, b in itertools.combinations(sorted(recs), 2):
            d = math.hypot(C[a][0] - C[b][0], C[a][1] - C[b][1])
            if 200 < d < 0.75 * (rad[a] + rad[b]):
                pairs.append((d, a, b))
        pairs.sort()
        pairs = pairs[:maxpairs]
        if not pairs:
            print(f"  {label}: no overlapping pairs found"); continue

        # a map grid big enough for any pair
        E = [C[r][0] for r in recs]; N = [C[r][1] for r in recs]
        gw = max(rad.values())
        minE = min(E) - gw; maxN = max(N) + gw
        W = int((max(E) + gw - minE) / mpp); H = int((maxN - (min(N) - gw)) / mpp)
        res = []
        for d, a, b in pairs:
            ims = []
            for r in (a, b):
                out = frameadjust.render_frame(r, sol[r], f"{SP}/fullres",
                                               minE, maxN, W, H, mpp, crop=0.92)
                ims.append(out)
            if any(o is None for o in ims):
                continue
            (ia, xa, ya), (ib, xb, yb) = ims
            # common window
            x0 = max(xa, xb); y0 = max(ya, yb)
            x1 = min(xa + ia.shape[1], xb + ib.shape[1])
            y1 = min(ya + ia.shape[0], yb + ib.shape[0])
            if x1 - x0 < 160 or y1 - y0 < 160:
                continue
            A = ia[y0 - ya:y1 - ya, x0 - xa:x1 - xa]
            B = ib[y0 - yb:y1 - yb, x0 - xb:x1 - xb]
            if (A > 0).mean() < 0.55 or (B > 0).mean() < 0.55:
                continue
            ra, va = gridval.ridge_full(A.astype(np.float32), mpp)
            rb, vb = gridval.ridge_full(B.astype(np.float32), mpp)
            if ra.std() < 1e-9 or rb.std() < 1e-9:
                continue
            m = gridval.match_cell((ra * va).astype(np.float32), va.astype(np.float32),
                                   (rb * vb).astype(np.float32), vb.astype(np.float32),
                                   mpp, 60.0, 'mncc')
            if m is None or m['ratio'] < 1.10:
                continue
            res.append((m['mag'], a, b, m['ratio']))
        if not res:
            print(f"  {label}: nothing measurable"); continue
        v = np.array([r[0] for r in res])
        print(f"  {label}: {len(res)} overlapping pairs   "
              f"median {np.median(v):5.1f}  p90 {np.percentile(v,90):5.1f}  "
              f"max {v.max():5.1f} m   over 10 m: {(v>10).sum()}", flush=True)
        for mag, a, b, rt in sorted(res, reverse=True)[:6]:
            print(f"      worst: {a} vs {b}  {mag:5.1f} m (ratio {rt:.2f})")


if __name__ == '__main__':
    main()
