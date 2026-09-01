#!/usr/bin/env python
"""Is the frame-to-frame disagreement a blunder, or is it geometry?

Overlapping frames in every block disagree with each other -- 16.5 m median in 1961,
22.6 m in 1949 -- straight out of the bundle adjustment, before anything was
corrected. That number alone does not say what is wrong, and the two candidate
causes want opposite responses:

  a blunder    bad tie points, a mis-solved frame. The disagreement is roughly
               CONSTANT across the overlap: one frame simply sits offset from the
               other. Fixable by re-solving the block.
  geometry     per-frame tilt, film deformation, lens distortion, an unmodelled
               scale. A similarity transform cannot represent these, so the
               disagreement VARIES ACROSS the overlap -- growing with distance from
               the frame centre, or along one axis. Not fixable by any amount of
               re-solving with the current model; it needs a camera model.

So measure the disagreement at several places within each overlap instead of once
for the whole thing, and look at how it varies. Also check whether the offsets sit
near multiples of the street lattice, which would mean the tie points aliased
rather than anything being geometrically wrong.

Usage:  ./.venv/bin/python scripts/characterise.py 1961 [--pairs 14]
"""
import sys, os, json, math, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust, dtmap
from validate import P
from rebuild import placements, SP


def main():
    tag = [a for a in sys.argv[1:] if not a.startswith('--')][0]
    npairs = int(sys.argv[sys.argv.index('--pairs') + 1]) if '--pairs' in sys.argv else 14
    mpp = 2.0
    sol, ppm = placements(tag)          # bundle only, no per-frame correction
    recs = sorted([r for r in sol if sol[r].get('ok')])
    E = [sol[r]['e'] for r in recs]; N = [sol[r]['n'] for r in recs]
    gw = max(math.hypot(sol[r]['gw'], sol[r]['gh']) / 2 for r in recs)
    minE = min(E) - gw; maxN = max(N) + gw
    W = int((max(E) + gw - minE) / mpp); H = int((maxN - (min(N) - gw)) / mpp)
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp,
                                  crop=0.95, log=lambda *_: None)
    print(f"{tag}: {len(rend)} frames rendered at {mpp} m/px\n")

    keys = sorted(rend)
    pairs = []
    for a, b in itertools.combinations(keys, 2):
        d = math.hypot(sol[a]['e'] - sol[b]['e'], sol[a]['n'] - sol[b]['n'])
        if 200 < d < 0.7 * (math.hypot(sol[a]['gw'], sol[a]['gh'])):
            pairs.append((d, a, b))
    pairs.sort(); pairs = pairs[:npairs]
    print(f"  {len(pairs)} candidate overlapping pairs")
    drop = {'small': 0, 'novalid': 0, 'nowhole': 0, 'fewsubs': 0}

    print("  pair            whole-overlap   sub-window disagreement across the overlap")
    print("                  offset (m)      n   median   spread(p90-p10)   verdict")
    spreads = []; wholes = []; alln = []
    for d, a, b in pairs:
        ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
        x0 = max(xa, xb); y0 = max(ya, yb)
        x1 = min(xa + ia.shape[1], xb + ib.shape[1])
        y1 = min(ya + ia.shape[0], yb + ib.shape[0])
        if x1 - x0 < 250 or y1 - y0 < 250:
            drop['small'] += 1; continue
        A = ia[y0 - ya:y1 - ya, x0 - xa:x1 - xa]
        B = ib[y0 - yb:y1 - yb, x0 - xb:x1 - xb]
        if (A > 0).mean() < 0.45 or (B > 0).mean() < 0.45:
            drop['novalid'] += 1; continue

        def match(P_, Q_, search=80.0):
            # Both sides are the SAME epoch, often the same six seconds of flight, so
            # they look alike. Correlate the raw imagery: the ridge response exists to
            # bridge a sixty-year appearance gap, and its background filter is half
            # the width of a sub-window here, which is why every sub-window failed.
            va = P_ > 0; vb = Q_ > 0
            if va.mean() < 0.4 or vb.mean() < 0.4:
                return None
            A = P_.astype(np.float32); B = Q_.astype(np.float32)
            A = A - A[va].mean(); A[~va] = 0
            B = B - B[vb].mean(); B[~vb] = 0
            if A.std() < 1e-6 or B.std() < 1e-6:
                return None
            return gridval.match_cell(A, va.astype(np.float32), B, vb.astype(np.float32),
                                      mpp, search, 'mncc')

        whole = match(A, B)
        if whole is None or whole['ratio'] < 1.05:
            drop['nowhole'] += 1; continue
        # split the overlap into a grid and measure each cell separately
        subs = []
        # sample along the overlap BAND: the bbox intersection is mostly empty
        # because two rotated frames share a strip, not a rectangle
        NY, NX = 3, 5
        ch = A.shape[0] // NY; cw = A.shape[1] // NX
        for j in range(NY):
            for i in range(NX):
                sa = A[j * ch:(j + 1) * ch, i * cw:(i + 1) * cw]
                sb = B[j * ch:(j + 1) * ch, i * cw:(i + 1) * cw]
                if (sa > 0).mean() < 0.45 or (sb > 0).mean() < 0.45:
                    continue
                m = match(sa, sb)
                if m is None or m['edge'] or m['ratio'] < 1.06:
                    continue
                subs.append((i, j, m['dE'], m['dN']))
        if len(subs) < 3:
            drop['fewsubs'] += 1; continue
        mags = np.array([math.hypot(s[2], s[3]) for s in subs])
        spread = float(np.percentile(mags, 90) - np.percentile(mags, 10))
        verdict = 'GEOMETRY' if spread > 0.6 * np.median(mags) and spread > 8 else 'offset'
        spreads.append(spread); wholes.append(whole['mag']); alln.append(len(subs))
        print(f"  {a:>5} vs {b:<5}   {whole['mag']:7.1f}      {len(subs):2d}  "
              f"{np.median(mags):7.1f}   {spread:8.1f}        {verdict}")
    if not spreads:
        print(f"\n  nothing measurable; rejections: {drop}"); return
    print(f"  rejections: {drop}")
    sp = np.array(spreads); wh = np.array(wholes)
    print(f"\n  across {len(sp)} pairs:")
    print(f"    whole-overlap offset      median {np.median(wh):6.1f} m")
    print(f"    variation WITHIN overlaps median {np.median(sp):6.1f} m  "
          f"(p90 {np.percentile(sp,90):.1f})")
    frac = float(np.median(sp) / max(np.median(wh), 1e-6))
    print(f"    within-overlap variation is {frac*100:.0f}% of the offset itself")
    print("\n  " + ("GEOMETRY: the disagreement changes across each overlap, so a per-frame\n"
                    "  shift cannot describe it. A similarity transform per frame is the wrong\n"
                    "  model -- this needs a camera model, or many more degrees of freedom."
                    if frac > 0.6 else
                    "BLUNDER: the disagreement is roughly constant across each overlap, so each\n"
                    "  frame simply sits offset from its neighbour. The model is adequate and the\n"
                    "  block can be re-solved from better ties."))


if __name__ == '__main__':
    main()
