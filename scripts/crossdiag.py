#!/usr/bin/env python
"""Why do sidelap pairs fail to enter Stage A?

225 of 241 geometric cross-line pairs on 1961 were dropped before the solve. Each
dropped pair is a constraint the block needed. This reports, for every cross-line
pair regardless of threshold: overlap size, valid fraction on each side, the match
ratio, whether the peak hit the search edge, and the offset found -- so the failure
is attributable to weak peaks, to disagreement beyond the search, or to the crop.

Usage:  ./.venv/bin/python scripts/crossdiag.py 1961 stageA [--search 150]
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
import seamclass as SC


def main():
    tag, stage = sys.argv[1], sys.argv[2]
    search = float(sys.argv[sys.argv.index('--search') + 1]) if '--search' in sys.argv else 150.0
    sol, ppm = placements(tag)
    saved = json.load(open(P('data', f'{stage}_{tag}.json')))
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
            if 'rot' in v: sol[r]['rot'] = v['rot']
            if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
    recs = sorted([r for r in sol if sol[r].get('ok')])
    lab = SC.lines_of(sol, recs)
    mpp = C.MPP_TIE
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp,
                                  crop=0.95, log=lambda *_: None)
    ext = {r: (x0, y0, x0 + im.shape[1], y0 + im.shape[0]) for r, (im, x0, y0) in rend.items()}
    rows = []
    for a, b in itertools.combinations(sorted(rend), 2):
        if lab[a] == lab[b]:
            continue
        ax0, ay0, ax1, ay1 = ext[a]; bx0, by0, bx1, by1 = ext[b]
        x0 = max(ax0, bx0); y0 = max(ay0, by0); x1 = min(ax1, bx1); y1 = min(ay1, by1)
        if x1 - x0 < 180 or y1 - y0 < 180:
            continue
        ia, xa, ya = rend[a]; ib, xb, yb = rend[b]
        A = ia[y0 - ya:y1 - ya, x0 - xa:x1 - xa]; B = ib[y0 - yb:y1 - yb, x0 - xb:x1 - xb]
        va = float((A > 0).mean()); vb = float((B > 0).mean())
        rec = dict(a=a, b=b, ow=(x1 - x0) * mpp, oh=(y1 - y0) * mpp, va=va, vb=vb,
                   n=0.5 * (sol[a]['n'] + sol[b]['n']))
        if va < 0.3 or vb < 0.3:
            rec['why'] = 'valid'; rows.append(rec); continue
        ra, ma = gridval.ridge_full(A.astype(np.float32), mpp)
        rb, mb = gridval.ridge_full(B.astype(np.float32), mpp)
        if ra.std() < 1e-9 or rb.std() < 1e-9:
            rec['why'] = 'flat'; rows.append(rec); continue
        m = gridval.match_cell((ra * ma).astype(np.float32), ma.astype(np.float32),
                               (rb * mb).astype(np.float32), mb.astype(np.float32),
                               mpp, search, 'mncc')
        if m is None:
            rec['why'] = 'nopeak'; rows.append(rec); continue
        rec.update(dE=m['dE'], dN=m['dN'], mag=m['mag'], ratio=m['ratio'], edge=m['edge'])
        rec['why'] = 'edge' if m['edge'] else ('weak' if m['ratio'] < 1.10 else 'ok')
        rows.append(rec)
    from collections import Counter
    print(f"{tag} {stage}: {len(rows)} cross-line pairs with >=360 m overlap, search +/-{search:.0f} m")
    print("  outcome:", dict(Counter(r['why'] for r in rows)))
    ok = [r for r in rows if r['why'] == 'ok']
    weak = [r for r in rows if r['why'] == 'weak']
    if ok:
        v = np.array([r['mag'] for r in ok])
        print(f"  ok   (ratio>=1.10): {len(ok)}  disagreement median {np.median(v):.1f}  p90 {np.percentile(v,90):.1f}  max {v.max():.1f}")
    if weak:
        v = np.array([r['mag'] for r in weak]); rt = np.array([r['ratio'] for r in weak])
        print(f"  weak (ratio<1.10):  {len(weak)}  ratio median {np.median(rt):.3f}  their offsets median {np.median(v):.1f}  p90 {np.percentile(v,90):.1f}")
    ov = np.array([min(r['ow'], r['oh']) for r in rows]); va = np.array([min(r['va'], r['vb']) for r in rows])
    print(f"  overlap narrow side: median {np.median(ov):.0f} m   valid fraction median {np.median(va):.2f}")
    print("\n  per-pair (N, overlap w x h m, valid a/b, ratio, |d|, why):")
    for r in sorted(rows, key=lambda r: -r['n']):
        print(f"    N {r['n']:+7.0f}  {r['ow']:5.0f}x{r['oh']:5.0f}  {r['va']:.2f}/{r['vb']:.2f}  "
              f"{r.get('ratio', 0):5.2f}  {r.get('mag', 0):6.1f}  {r['why']}")
    json.dump(rows, open(P('data', f'crossdiag_{tag}_{stage}.json'), 'w'))


if __name__ == '__main__':
    main()
