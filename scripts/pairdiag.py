#!/usr/bin/env python
"""Per-pair tie yield and disagreement on a closed block, worst first.

1949 closes along-track to 2 m and stays at 20 m cross-line. Is that a handful of
bad pairs, or the whole sidelap? And do the bad pairs yield few windows -- so the
solver never had the constraint -- or many windows that disagree?

Usage:  ./.venv/bin/python scripts/pairdiag.py 1949 stageA3
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import blockadjust
from validate import P
from rebuild import placements, SP
import close as C, close3, seamclass as SC


def main():
    tag, stage = sys.argv[1], sys.argv[2]
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
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp, crop=0.95, log=lambda *_: None)
    # loose consensus so we also see what a strict filter would have thrown away
    ties = close3.tie_windows(rend, lab, mpp, agree_m=60.0, min_windows=2, log=print)
    del rend
    from collections import defaultdict
    per = defaultdict(list)
    for t in ties:
        per[(t['a'], t['b'])].append((t['dE'], t['dN']))
    rows = []
    for (a, b), v in per.items():
        v = np.array(v); c = np.median(v, axis=0); d = np.hypot(v[:, 0] - c[0], v[:, 1] - c[1])
        rows.append(dict(a=a, b=b, cross=lab[a] != lab[b], n=len(v), off=float(math.hypot(*c)),
                         spread=float(np.percentile(d, 75)), nN=(sol[a]['n'] + sol[b]['n']) / 2))
    for cls in (True, False):
        rr = sorted([r for r in rows if r['cross'] == cls], key=lambda r: -r['off'])
        v = np.array([r['off'] for r in rr]); n = np.array([r['n'] for r in rr])
        print(f"\n{'CROSS-LINE' if cls else 'ALONG-TRACK'} pairs: {len(rr)}  offset median {np.median(v):.1f}  p90 {np.percentile(v,90):.1f} m   windows/pair median {np.median(n):.0f}")
        print(f"  {'pair':14s} {'lines':6s} {'N':>7s} {'windows':>7s} {'offset':>7s} {'spread':>7s}")
        for r in rr[:14]:
            print(f"  {r['a']:>6s}/{r['b']:<6s} {lab[r['a']]}-{lab[r['b']]}   {r['nN']:+7.0f} {r['n']:7d} {r['off']:7.1f} {r['spread']:7.1f}")
    json.dump(rows, open(P('data', f'pairdiag_{tag}_{stage}.json'), 'w'))


if __name__ == '__main__':
    main()
