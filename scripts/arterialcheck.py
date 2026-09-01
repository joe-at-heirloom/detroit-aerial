#!/usr/bin/env python
"""Cross-check a measured displacement field against the mile-grid arterials.

The field was measured by ridge-vs-modern-imagery correlation, which can alias on
Detroit's 97.5 m residential lattice wherever the regional prior is more than 40 m
off. Adjacent cells jumping by ~90 m is what that looks like. The arterials
repeat every 1609 m, so correlating against them cannot alias inside +/-300 m --
it is blunter (a ~5 m floor, a known systematic bias, and it needs arterials in
two directions), but it answers the one question that matters: is the ridge field
real, or is part of it lattice aliasing?

The arterial reference's own bias is measured here too, by running the identical
check on the MODERN raster masked to the historical footprint. Modern imagery is
correct by definition, so whatever that reports is the arterial method's error,
and it is subtracted.

Usage:  ./.venv/bin/python scripts/arterialcheck.py 1961 closedA
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, dtmap
import validate
from validate import P, MPP


def main():
    tag, label = sys.argv[1], sys.argv[2]
    F = json.load(open(P('data', f'field_{tag}_{label}.json')))
    arr, mod, bbox = validate.load(tag, label)
    minE, maxN = F['minE'], F['maxN']
    mask = arr > 0
    NY, NX = 16, 3
    print(f"{tag} {label}: arterial cross-check on the {NY}x{NX} grid (search +/-300 m)", flush=True)
    ctl = gridval.arterial_grid((mod * mask).astype(np.uint8), MPP, minE, maxN,
                                NY=NY, NX=NX, search_m=300.0, min_div=0.20)
    C = {(c['row'], c['col']): c for c in ctl if 'skip' not in c}
    good = [c for c in C.values() if c['ratio'] >= 1.15 and c['mag'] < 40]
    if len(good) < 4:
        print("  arterial control did not lock enough cells on modern imagery; cannot use it")
        return
    bias = np.median([[c['dE'], c['dN']] for c in good], axis=0)
    print(f"  arterial reference bias (from modern, {len(good)} cells): "
          f"dE {bias[0]:+.1f}  dN {bias[1]:+.1f} m", flush=True)
    art = gridval.arterial_grid(arr, MPP, minE, maxN, NY=NY, NX=NX, search_m=300.0, min_div=0.20)
    A = {(c['row'], c['col']): c for c in art if 'skip' not in c}
    R = {(c['row'], c['col']): c for c in F['cells_16x3'] if 'skip' not in c and not c.get('pegged')}
    print(f"\n  cell     ridge-vs-modern        arterial (bias-removed)    diff    art.ratio div   verdict")
    agree = []; rows = []
    for k in sorted(set(A) | set(R)):
        a = A.get(k); r = R.get(k)
        ctl_ok = k in C and C[k]['ratio'] >= 1.15 and C[k]['mag'] < 40
        if not a or a['ratio'] < 1.12 or not ctl_ok:
            s_art = '     (no lock)     '
        else:
            adE, adN = a['dE'] - bias[0], a['dN'] - bias[1]
            s_art = f"({adE:+6.1f},{adN:+6.1f})     "
        s_r = f"({r['dE']:+6.1f},{r['dN']:+6.1f})" if r and r['ratio'] >= 1.15 else "   (no lock)   "
        verdict = ''
        if a and r and a['ratio'] >= 1.12 and r['ratio'] >= 1.15 and ctl_ok:
            d = math.hypot(a['dE'] - bias[0] - r['dE'], a['dN'] - bias[1] - r['dN'])
            verdict = 'AGREE' if d < 25 else ('ridge ALIASED?' if abs(d - 97.5) < 30 or abs(d-138)<30 else 'DISAGREE')
            agree.append(d < 25)
            rows.append((k, r['dE'], r['dN'], a['dE'] - bias[0], a['dN'] - bias[1], d))
            print(f"  {k[0]:2d},{k[1]}   {s_r:16s}   {s_art}  {d:5.1f}     {a['ratio']:4.2f}  {a['div']:.2f}  {verdict}")
        else:
            print(f"  {k[0]:2d},{k[1]}   {s_r:16s}   {s_art}          {a['ratio'] if a else 0:4.2f}  {a['div'] if a else 0:.2f}")
    if agree:
        print(f"\n  {sum(agree)}/{len(agree)} comparable cells agree within 25 m")
        m = np.array([x[5] for x in rows])
        print(f"  ridge minus arterial: median {np.median(m):.1f}  p90 {np.percentile(m,90):.1f} m")
    json.dump(dict(bias=bias.tolist(), arterial=art, control=ctl),
              open(P('data', f'artcheck_{tag}_{label}.json'), 'w'))


if __name__ == '__main__':
    main()
