#!/usr/bin/env python
"""Build a solve plan for a whole roll of film, including frames nobody has placed.

The old catalogue only ever held the frames a previous pass had already located,
so extending a block means solving negatives with no position at all. Two things
make that tractable, both from the roll itself:

  * Consecutive exposure numbers are consecutive shots along a flight line and
    overlap about 60%. That is a reliable pair without any position.
  * Within one run of consecutive numbers, position is very nearly linear in
    number -- on 1949 roll ha-17 the spacing holds at 1200-1400 m a frame -- so a
    frame a few numbers off a solved run can be extrapolated to a few hundred
    metres, which is far inside the pair-plausibility radius.

Frames too far from any solved run get no position; they are still extracted and
matched through the adjacency pairs, and the reconstruction places them.

Usage:
  ./.venv/bin/python scripts/rollplan.py 1949 --roll ha-17 [--adj 3] [--reach 6]
                                         [--out /tmp/plan_1949_ha17.json]
"""
import sys, os, json
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts')); sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
from validate import P
from rewarp import SP, BUNDLE
import dtmap


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def main():
    tag = sys.argv[1]
    roll = arg('--roll'); adj = int(arg('--adj', 3)); reach = int(arg('--reach', 6))
    county = arg('--county', 'Wayne County')
    out = arg('--out', f"/tmp/plan_{tag}_{roll}.json")

    cat = [i for i in json.load(open(P('data', 'dte_catalogue.json')))['items']
           if i['county'] == county and i['year'] == tag and i['section']
           and i['section'].rsplit('-', 1)[0] == roll]
    num_of = {i['pointer']: int(i['section'].rsplit('-', 1)[1]) for i in cat}
    ptr_of = {v: k for k, v in num_of.items()}
    print(f"{tag} roll {roll}: {len(cat)} frames in the catalogue, numbers "
          f"{min(num_of.values())}..{max(num_of.values())}")

    # positions already solved, in the dtmap frame
    b = json.load(open(P('data', BUNDLE[tag][0])))[tag]
    c = json.load(open(P('data', f'calib{tag}.json')))
    solved = {}
    for r in b['frames']:
        if int(r) not in num_of:
            continue
        e, n = b['pos'][r]
        lat = c['lat0'] + n / c['MLAT']; lon = c['lon0'] + e / c['MLON']
        solved[num_of[int(r)]] = ((lon - dtmap.LON0) * dtmap.MLON, (lat - dtmap.LAT0) * dtmap.MLAT)
    print(f"  {len(solved)} of them already have a solved position")

    # runs of consecutive numbers whose spacing looks like one flight line
    runs = []
    if solved:
        ks = sorted(solved); cur = [ks[0]]
        for a, b_ in zip(ks, ks[1:]):
            step = float(np.hypot(*(np.array(solved[b_]) - np.array(solved[a])))) / (b_ - a)
            if b_ - a > reach or not (900 < step < 1800):
                runs.append(cur); cur = []
            cur.append(b_)
        runs.append(cur)
        print("  runs: " + ', '.join(f"{r[0]}-{r[-1]}({len(r)})" for r in runs))

    pos, interp = {}, 0
    for num, ptr in sorted(ptr_of.items()):
        if num in solved:
            pos[str(ptr)] = [round(solved[num][0], 1), round(solved[num][1], 1)]
            continue
        best = None
        for r in runs:
            near = [x for x in r if abs(x - num) <= reach]
            if len(near) < 2:
                continue
            A = np.array([[x, 1.0] for x in near])
            fe = np.linalg.lstsq(A, np.array([solved[x][0] for x in near]), rcond=None)[0]
            fn = np.linalg.lstsq(A, np.array([solved[x][1] for x in near]), rcond=None)[0]
            d = min(abs(x - num) for x in near)
            if best is None or d < best[0]:
                best = (d, (fe[0] * num + fe[1], fn[0] * num + fn[1]))
        if best:
            pos[str(ptr)] = [round(best[1][0], 1), round(best[1][1], 1)]; interp += 1
    print(f"  {len(solved)} solved + {interp} extrapolated = {len(pos)} positions; "
          f"{len(cat)-len(pos)} frames without one")

    nums = sorted(ptr_of)
    pairs = [[ptr_of[a], ptr_of[b]] for i, a in enumerate(nums) for b in nums[i + 1:]
             if b - a <= adj]
    have = [i['pointer'] for i in cat if os.path.exists(f"{SP}/fullres/{i['pointer']}.jpg")]
    print(f"  {len(have)} scans present in {SP}/fullres; {len(pairs)} adjacency pairs (|dnum| <= {adj})")
    json.dump(dict(recs=sorted(have), pos={k: v for k, v in pos.items() if int(k) in set(have)},
                   pairs=[p for p in pairs if p[0] in set(have) and p[1] in set(have)]),
              open(out, 'w'))
    print(f"wrote {out}")


if __name__ == '__main__':
    main()
