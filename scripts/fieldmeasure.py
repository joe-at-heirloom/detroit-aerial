#!/usr/bin/env python
"""Measure the absolute displacement field of a composited block, and save it.

The closed block has perfect seams and a large, smooth absolute error. Before
fitting anything to that error, look at it: how many parameters does it actually
have? Every failed placement so far fitted a model without first checking the
model was the right shape for the data.

Saves per-cell measurements (both grid scales) and per-frame observations to
data/field_<tag>_<label>.json for analysis.

Usage:  ./.venv/bin/python scripts/fieldmeasure.py 1961 closedA
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust, dtmap
import validate
from validate import P, MPP, report
from rebuild import placements, SP


def main():
    tag, label = sys.argv[1], sys.argv[2]
    ref = sys.argv[sys.argv.index('--ref') + 1] if '--ref' in sys.argv else 'modern'
    label_out = sys.argv[sys.argv.index('--label-out') + 1] if '--label-out' in sys.argv else label
    t0 = time.time()
    arr, mod, bbox = validate.load(tag, label, ref)
    t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_{label}_{ref.replace(":", "_")}', log=lambda *_: None)
    print(f"ridge maps ready [{time.time()-t0:.0f}s]", flush=True)
    # --coarse-search M : first-stage search radius (default 250). Against the
    # same film (e.g. a COLMAP block vs the old build of the same year) content is
    # identical, so a wide search is unambiguous, and a block seeded from catalogue
    # positions that scatter 361 m can sit outside +/-250.
    cs = float(sys.argv[sys.argv.index('--coarse-search') + 1]) if '--coarse-search' in sys.argv else 250.0
    pri = gridval.prior_for(t, MPP, search_m=cs)
    out = {'tag': tag, 'label': label, 'bbox': bbox, 'mpp': MPP,
           'shape': list(t['shape']),
           'minE': (bbox[1] - dtmap.LON0) * dtmap.MLON,
           'maxN': (bbox[2] - dtmap.LAT0) * dtmap.MLAT}
    for NY, NX in ((16, 3), (32, 6)):
        cells = gridval.grid_hier_prep(t, NY=NY, NX=NX, prior=pri)
        H, W = t['shape']; ch = H // NY; cw = W // NX
        for c in cells:
            c['cE'] = out['minE'] + (c['col'] + 0.5) * cw * MPP
            c['cN'] = out['maxN'] - (c['row'] + 0.5) * ch * MPP
        out[f'cells_{NY}x{NX}'] = cells
        report(cells, f'  {tag} {label} {NY}x{NX}')
    # per-frame observations too: more samples, and they say which flight line
    # which saved placement was this composite built from
    stage = {'closedA': 'stageA', 'closedA2': 'stageA2', 'closedA3': 'stageA3', 'closed': 'closed', 'placed': 'placed', 'placedW': 'stageA'}.get(label)
    if stage and not os.path.exists(P('data', f'{stage}_{tag}.json')):
        stage = None
    if stage:
        sol, ppm = placements(tag)
        saved = json.load(open(P('data', f'{stage}_{tag}.json')))
        for r in list(sol):
            v = saved.get(str(r)) or saved.get(r)
            if v:
                sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
                if 'rot' in v: sol[r]['rot'] = v['rot']
                if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
        abs_obs, _ = blockadjust.observe(sol, f"{SP}/fullres", t, out['minE'], out['maxN'],
                                         MPP, prior=pri, log=lambda *_: None)
        fr = []
        for r, o in abs_obs.items():
            fr.append(dict(rec=r, e=sol[r]['e'] - sol[r]['dE'], n=sol[r]['n'] - sol[r]['dN'],
                           dE=o['dE'], dN=o['dN'], ratio=o['ratio']))
        out['frames'] = fr
        v = np.array([math.hypot(f['dE'], f['dN']) for f in fr]) if fr else np.array([0.])
        print(f"  {len(fr)} frames matched: median {np.median(v):.1f}  p90 {np.percentile(v,90):.1f} m",
              flush=True)
    out['ref'] = ref
    json.dump(out, open(P('data', f'field_{tag}_{label_out}.json'), 'w'))
    print(f"saved data/field_{tag}_{label_out}.json (reference {ref}) [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
