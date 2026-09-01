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


def run(tag, stage='closed', label=None):
    """Composite a saved placement `stage` (data/<stage>_<tag>.json) to
    mosaics/detroit_<tag>_<label>.tif. Measurement is left to fieldmeasure.py,
    which does it with an alias-proof prior; the old zero-prior report here was
    what produced the 46 m number on a block that was 72 m off."""
    t0 = time.time()
    label = label or {'stageA': 'closedA', 'stageA2': 'closedA2'}.get(stage, stage)
    sol, ppm = placements(tag)
    saved = json.load(open(P('data', f'{stage}_{tag}.json')))
    n = 0
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
            if 'rot' in v: sol[r]['rot'] = v['rot']
            if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
            n += 1
    print(f"\n===== {tag}: compositing {stage} -> {label} ({n} frames) =====", flush=True)
    composite(tag, sol, ppm, label)
    print(f"  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    stage = sys.argv[sys.argv.index('--stage') + 1] if '--stage' in sys.argv else 'closed'
    label = sys.argv[sys.argv.index('--label') + 1] if '--label' in sys.argv else None
    for tg in [a for a in sys.argv[1:] if not a.startswith('--') and a not in (stage, label)]:
        run(tg, stage, label)
