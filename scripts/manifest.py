#!/usr/bin/env python
"""Regenerate data/manifest.json from whatever mosaics actually exist.

Each block is served from the best build present -- `final` (per-frame corrected
plus residual field) if it is there, else `v2`, else the original `rbf` -- so the
viewer follows the pipeline forward without anything being edited by hand.
"""
import os, json, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def P(*a): return os.path.join(ROOT, *a)

TAGS = ['1949', '1956', '1961', '1967']
# placed3 = closed by the per-frame tie-point bundle (close3), then placed by one
# continuous low-order warp with seams verified before and after. It is preferred
# wherever it exists; `final` is the earlier per-frame build whose seams were torn.
PREF = ['placedC', 'placed3', 'final', 'v2', 'rbf']


def best(tag):
    # --placed3 1961,1956 : only these tags may be served from placed3; a block
    # is switched when its seams AND absolute placement have been verified, not
    # when its file appears.
    allow = sys.argv[sys.argv.index('--placed3') + 1].split(',') if '--placed3' in sys.argv else []
    for s in PREF:
        if s in ('placed3', 'placedC') and tag not in allow:
            continue
        g = P('data', f'{tag}_{s}_geo.json')
        m = P('mosaics', f'detroit_{tag}_{s}.tif')
        if os.path.exists(g) and os.path.exists(m):
            return s, json.load(open(g))
    return None, None


def downtown_layers(old):
    """Serve the corrected downtown rasters where they exist.

    Downtown was 20-42 m out and nobody knew, because it sits east of the modern
    reference and every measurement silently returned "no peak". The corrected
    versions are GeoTIFFs, which serve.py reads the same way it reads the block
    mosaics."""
    out = []
    for l in old['layers']:
        if l['group'] != 'downtown':
            continue
        l = dict(l)
        fix = P('mosaics', f"downtown_{l['label']}_fix.tif")
        if l['id'] != 'dtnow' and os.path.exists(fix):
            l['file'] = f"mosaics/downtown_{l['label']}_fix.tif"
            g = P('data', f"dt{l['label']}_fix_geo.json")
            if os.path.exists(g):
                l['bbox'] = json.load(open(g))['bbox']
            print(f"  downtown {l['label']}: corrected")
        out.append(l)
    return out


def main():
    old = json.load(open(P('data', 'manifest.json')))
    layers = downtown_layers(old)
    ids = []
    boxes = []
    for tag in TAGS:
        s, g = best(tag)
        if not g:
            print(f"  {tag}: no mosaic"); continue
        lid = f'b{tag}'
        layers.append(dict(id=lid, group='west', label=tag,
                           file=f'mosaics/detroit_{tag}_{s}.tif',
                           bbox=g['bbox'], gray=True))
        ids.append(lid); boxes.append(g['bbox'])
        print(f"  {tag}: {s}")
    mw = json.load(open(P('data', 'modern_west_geo.json')))
    layers.append(dict(id='modwest', group='west', label='Today',
                       file='mosaics/modern_west.png', bbox=mw['bbox'], gray=True))
    ids.append('modwest'); boxes.append(mw['bbox'])
    union = [min(b[0] for b in boxes), min(b[1] for b in boxes),
             max(b[2] for b in boxes), max(b[3] for b in boxes)]
    groups = {'downtown': old['groups']['downtown'],
              'west': dict(label='West Detroit', layers=ids, bbox=union,
                           note='1949, 1956, 1961 and 1967 blocks. Each block closed '
                                'by a tie-point bundle adjustment -- a translation, '
                                'scale and rotation per negative from thousands of '
                                'phase-correlation windows, along-track and across '
                                'flight lines -- then placed by one continuous '
                                'low-order warp against modern imagery.')}
    json.dump(dict(layers=layers, groups=groups), open(P('data', 'manifest.json'), 'w'), indent=1)
    print(f"  wrote manifest: {len(layers)} layers")


if __name__ == '__main__':
    main()
