#!/usr/bin/env python
"""Regenerate data/manifest.json from whatever mosaics actually exist.

Each block is served from the best build present -- `final` (per-frame corrected
plus residual field) if it is there, else `v2`, else the original `rbf` -- so the
viewer follows the pipeline forward without anything being edited by hand.
"""
import os, json, sys, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def P(*a): return os.path.join(ROOT, *a)

TAGS = ['1949', '1956', '1961', '1967']
# placed3 = closed by the per-frame tie-point bundle (close3), then placed by one
# continuous low-order warp with seams verified before and after. It is preferred
# wherever it exists; `final` is the earlier per-frame build whose seams were torn.
PREF = ['placedC3', 'placedC2', 'placedC', 'placed3', 'final', 'v2', 'rbf']


def best(tag, recorded=None):
    # --serve 1961:placedC,1949:placedC : the build to serve, per block, chosen
    # explicitly after its seams AND absolute placement verified. Nothing is
    # picked by label order -- that once served a failed placement because its
    # label sorted "newer".
    #
    # Without --serve, the choice recorded in the manifest by the last run is
    # kept, so a bare regeneration (to bump `ver`, to pick up a new layer) is
    # idempotent. It used to fall straight through to PREF, which silently
    # downgraded every block to its torn `final` build on 2026-09-04.
    chosen = dict(recorded or {})
    if '--serve' in sys.argv:
        for kv in sys.argv[sys.argv.index('--serve') + 1].split(','):
            k, v = kv.split(':'); chosen[k] = v
    allow = sys.argv[sys.argv.index('--placed3') + 1].split(',') if '--placed3' in sys.argv else []
    if tag in chosen:
        g = P('data', f'{tag}_{chosen[tag]}_geo.json'); m = P('mosaics', f'detroit_{tag}_{chosen[tag]}.tif')
        if os.path.exists(g) and os.path.exists(m):
            return chosen[tag], json.load(open(g))
        # The preference order below ends at `final`, the per-frame build whose
        # seams are torn, so falling into it is serving a known-bad mosaic. A
        # missing verified build therefore drops the block from the viewer unless
        # --fallback says otherwise; a hole is honest, a torn block is not.
        if '--fallback' not in sys.argv:
            print(f"  {tag}: chosen build {chosen[tag]} MISSING -- block left out (pass --fallback to serve the older builds)")
            return None, None
        print(f"  {tag}: chosen build {chosen[tag]} missing; falling back")
    for s in PREF:
        if s in ('placed3', 'placedC', 'placedC2', 'placedC3') and tag not in allow:
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
    for l in layers:
        fp = P(l['file'])
        if os.path.exists(fp):
            l['ver'] = int(os.path.getmtime(fp))
    ids = []
    boxes = []
    for tag in TAGS:
        s, g = best(tag, old.get('serve'))
        if not g:
            print(f"  {tag}: no mosaic"); continue
        lid = f'b{tag}'
        layers.append(dict(id=lid, group='west', label=tag,
                           file=f'mosaics/detroit_{tag}_{s}.tif',
                           bbox=g['bbox'], gray=True,
                           ver=int(os.path.getmtime(P('mosaics', f'detroit_{tag}_{s}.tif')))))
        ids.append(lid); boxes.append(g['bbox'])
        print(f"  {tag}: {s}")
    # Already-orthorectified layers fetched by fetchlayer.py -- the City of
    # Detroit's 1998 / 2005 / 2010 orthos and USDA NAIP 2012-2022 -- sit between
    # the film and Today, in year order. They are somebody else's georeferencing:
    # the film blocks were placed against modern imagery, these were not measured
    # against anything here, so a wipe between them and the film shows both.
    extra = []
    for g in sorted(glob.glob(P('data', 'layer_*_geo.json'))):
        name = os.path.basename(g)[len('layer_'):-len('_geo.json')]
        m = P('mosaics', f'layer_{name}.tif')
        if os.path.exists(m):
            extra.append((name, json.load(open(g))))
    for name, g in sorted(extra, key=lambda x: str(x[1].get('label', x[0]))):
        lid = f'l{name}'
        layers.append(dict(id=lid, group='west', label=str(g.get('label', name)),
                           file=f'mosaics/layer_{name}.tif', bbox=g['bbox'], gray=False,
                           ver=int(os.path.getmtime(P('mosaics', f'layer_{name}.tif')))))
        ids.append(lid)            # not in the union: these cover the blocks, never more
        print(f"  layer {name}: {g.get('source', {}).get('kind', '?')}")
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
                                'low-order warp against modern imagery. '
                                + ('1998, 2005 and 2010 are the City of Detroit\'s '
                                   'orthophotos and 2012-2022 are USDA NAIP, served '
                                   'as published; their georeferencing is theirs.'
                                   if extra else ''))}
    # keep a record for a block that was left out this run too, so the next bare
    # run still knows which build it is supposed to serve
    served = dict(old.get('serve') or {})
    served.update({l['label']: os.path.basename(l['file'])[len('detroit_') + len(l['label']) + 1:-len('.tif')]
                   for l in layers if l['id'].startswith('b')})
    json.dump(dict(layers=layers, groups=groups, serve=served),
              open(P('data', 'manifest.json'), 'w'), indent=1)
    print(f"  wrote manifest: {len(layers)} layers")


if __name__ == '__main__':
    main()
