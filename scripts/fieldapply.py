#!/usr/bin/env python
"""Apply a fitted low-order field to the closed composite as a CONTINUOUS warp,
then verify both numbers -- seams and absolute -- on the result.

Why a warp of the composite and not per-frame transforms: the closed block has
0.19 m seams. A continuous deformation cannot open a seam; the disagreement
between two overlapping frames after warping is the disagreement before, scaled
by the field's local gradient (under 1%). Per-frame transforms can only reproduce
a field exactly when it is a similarity; anything with shear or curvature leaves
per-frame residue at every seam. So: composite first, warp second.

This is NOT the RBF warp the panel condemned. That was hundreds of parameters
fitted to aliasing control. This is 4-12 parameters fitted to alias-proof cell
measurements, chosen by held-out error, and it is judged afterwards by the seam
metric -- which touches no external reference -- and rejected if seams degrade.

Seam verification after a continuous warp is done honestly: every frame is
rendered alone into map space exactly as Stage A does, then each rendered frame
is pushed through the same field, and overlapping pairs are compared. That is
what happens to the composite.

Usage:  ./.venv/bin/python scripts/fieldapply.py 1961 closedA [--model affine] [--out placedW]
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import gridval, blockadjust, rbfapply, dtmap
import validate
from validate import P, MPP, report
from rebuild import placements, SP
import close as C
import fieldfit as FF


def make_warp(models):
    """Sum of several fitted increments. Fields here are under 1%, so the
    composition of two is their sum to second order (centimetres)."""
    if isinstance(models, dict):
        models = [models]
    def warp_at(E, N):
        dE = dN = 0.0
        for m in models:
            pe, pn = FF.predict_model(m, [E + dE], [N + dN])
            dE += float(pe[0]); dN += float(pn[0])
        return (dE, dN, 0.0)
    return warp_at


def warp_render(img, x0, y0, warp_at, minE, maxN, mpp, step_m=400.0):
    """Push one rendered frame through the field: out(p) = img(p + d(p))."""
    h, w = img.shape
    # field on a coarse lattice over this frame's extent, then bilinear
    nx = max(2, int(w * mpp / step_m) + 2); ny = max(2, int(h * mpp / step_m) + 2)
    gx = np.linspace(0, w, nx); gy = np.linspace(0, h, ny)
    DE = np.zeros((ny, nx), np.float32); DN = np.zeros((ny, nx), np.float32)
    for j, oy in enumerate(gy):
        for i, ox in enumerate(gx):
            E = minE + (x0 + ox) * mpp; N = maxN - (y0 + oy) * mpp
            d = warp_at(E, N); DE[j, i] = d[0]; DN[j, i] = d[1]
    xs = np.arange(w) + 0.5; ys = np.arange(h) + 0.5
    de = rbfapply._bilerp(DE, gx, gy, xs, ys); dn = rbfapply._bilerp(DN, gx, gy, xs, ys)
    sx = xs[None, :] + de / mpp - 0.5
    sy = ys[:, None] - dn / mpp - 0.5
    inb = (sx >= 0) & (sx < w - 1) & (sy >= 0) & (sy < h - 1)
    xi = np.clip(sx, 0, w - 1.001); yi = np.clip(sy, 0, h - 1.001)
    xa = xi.astype(np.int32); ya = yi.astype(np.int32)
    fx = (xi - xa).astype(np.float32); fy = (yi - ya).astype(np.float32)
    v = (img[ya, xa] * (1 - fx) * (1 - fy) + img[ya, xa + 1] * fx * (1 - fy)
         + img[ya + 1, xa] * (1 - fx) * fy + img[ya + 1, xa + 1] * fx * fy)
    return np.where(inb, v, 0).astype(np.uint8)


def seams_after_warp(tag, stage, warp_at, log=print):
    sol, ppm = placements(tag)
    saved = json.load(open(P('data', f'{stage}_{tag}.json')))
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r]['dE'] = v['dE']; sol[r]['dN'] = v['dN']
            if 'rot' in v: sol[r]['rot'] = v['rot']
            if 'gw' in v: sol[r]['gw'] = v['gw']; sol[r]['gh'] = v['gh']
    mpp = C.MPP_TIE
    minE, maxN, W, H = C.canvas(sol, mpp)
    rend = blockadjust.render_all(sol, f"{SP}/fullres", minE, maxN, W, H, mpp,
                                  crop=0.95, log=lambda *_: None)
    # Along-track and cross-line pairs are reported separately. A block of parallel
    # lines can be closed along each line and still have the lines disagree by
    # 100 m; one median over all pairs hides that because along-track pairs
    # outnumber sidelap pairs three to one. Cross-line is the number that says
    # whether a per-line correction did its job.
    import seamclass as SC
    recs = sorted(rend)
    lab = SC.lines_of(sol, recs)
    def stats(rd, nm):
        obs = blockadjust.relative_observations(rd, mpp, search_m=150.0, min_ratio=1.08,
                                                log=lambda *_: None)
        cls = {'along': [], 'cross': []}
        for a, b, de, dn, rt in obs:
            cls['along' if lab[a] == lab[b] else 'cross'].append(math.hypot(de, dn))
        out = {}
        for k in ('along', 'cross'):
            v = np.array(cls[k])
            if len(v):
                log(f"    {nm} {k:5s}: {len(v):3d} pairs  median {np.median(v):5.1f}  "
                    f"p90 {np.percentile(v,90):5.1f}  max {v.max():5.1f} m  over 10 m {(v>10).sum()}")
                out[k] = float(np.median(v))
            else:
                log(f"    {nm} {k:5s}:   0 pairs"); out[k] = None
        return out
    s0 = stats(rend, 'before warp')
    wr = {r: (warp_render(img, x0, y0, warp_at, minE, maxN, mpp), x0, y0)
          for r, (img, x0, y0) in rend.items()}
    del rend
    s1 = stats(wr, 'after  warp')
    del wr
    return s0, s1


def main():
    tag, label = sys.argv[1], sys.argv[2]
    kind = sys.argv[sys.argv.index('--model') + 1] if '--model' in sys.argv else None
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else 'placedW'
    # --fits a,b,c : fieldfit labels whose chosen/--model increments are summed,
    # e.g. --fits closedA,placedW  applies fit(closedA) + fit(placedW) to closedA
    fits = sys.argv[sys.argv.index('--fits') + 1].split(',') if '--fits' in sys.argv else [label]
    t0 = time.time()
    models = []
    for fl in fits:
        FIT = json.load(open(P('data', f'fieldfit_{tag}_{fl}.json')))
        k = kind or FIT['choice']
        models.append(FIT['models'][k])
    model = models[-1]; kind = model['kind']
    warp_at = make_warp(models)
    npar = sum(4 * len(m['lines']) if m['kind'] == 'perline' else len(m['p']) for m in models)
    print(f"\n===== {tag}: warp {label} with {'+'.join(fits)} {kind} field ({npar} params, "
          f"last held-out {model['loo_median']:.1f} m) =====", flush=True)
    json.dump(models, open(P('data', f'fieldmodel_{tag}_{out}.json'), 'w'))

    stage = {'closedA': 'stageA', 'closed': 'closed', 'placed': 'placed'}[label]
    print("  internal consistency, before and after the same warp:", flush=True)
    s0, s1 = seams_after_warp(tag, stage, warp_at)
    a0, a1 = s0.get('along'), s1.get('along')
    if a0 is not None and a1 is not None and a1 > a0 + 1.5:
        print(f"  the warp degrades along-track seams ({a0:.2f} -> {a1:.2f} m); refusing to apply it")
        return

    geo = json.load(open(P('data', f'{tag}_{label}_geo.json')))
    g = dict(minE=geo['minE'], maxN=geo['maxN'], W=geo['W'], H=geo['H'], mpp=geo['mpp'])
    res = rbfapply.apply(P('mosaics', f'detroit_{tag}_{label}.tif'), g, warp_at,
                         P('mosaics', f'detroit_{tag}_{out}.tif'), f'/tmp/fw_{tag}.raw',
                         step_m=250.0, log=lambda s: None)
    json.dump(res, open(P('data', f'{tag}_{out}_geo.json'), 'w'))
    print(f"  wrote mosaics/detroit_{tag}_{out}.tif [{time.time()-t0:.0f}s]", flush=True)

    print("  absolute placement against modern imagery:", flush=True)
    arr, mod, bbox = validate.load(tag, out)
    t = gridval.prepare_cached(arr, mod, MPP, f'{tag}_{out}_modern', log=lambda *_: None)
    pri = gridval.prior_for(t, MPP, log=lambda *_: None)
    a = gridval.grid_hier_prep(t, prior=pri); report(a, f'    {tag} 3.6 x 2.0 km')
    b = gridval.grid_hier_prep(t, NY=32, NX=6, prior=pri); report(b, f'    {tag} 1.8 x 1.0 km')
    json.dump(dict(coarse=a, fine=b), open(P('data', f'gridval_{tag}_{out}.json'), 'w'))
    print(f"  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
