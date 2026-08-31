"""Per-frame absolute adjustment.

A smooth warp field is the wrong shape for the error that is left. The mosaic is a
Voronoi composite of individual negatives, and each negative carries its own
residual position and crab angle; where two frames meet, the error can step
discontinuously. No continuous RBF can represent a step, so it splits the
difference and leaves roughly half the disagreement on both sides. That is what
keeps the worst cells at 30-40 m after the field has converged.

The frame is the physical unit, so correct the frame: render each negative alone
into map space, match it against modern imagery with the same alias-proof
coarse-to-fine search used everywhere else, and feed the result back as that
frame's own (dE, dN, rot) before re-compositing.
"""
import math, os, json
import numpy as np
from PIL import Image
import gridval
Image.MAX_IMAGE_PIXELS = None


def render_frame(rec, p, imgdir, minE, maxN, W, H, mpp, drot=0.0, crop=0.78):
    """One negative, alone, on the block's map grid. Returns (img, x0, y0) or None.

    Same inverse mapping mosaic_sim.build uses, minus the Voronoi selection, and
    keeping only the middle `crop` of the negative. The outer margin is where
    vignetting is worst and where relief displacement -- which no rigid placement
    can absorb -- grows with distance from the principal point, so including it
    blurs the correlation peak that decides this frame's position."""
    try:
        I = np.asarray(Image.open(f"{imgdir}/{rec}.jpg").convert('L'))
    except Exception:
        return None
    ih, iw = I.shape
    mpp_f = p['gw'] / iw
    e = p['e'] - p.get('dE', 0.0); n = p['n'] - p.get('dN', 0.0)
    rad = math.hypot(p['gw'], p['gh']) / 2 * (crop if crop < 1.0 else 1.0)
    x0 = max(0, int((e - rad - minE) / mpp)); x1 = min(W, int((e + rad - minE) / mpp))
    y0 = max(0, int((maxN - (n + rad)) / mpp)); y1 = min(H, int((maxN - (n - rad)) / mpp))
    if x1 - x0 < 32 or y1 - y0 < 32:
        return None
    th = math.radians(p['rot'] + drot); ct, st = math.cos(th), math.sin(th)
    xs = (minE + (np.arange(x0, x1) + 0.5) * mpp).astype(np.float32)
    ns = (maxN - (np.arange(y0, y1) + 0.5) * mpp).astype(np.float32)[:, None]
    dE = xs[None, :] - e; dN = ns - n
    lx = ct * dE + st * dN
    ly = -st * dE + ct * dN
    px = lx / mpp_f + iw / 2.0
    py = -ly / mpp_f + ih / 2.0
    m = (px >= 0) & (px < iw - 1) & (py >= 0) & (py < ih - 1)
    if crop < 1.0:
        m &= (np.abs(px - iw / 2.0) < crop * iw / 2.0) & \
             (np.abs(py - ih / 2.0) < crop * ih / 2.0)
    out = np.zeros(m.shape, np.uint8)
    if not m.any():
        return None
    xi = np.clip(px, 0, iw - 2); yi = np.clip(py, 0, ih - 2)
    xa = xi.astype(np.int32); ya = yi.astype(np.int32)
    fx = (xi - xa).astype(np.float32); fy = (yi - ya).astype(np.float32)
    v = (I[ya, xa] * (1 - fx) * (1 - fy) + I[ya, xa + 1] * fx * (1 - fy)
         + I[ya + 1, xa] * (1 - fx) * fy + I[ya + 1, xa + 1] * fx * fy)
    np.copyto(out, v.astype(np.uint8), where=m)
    return out, x0, y0


def solve_frames(sol, imgdir, t, minE, maxN, mpp, prior=None,
                 rots=(0.0,), min_valid=0.25, min_ratio=1.10, crop=0.78, log=print):
    """Per-frame (dE, dN, rot) against the modern ridge map already in `t`.

    `rots` is searched by re-rendering the frame at each angle -- rotating the
    negative smears the correlation peak, so the angle that sharpens it most is the
    frame's residual crab."""
    H, W = t['shape']
    out = {}
    for i, (rec, p) in enumerate(sorted(sol.items())):
        if not p.get('ok'):
            continue
        best = None
        for dr in rots:
            r = render_frame(rec, p, imgdir, minE, maxN, W, H, mpp, drot=dr, crop=crop)
            if r is None:
                continue
            img, x0, y0 = r
            if (img > 0).mean() < min_valid:
                continue
            rf, vf = gridval.ridge_full(img.astype(np.float32), mpp)
            y1 = y0 + img.shape[0]; x1 = x0 + img.shape[1]
            pr = prior.at((y0 + y1) / 2, (x0 + x1) / 2) if prior is not None else (0.0, 0.0)
            m = _match_local(rf, vf, t, y0, y1, x0, x1, mpp, pr)
            if m is None:
                continue
            m['rot'] = dr
            if best is None or m['ratio'] > best['ratio']:
                best = m
        if best is None or best['ratio'] < min_ratio:
            continue
        out[rec] = best
        if (i + 1) % 15 == 0:
            log(f"    {i+1}/{len(sol)} frames")
    return out


def _match_local(rf, vf, t, y0, y1, x0, x1, mpp, prior):
    """Correlate a frame-local ridge raster against the block's modern ridge map."""
    tf = {'mpp': mpp, 'shape': (y1 - y0, x1 - x0), 'ds': 1,
          'hf': rf, 'hv': vf, 'hc': rf, 'hvc': vf,
          'mf': gridval._cut(t['mf'], y0, y1, x0, x1),
          'mv': gridval._cut(t['mv'], y0, y1, x0, x1),
          'mc': gridval._cut(t['mf'], y0, y1, x0, x1),
          'mvc': gridval._cut(t['mv'], y0, y1, x0, x1)}
    return gridval.match_hier_prep(tf, 0, y1 - y0, 0, x1 - x0,
                                   fine_search=gridval.FINE_SEARCH,
                                   prior=prior, fine_iters=4)


def regularise(fx, sol, radius=4000.0, k=2.5, log=print):
    """Make the per-frame corrections cover every frame, consistently.

    Two failure modes to close. A frame that locks badly would be dragged away from
    its neighbours and open a visible step at the seam. And a frame that does not
    lock at all is worse: left uncorrected it stays where it was while the frames
    around it move by tens of metres, which tears the mosaic exactly where it used
    to be continuous. So outliers are replaced by their neighbourhood's consensus,
    and frames with no solution of their own borrow it."""
    import numpy as _np, math as _m
    ok = [r for r in sol if sol[r].get('ok')]
    locked = [r for r in ok if r in fx]
    if len(locked) < 6:
        log("    too few frames locked to regularise; leaving them alone")
        return {r: fx[r] for r in locked}
    P = _np.array([[sol[r]['e'], sol[r]['n']] for r in locked], float)
    V = _np.array([[fx[r]['dE'], fx[r]['dN']] for r in locked], float)
    W = _np.array([min(fx[r].get('ratio', 1.0), 3.0) for r in locked], float)

    def consensus(e, n, skip=None):
        d2 = ((P[:, 0] - e) ** 2 + (P[:, 1] - n) ** 2)
        w = W * _np.exp(-d2 / (2 * radius ** 2))
        if skip is not None:
            w[skip] = 0.0
        if w.sum() < 1e-9:
            j = int(_np.argmin(d2)); return V[j], 0.0
        loc = (V * w[:, None]).sum(0) / w.sum()
        spread = _m.sqrt(max((w * ((V - loc) ** 2).sum(1)).sum() / w.sum(), 25.0))
        return loc, spread

    out = {}; replaced = 0; borrowed = 0
    for i, r in enumerate(locked):
        loc, spread = consensus(P[i, 0], P[i, 1], skip=i)
        if spread > 0 and _m.hypot(*(V[i] - loc)) > k * spread:
            out[r] = dict(fx[r], dE=float(loc[0]), dN=float(loc[1]), borrowed=True)
            replaced += 1
        else:
            out[r] = dict(fx[r], borrowed=False)
    for r in ok:
        if r in out:
            continue
        loc, _ = consensus(sol[r]['e'], sol[r]['n'])
        out[r] = dict(dE=float(loc[0]), dN=float(loc[1]), rot=0.0, ratio=1.0,
                      borrowed=True)
        borrowed += 1
    log(f"    {replaced} corrections replaced by neighbour consensus, "
        f"{borrowed} frames had none of their own and borrowed one "
        f"({len(out)}/{len(ok)} frames now corrected)")
    return out
