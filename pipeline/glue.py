"""Learned feature matching for frames that intensity correlation cannot place.

About 14 of 51 frames in the 1967 block and 16 of 50 in 1949 never lock against
modern imagery. That is not a search problem -- widening the search only makes it
ambiguous -- it is the appearance gap between mid-century panchromatic film and
present-day satellite imagery. Ridge correlation compares brightness structure;
where the neighbourhood has been rebuilt, or the film is hazy, there is no
brightness structure in common to correlate.

Learned local features are trained for exactly that gap, so this uses LightGlue
(over DISK or ALIKED keypoints) via kornia.

**It does not replace the alias-proof search.** A keypoint matcher is not immune to
Detroit's 97.5 m street lattice either -- repeated structure produces repeated
descriptors. So LightGlue only ever supplies a *candidate* displacement for a frame
that could not otherwise be placed; the existing +/-40 m fine refinement still
decides the final number, and still rejects a candidate that is a block off.

MEASURED RESULT, 1967 block, DISK + LightGlue at side=1024 on 2.5 m/px renders:

    ridge correlation locked   39 / 51 frames
    LightGlue locked           13 / 51 frames
    both locked                 9, agreeing within 20 m in only 3
    frames rescued by LightGlue 4

So on this data, as configured, it is markedly *worse* than the ridge correlation
it was meant to rescue, and it does not agree with it where both fire. It verifies
perfectly on synthetic pairs (five planted shifts up to 120 m recovered exactly,
~1700 inliers), so the implementation is sound -- the failure is distributional.
DISK is trained on close-range natural scenes; nadir aerial imagery of a repeating
street grid, resampled to metres per pixel, is not that. Feeding it at 2.5 m/px
also leaves a house about three pixels across, which is very little for a keypoint
detector.

Worth trying before writing the approach off: render frames nearer 1 m/px, raise
`side`, swap DISK for ALIKED, or feed the ridge response instead of raw film so the
two epochs look alike. Keep it as an assist for frames that fail, not a
replacement -- and keep measuring, because the honest comparison above is the only
reason to prefer one matcher over the other.
"""
import math
import numpy as np

_MODEL = None


def available():
    try:
        import torch, kornia  # noqa: F401
        return True
    except Exception:
        return False


def _load(device=None):
    """Detector + matcher, built once. Falls back through kornia's options so a
    version that names things differently still works."""
    global _MODEL
    if _MODEL is not None:
        return _MODEL
    import torch, kornia
    from kornia.feature import LightGlueMatcher, laf_from_center_scale_ori
    dev = torch.device(device or ('cuda' if torch.cuda.is_available()
                                  else 'mps' if torch.backends.mps.is_available()
                                  else 'cpu'))
    det = None; name = None
    for cand, nm in ((lambda: kornia.feature.DISK.from_pretrained('depth'), 'disk'),
                     (lambda: kornia.feature.KeyNetAffNetHardNet(upright=True), 'sift')):
        try:
            det = cand().to(dev).eval(); name = nm; break
        except Exception:
            continue
    if det is None:
        raise RuntimeError('no kornia detector available')
    m = LightGlueMatcher(name).to(dev).eval()
    _MODEL = (det, m, name, dev, torch, kornia, laf_from_center_scale_ori)
    return _MODEL


def _prep(a, side, torch, dev):
    """Contrast-normalise and scale to `side` px. Nodata stays black; the network
    sees it as texture-free and finds nothing there, which is what we want."""
    import torch.nn.functional as F
    v = a > 0
    x = a.astype(np.float32)
    if v.sum() > 100:
        lo, hi = np.percentile(x[v], [2, 98])
        x = np.clip((x - lo) / max(float(hi - lo), 1e-6), 0, 1).astype(np.float32)
    x[~v] = 0.0
    t = torch.from_numpy(np.ascontiguousarray(x, dtype=np.float32))[None, None].to(dev)
    sc = side / max(a.shape)
    if sc < 1.0:
        t = F.interpolate(t, scale_factor=sc, mode='bilinear', align_corners=False)
    else:
        sc = 1.0
    return t, sc


def match(hist, modern, mpp, side=1024, min_inliers=15, tol_m=25.0, device=None):
    """Displacement of `hist` content relative to `modern`, both on the same grid.

    Returns dict(dE, dN, inliers, spread, ratio) or None. Sign convention matches
    gridval: positive dE means the historical content sits that far east of where
    the modern imagery puts it."""
    det, matcher, name, dev, torch, kornia, laf_from = _load(device)
    with torch.inference_mode():
        ta, sa = _prep(hist, side, torch, dev)
        tb, sb = _prep(modern, side, torch, dev)
        if name == 'disk':
            # DISK is an RGB network; the negatives are panchromatic, so the single
            # channel is repeated rather than colourised
            fa = det(ta.repeat(1, 3, 1, 1), 2048, pad_if_not_divisible=True)[0]
            fb = det(tb.repeat(1, 3, 1, 1), 2048, pad_if_not_divisible=True)[0]
            ka, da = fa.keypoints, fa.descriptors
            kb, db = fb.keypoints, fb.descriptors
        else:
            la, ra, da = det(ta)
            lb, rb, db = det(tb)
            ka = la[0, :, :, 2]; kb = lb[0, :, :, 2]
            da = da[0]; db = db[0]
        if len(ka) < min_inliers or len(kb) < min_inliers:
            return None
        lafa = laf_from(ka[None], torch.ones(1, len(ka), 1, 1, device=dev) * 8.0)
        lafb = laf_from(kb[None], torch.ones(1, len(kb), 1, 1, device=dev) * 8.0)
        hw1 = (ta.shape[-2], ta.shape[-1]); hw2 = (tb.shape[-2], tb.shape[-1])
        _, idx = matcher(da, db, lafa, lafb, hw1, hw2)
        if idx is None or len(idx) < min_inliers:
            return None
        pa = ka[idx[:, 0]].cpu().numpy() / sa
        pb = kb[idx[:, 1]].cpu().numpy() / sb

    d = pa - pb                              # historical minus modern, in pixels
    if len(d) < min_inliers:
        return None
    # robust translation: the matcher will still produce some wrong pairs, and on a
    # repeating street grid the wrong ones cluster one block away, so take the
    # densest agreement rather than a mean
    c = np.median(d, axis=0)
    for _ in range(6):
        r = np.hypot(d[:, 0] - c[0], d[:, 1] - c[1])
        s = max(np.median(r) * 1.4826, tol_m / mpp * 0.5)
        w = 1.0 / (1.0 + (r / (1.5 * s)) ** 2)
        c = (d * w[:, None]).sum(0) / w.sum()
    r = np.hypot(d[:, 0] - c[0], d[:, 1] - c[1])
    inl = int((r <= tol_m / mpp).sum())
    if inl < min_inliers:
        return None
    dE = float(c[0] * mpp); dN = float(-c[1] * mpp)
    return dict(dE=dE, dN=dN, inliers=inl, total=int(len(d)),
                spread=float(np.median(r) * mpp),
                ratio=float(inl / max(len(d) - inl, 1)))
