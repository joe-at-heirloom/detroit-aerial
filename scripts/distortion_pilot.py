#!/usr/bin/env python
"""Is there lens or scan distortion the fixed camera leaves behind? (plan step 3)

The COLMAP builds use one fixed SIMPLE_RADIAL camera (focal 3602 px, k -0.00027)
because self-calibration on flat ground is degenerate. If the real lens bends rays
differently by dk, every negative's content lands displaced on the ground by

    d(rho) = dk * |rho|^2 * rho / H^2          (+ dk2 * |rho|^4 * rho / H^4)

where rho is the ground vector from that negative's nadir and H the flying height.
A whole-map warp cannot remove that: it repeats in every frame, about every frame's
own centre. Two tests, one without any reference:

  seams     every tie window in every overlap sees one ground point at two
            different rho. The pair's disagreement there is d(rho_a) - d(rho_b),
            which VARIES across the overlap; a per-frame shift, scale or rotation
            is constant across it and is carried as a per-pair nuisance. dk is
            fitted on 4/5 of the pairs and scored on the fifth, on the within-pair
            scatter it removes -- so it has to generalise to earn anything.
  absolute  each rendered negative against the reference in ~500 m windows, with a
            per-frame similarity removed (placement and the block's smooth bow),
            and what remains stacked in IMAGE coordinates about the principal
            point: a radial pattern shared by every frame is the camera; a pattern
            that differs by scan batch is the scan; a constant is placement.

Scan affinity per batch is also tested on the seams: frames from different scan
windows could carry different x/y scale, which only shows where batches meet.

Usage:  ./.venv/bin/python scripts/distortion_pilot.py runs/colmap/1961 --tag 1961
            [--ref l2016] [--surface 2] [--mpp 2.0] [--out runs/pilot_1961]
"""
import sys, os, json, math, time, argparse
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image
import dtmap, gridval
import colmap_render as CR
import close3, seamclass as SC
import crossmatrix as XM
Image.MAX_IMAGE_PIXELS = None


def subpixel_pc(A, B, pad=2):
    """close2._pc with a parabolic sub-pixel peak. The stock correlator returns whole
    pixels, so at 2 m/px every seam window reads 0, 2, 2.8, 4... m and a median of
    2.0 means 'within one pixel' -- too coarse to see a metre of lens error, and
    quantised enough that any smooth model can appear to help."""
    import close2
    h, w = A.shape
    wn = np.outer(np.hanning(h), np.hanning(w)).astype(np.float32)
    a = (A - A.mean()) * wn; b = (B - B.mean()) * wn
    H, W = h * pad, w * pad
    F = np.fft.rfft2(a, s=(H, W)) * np.conj(np.fft.rfft2(b, s=(H, W)))
    F /= (np.abs(F) + 1e-9)
    c = np.fft.fftshift(np.fft.irfft2(F, s=(H, W)))
    pk = np.unravel_index(np.argmax(c), c.shape)
    peak = float(c[pk]); sharp = peak / (float(c.std()) + 1e-12)
    oy, ox = gridval._parabolic(c, pk)
    dy = pk[0] - H // 2 + oy; dx = pk[1] - W // 2 + ox
    at_edge = abs(dx) > w * 0.60 or abs(dy) > h * 0.60
    return dx, dy, sharp, at_edge


def render(work, surface, mpp):
    meta, cams, imgs, zg = CR.load_block(work, surface=surface)
    E0, N0 = meta['E0'], meta['N0']
    C = np.array([imgs[n]['C'] for n in imgs]); half = 2600.0
    minE = C[:, 0].min() + E0 - half; maxE = C[:, 0].max() + E0 + half
    minN = C[:, 1].min() + N0 - half; maxN = C[:, 1].max() + N0 + half
    W = int((maxE - minE) / mpp); H = int((maxN - minN) / mpp)
    rend, fr = {}, {}
    for name, img in imgs.items():
        rec = name.replace('.jpg', '')
        im = np.asarray(Image.open(f"{work}/images/{name}").convert('L'))
        cam = cams[img['cam']]
        r = CR.render_frame(name, img, cam, im, zg, E0, N0, minE, maxN, W, H, mpp)
        if r is None:
            continue
        rend[rec] = r
        C_ = img['C']
        z = float(zg.z(C_[0], C_[1])) if isinstance(zg, CR.Surface) else float(zg)
        # image x/y axes on the ground, for stacking residuals in image coordinates
        ax = img['R'].T @ np.array([1.0, 0, 0]); ay = img['R'].T @ np.array([0, 1.0, 0])
        # nadir: where the camera's plumb line meets the ground -- approximated by
        # the centre's E/N, since tilts are a degree or two over flat ground
        fr[rec] = dict(E=C_[0] + E0, N=C_[1] + N0, H=float(C_[2] - z),
                       ax=ax[:2] / (np.linalg.norm(ax[:2]) + 1e-12),
                       ay=ay[:2] / (np.linalg.norm(ay[:2]) + 1e-12),
                       size=Image.open(f"{work}/images/{name}").size)
    grid = dict(minE=minE, maxN=maxN, W=W, H=H, mpp=mpp)
    return rend, fr, grid, cams


def cubic(rho, H, p=2):
    n2 = (rho ** 2).sum(-1, keepdims=True)
    return rho * n2 ** (p // 2) / H ** p


def aniso(rho, f):
    """x/y differential scale and shear in the frame's image axes, back on the ground."""
    x = rho @ f['ax']; y = rho @ f['ay']
    return (np.outer(x, f['ax']) - np.outer(y, f['ay']),      # stretch x, squeeze y
            np.outer(y, f['ax']) + np.outer(x, f['ay']))      # shear


def seam_design(ties, fr, grid, batch, ref_batch):
    """Rows (2 per tie) of [cubic, quintic, per-batch aniso..., per-pair translation...]."""
    pairs = sorted({(t['a'], t['b']) for t in ties}); pid = {p: i for i, p in enumerate(pairs)}
    bats = sorted(set(batch.values()) - {ref_batch})
    nb = len(bats) * 2
    NG = 3                                   # cubic, quintic, null cubic
    # the control: the same cubic about a deliberately wrong centre per frame
    # (1.5 km off, random direction). Real lens geometry must beat it; if it helps
    # as much, the 'improvement' is the model soaking up measurement structure.
    rng = np.random.default_rng(7)
    off = {r: 1500.0 * np.array([math.cos(t), math.sin(t)])
           for r, t in zip(sorted(fr), rng.uniform(0, 2 * math.pi, len(fr)))}
    X = np.zeros((2 * len(ties), NG + nb + 2 * len(pairs))); y = np.zeros(2 * len(ties))
    g = np.zeros(2 * len(ties), int)
    for i, t in enumerate(ties):
        E = grid['minE'] + t['qx'] * grid['mpp']; N = grid['maxN'] - t['qy'] * grid['mpp']
        fa, fb = fr[t['a']], fr[t['b']]
        ra = np.array([[E - fa['E'], N - fa['N']]]); rb = np.array([[E - fb['E'], N - fb['N']]])
        c3 = cubic(ra, fa['H']) - cubic(rb, fb['H'])
        c5 = cubic(ra, fa['H'], 4) - cubic(rb, fb['H'], 4)
        X[2*i:2*i+2, 0] = c3[0]; X[2*i:2*i+2, 1] = c5[0]
        cn = cubic(ra + off[t['a']], fa['H']) - cubic(rb + off[t['b']], fb['H'])
        X[2*i:2*i+2, 2] = cn[0]
        for side, f, rr, sgn in (('a', fa, ra, 1), ('b', fb, rb, -1)):
            bt = batch[t[side]]
            if bt in bats:
                k = NG + 2 * bats.index(bt)
                s1, s2 = aniso(rr, f)
                X[2*i:2*i+2, k] += sgn * s1[0]; X[2*i:2*i+2, k+1] += sgn * s2[0]
        j = NG + nb + 2 * pid[(t['a'], t['b'])]
        X[2*i, j] = 1; X[2*i+1, j+1] = 1
        y[2*i] = t['dE']; y[2*i+1] = t['dN']; g[2*i:2*i+2] = pid[(t['a'], t['b'])]
    return X, y, g, pairs, bats


def robust_lstsq(X, y, iters=5):
    w = np.ones(len(y))
    for _ in range(iters):
        b = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)[0]
        r = y - X @ b; s = max(np.median(np.abs(r)) * 1.4826, 0.3)
        w = 1.0 / np.sqrt(1.0 + (r / (2 * s)) ** 2)
    return b


def cv_seams(X, y, g, cols, folds=5, seed=0):
    """Held-out within-pair scatter. The model's global terms (cols) come from the
    training pairs; on a held-out pair only its own translation is re-fitted, so all
    the model can do there is explain how the disagreement VARIES across the overlap."""
    npair = g.max() + 1
    rng = np.random.default_rng(seed); fold = rng.integers(0, folds, npair)
    nglob = X.shape[1] - 2 * npair
    out = []
    for k in range(folds):
        tr = fold[g] != k; te = ~tr
        keep = [c for c in range(nglob) if c in cols] + list(range(nglob, X.shape[1]))
        b = np.zeros(X.shape[1])
        if tr.any():
            bb = robust_lstsq(X[tr][:, keep], y[tr]); b[keep] = bb
        for p in np.unique(g[te]):
            m = g == p
            r = y[m] - X[m][:, :nglob] @ b[:nglob]
            r2 = r.reshape(-1, 2); r2 = r2 - np.median(r2, axis=0)
            out.extend(np.hypot(r2[:, 0], r2[:, 1]))
    return np.array(out)


def absolute(rend, fr, grid, ref, win_m=520.0):
    """Per frame: coarse offset against the reference, then windows; per-frame
    similarity removed; residuals returned in image axes about the frame centre."""
    s = dtmap.LAT0 + (grid['maxN'] - grid['H'] * grid['mpp']) / dtmap.MLAT
    n = dtmap.LAT0 + grid['maxN'] / dtmap.MLAT
    w = dtmap.LON0 + grid['minE'] / dtmap.MLON
    e = dtmap.LON0 + (grid['minE'] + grid['W'] * grid['mpp']) / dtmap.MLON
    g = dict(bbox=[s, w, n, e], W=grid['W'], H=grid['H'], mpp=grid['mpp'])
    L = {l['id']: l for l in XM.load_layers('west')}
    t0 = time.time()
    R = XM.read_onto(L[ref], g, None)
    print(f"  reference {ref} on the render grid, coverage {(R>0).mean():.2f} [{time.time()-t0:.0f}s]", flush=True)
    mpp = grid['mpp']
    Rf, Rv = gridval.ridge_full(R.astype(np.float32), mpp)
    Rc, Rvc = gridval.ridge_full(R.astype(np.float32), mpp, gridval.COARSE_BAND, gridval.COARSE_FLAT)
    wp = int(win_m / mpp)
    rows = []
    for rec, (img, x0, y0) in sorted(rend.items()):
        h, w_ = img.shape
        A = img.astype(np.float32)
        Af, Av = gridval.ridge_full(A, mpp)
        Ac, Avc = gridval.ridge_full(A, mpp, gridval.COARSE_BAND, gridval.COARSE_FLAT)
        sl = (slice(y0, y0 + h), slice(x0, x0 + w_))
        c = gridval.match_cell(Ac, Avc.astype(np.float32), Rc[sl], Rvc[sl].astype(np.float32), mpp, 300.0)
        if c is None or c['edge'] or c['ratio'] < 1.15:
            continue
        oy = int(round(c['dN'] / mpp)); ox = int(round(c['dE'] / mpp))
        f = fr[rec]
        for yy in range(0, h - wp + 1, wp // 2):
            for xx in range(0, w_ - wp + 1, wp // 2):
                a_ = Af[yy:yy + wp, xx:xx + wp]; ma = Av[yy:yy + wp, xx:xx + wp].astype(np.float32)
                if ma.mean() < 0.9:
                    continue
                Y0, X0 = y0 + yy + oy, x0 + xx - ox       # reference window under the coarse offset
                b_ = gridval._cut(Rf, Y0, Y0 + wp, X0, X0 + wp)
                mb = gridval._cut(Rv, Y0, Y0 + wp, X0, X0 + wp).astype(np.float32)
                if mb.mean() < 0.9 or a_.std() < 1e-6 or b_.std() < 1e-6:
                    continue
                r = gridval.match_cell(a_, ma, b_, mb, mpp, 40.0)
                if r is None or r['edge'] or r['ratio'] < 1.2:
                    continue
                E = grid['minE'] + (x0 + xx + wp / 2) * mpp; N = grid['maxN'] - (y0 + yy + wp / 2) * mpp
                rows.append(dict(rec=rec, rE=E - f['E'], rN=N - f['N'],
                                 dE=c['dE'] + r['dE'], dN=c['dN'] + r['dN'], ratio=r['ratio']))
    return rows


def absolute_fit(rows, fr):
    """Per-frame similarity (4 params) removed; the shared cubic fitted on top.
    Returns residuals in image axes and the fitted dk."""
    recs = sorted({r['rec'] for r in rows}); rid = {r: i for i, r in enumerate(recs)}
    X = np.zeros((2 * len(rows), 1 + 4 * len(recs))); y = np.zeros(2 * len(rows))
    for i, r in enumerate(rows):
        f = fr[r['rec']]; rho = np.array([[r['rE'], r['rN']]])
        X[2*i:2*i+2, 0] = cubic(rho, f['H'])[0]
        j = 1 + 4 * rid[r['rec']]
        X[2*i, j] = 1; X[2*i+1, j+1] = 1
        X[2*i, j+2] = r['rE']; X[2*i+1, j+2] = r['rN']          # scale
        X[2*i, j+3] = -r['rN']; X[2*i+1, j+3] = r['rE']         # rotation
        y[2*i] = r['dE']; y[2*i+1] = r['dN']
    b0 = robust_lstsq(X[:, 1:], y)
    res0 = y - X[:, 1:] @ b0
    b1 = robust_lstsq(X, y)
    res1 = y - X @ b1
    out = []
    for i, r in enumerate(rows):
        f = fr[r['rec']]
        v = res0[2*i:2*i+2]; rho = np.array([r['rE'], r['rN']])
        out.append(dict(rec=r['rec'], x=float(rho @ f['ax']), y=float(rho @ f['ay']),
                        vx=float(v @ f['ax']), vy=float(v @ f['ay']),
                        res_sim=float(np.hypot(*v)), res_cubic=float(np.hypot(*res1[2*i:2*i+2]))))
    return out, float(b1[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('work'); ap.add_argument('--tag', default='1961')
    ap.add_argument('--ref', default='l2016'); ap.add_argument('--surface', type=int, default=2)
    ap.add_argument('--mpp', type=float, default=2.0); ap.add_argument('--out')
    ap.add_argument('--no-absolute', action='store_true')
    ap.add_argument('--subpixel', action='store_true', help='parabolic peaks in the seam ties')
    a = ap.parse_args()
    out = a.out or os.path.join(ROOT, 'runs', f'pilot_{a.tag}')
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    rend, fr, grid, cams = render(a.work, a.surface, a.mpp)
    cam = cams[list(cams)[0]]
    Hs = np.array([f['H'] for f in fr.values()])
    print(f"{a.tag}: {len(rend)} frames rendered, camera {cam['model']} {np.round(cam['params'], 5).tolist()}, "
          f"flying height {np.median(Hs):.0f} m (range {Hs.min():.0f}-{Hs.max():.0f}) [{time.time()-t0:.0f}s]", flush=True)
    recs = sorted(rend)
    lab = SC.lines_of({r: dict(e=fr[r]['E'], n=fr[r]['N']) for r in recs}, recs)
    if a.subpixel:
        import close2
        close2._pc = subpixel_pc
    ties = close3.tie_windows(rend, lab, a.mpp, log=print)
    close3.seam_stats(ties, 'seams')
    # scan batch = the crop's width before padding (colmap_block pads to one canvas;
    # the scan window, not the negative, sets the width)
    sizes = {r: Image.open(os.path.join(ROOT, 'scans', 'fullres', f'{r}.jpg')).size for r in recs}
    batch = {r: f"{sizes[r][0]}x{sizes[r][1]}" for r in recs}
    from collections import Counter
    bc = Counter(batch.values()); ref_batch = bc.most_common(1)[0][0]
    print(f"  scan batches: {dict(bc)}", flush=True)
    X, y, g, pairs, bats = seam_design(ties, fr, grid, batch, ref_batch)
    nb = 2 * len(bats)
    NG = 3
    models = {'none (fixed camera as solved)': [],
              'CONTROL: r^3 about wrong centres': [2],
              'radial r^3': [0],
              'radial r^3 + r^5': [0, 1],
              'scan affinity per batch': list(range(NG, NG + nb)),
              'radial r^3 + scan affinity': [0] + list(range(NG, NG + nb))}
    rep = dict(tag=a.tag, frames=len(rend), ties=len(ties), pairs=len(pairs), batches=dict(bc),
               camera=dict(model=cam['model'], params=list(map(float, cam['params']))),
               H=float(np.median(Hs)), subpixel=a.subpixel, seams={}, when=time.strftime('%Y-%m-%d %H:%M'))
    print("\n  SEAMS: held-out within-pair scatter (5-fold by pair; lower is better)")
    base = None
    for nm, cols in models.items():
        r = cv_seams(X, y, g, cols)
        med, p90 = float(np.median(r)), float(np.percentile(r, 90))
        base = base or med
        full = robust_lstsq(X[:, cols + list(range(NG + nb, X.shape[1]))], y) if cols else None
        coef = {('dk' if c == 0 else 'dk2' if c == 1 else 'null' if c == 2 else f'aniso{c-NG}'): float(full[i])
                for i, c in enumerate(cols)} if cols else {}
        rep['seams'][nm] = dict(median=med, p90=p90, coef=coef)
        cs = '  '.join(f"{k} {v:+.2e}" for k, v in coef.items() if not k.startswith('aniso'))
        print(f"    {nm:32s} median {med:5.2f}  p90 {p90:5.2f} m   {cs}", flush=True)
    # effect size: what the fitted radial term would move at the frame corner
    dk = rep['seams']['radial r^3']['coef']['dk']
    corner = math.hypot(*[s / 2 for s in fr[recs[0]]['size']]) / cam['params'][0] * np.median(Hs)
    rep['radial_at_corner_m'] = float(abs(dk) * corner ** 3 / np.median(Hs) ** 2)
    print(f"    fitted dk would move a frame corner ({corner:.0f} m out) by {rep['radial_at_corner_m']:.2f} m", flush=True)

    if not a.no_absolute:
        print(f"\n  ABSOLUTE against {a.ref}, per-frame similarity removed", flush=True)
        rows = absolute(rend, fr, grid, a.ref)
        res, dk_abs = absolute_fit(rows, fr)
        v0 = np.array([r['res_sim'] for r in res]); v1 = np.array([r['res_cubic'] for r in res])
        print(f"    {len(rows)} windows on {len({r['rec'] for r in rows})} frames; residual after similarity "
              f"median {np.median(v0):.2f} p90 {np.percentile(v0,90):.2f}; with shared cubic "
              f"{np.median(v1):.2f} / {np.percentile(v1,90):.2f}; dk {dk_abs:+.2e}", flush=True)
        rep['absolute'] = dict(ref=a.ref, windows=len(rows), sim_median=float(np.median(v0)),
                               sim_p90=float(np.percentile(v0, 90)), cubic_median=float(np.median(v1)),
                               cubic_p90=float(np.percentile(v1, 90)), dk=dk_abs)
        # radial profile in image coordinates, all frames stacked
        rr = np.array([math.hypot(r['x'], r['y']) for r in res])
        rad = np.array([(r['vx'] * r['x'] + r['vy'] * r['y']) / max(math.hypot(r['x'], r['y']), 1) for r in res])
        tan = np.array([(-r['vx'] * r['y'] + r['vy'] * r['x']) / max(math.hypot(r['x'], r['y']), 1) for r in res])
        bins = np.arange(0, rr.max() + 300, 300)
        prof = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = (rr >= lo) & (rr < hi)
            if m.sum() >= 10:
                prof.append(dict(r=float((lo + hi) / 2), n=int(m.sum()), radial=float(np.median(rad[m])),
                                 tangential=float(np.median(tan[m]))))
        rep['absolute']['profile'] = prof
        print("    radial profile (image coordinates, all frames): r m -> median radial / tangential m")
        for p in prof:
            print(f"      {p['r']:6.0f}  n={p['n']:4d}   {p['radial']:+5.2f}  {p['tangential']:+5.2f}")
        json.dump(res, open(os.path.join(out, 'absolute_rows.json'), 'w'))
        picture(res, prof, rep, os.path.join(out, 'residuals.png'))
    json.dump(rep, open(os.path.join(out, 'report.json'), 'w'), indent=1)
    print(f"\nwrote {out}/report.json  [{time.time()-t0:.0f}s]")


def picture(res, prof, rep, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.4))
    x = np.array([r['x'] for r in res]); y = np.array([r['y'] for r in res])
    u = np.array([r['vx'] for r in res]); v = np.array([r['vy'] for r in res])
    # 4000 raw arrows are noise at any exaggeration; the systematic part is what
    # survives averaging over every frame, so draw the mean per 300 m bin
    B = 300.0; bx = np.floor(x / B); by = np.floor(y / B)
    X_, Y_, U_, V_ = [], [], [], []
    for kx, ky in set(zip(bx, by)):
        m = (bx == kx) & (by == ky)
        if m.sum() >= 8:
            X_.append((kx + 0.5) * B); Y_.append((ky + 0.5) * B)
            U_.append(np.median(u[m])); V_.append(np.median(v[m]))
    X_, Y_, U_, V_ = map(np.array, (X_, Y_, U_, V_))
    ax[0].quiver(X_, -Y_, U_ * 150, -V_ * 150, np.hypot(U_, V_), angles='xy', scale_units='xy',
                 scale=1, cmap='plasma', clim=(0, 2), width=0.006)
    ax[0].set_aspect('equal'); ax[0].set_title(
        f"{rep['tag']}: median residual per 300 m of image, all frames stacked\n"
        f"(per-frame similarity removed; arrows x150, colour 0-2 m)", fontsize=9)
    ax[0].set_xlabel('image x from centre, ground m'); ax[0].set_ylabel('image y from centre, ground m')
    rr = [p['r'] for p in prof]
    ax[1].plot(rr, [p['radial'] for p in prof], 'o-', label='radial (outward +)')
    ax[1].plot(rr, [p['tangential'] for p in prof], 's--', label='tangential')
    ax[1].axhline(0, color='k', lw=0.5)
    ax[1].set_xlabel('distance from frame centre, m'); ax[1].set_ylabel('median residual, m')
    ax[1].set_title('radial profile: a lens error rises with r^3; placement is flat', fontsize=9)
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=110)
    print(f"  wrote {path}")


if __name__ == '__main__':
    main()
