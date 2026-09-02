#!/usr/bin/env python
"""Render a COLMAP-solved block onto the ground plane, measure its seams, and
composite it -- the acceptance test for the simpler path.

Each negative gets projected through its solved camera (focal, principal point,
radial distortion, full 3-axis pose) onto one flat plane at the median height of
COLMAP's 3D points. Detroit is flat, so that is orthorectification. The same two
numbers as everywhere else in this project -- along-track and cross-line seam
disagreement, from dense tie windows -- are computed on the rendered frames, so
COLMAP's geometry is judged by the instrument that touches no reference.

Usage:  ./.venv/bin/python scripts/colmap_render.py /tmp/colmap_pilot [--mpp 2.0] [--tag 1961] [--out /tmp]
"""
import sys, os, json, math, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image, ImageDraw
import rasterio
from rasterio.transform import from_origin
import dtmap, close3, seamclass as SC
Image.MAX_IMAGE_PIXELS = None


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def quat_to_R(qw, qx, qy, qz):
    return np.array([[1 - 2*(qy*qy + qz*qz), 2*(qx*qy - qz*qw), 2*(qx*qz + qy*qw)],
                     [2*(qx*qy + qz*qw), 1 - 2*(qx*qx + qz*qz), 2*(qy*qz - qx*qw)],
                     [2*(qx*qz - qy*qw), 2*(qy*qz + qx*qw), 1 - 2*(qx*qx + qy*qy)]])


def self_align(cams, imgs, pts, priors):
    """Put an unaligned COLMAP model into the local metric frame ourselves.

    model_aligner needs a 3D similarity from camera centres, which is degenerate
    for cameras along one flight line (roll about the line is free). We do not
    need it: the ground is a plane, so fit that plane to the 3D points and take
    its normal as up; then a 2D similarity from the cameras' plane coordinates to
    the catalogue positions fixes scale, heading and translation, which six
    points along a line determine perfectly well. Ground goes to z = 0."""
    P = np.array(pts); c = P.mean(0); U, S, Vt = np.linalg.svd(P - c, full_matrices=False)
    nrm = Vt[2]
    C = np.array([imgs[n]['C'] for n in imgs])
    if np.dot(C.mean(0) - c, nrm) < 0: nrm = -nrm          # cameras above the ground
    # rotation taking nrm -> +z
    z = np.array([0, 0, 1.0]); v = np.cross(nrm, z); sn = np.linalg.norm(v); cs = float(np.dot(nrm, z))
    if sn < 1e-9: Rw = np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        Rw = np.eye(3) + vx + vx @ vx * ((1 - cs) / sn**2)
    zplane = float((Rw @ c)[2])
    names = [n for n in imgs if n in priors]
    A = np.array([(Rw @ imgs[n]['C'])[:2] for n in names]); B = np.array([priors[n][:2] for n in names])
    # Umeyama 2D similarity A -> B
    ma, mb = A.mean(0), B.mean(0); Ac, Bc = A - ma, B - mb
    H = Ac.T @ Bc / len(A); U2, S2, Vt2 = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt2.T @ U2.T)); D = np.diag([1, d])
    R2 = Vt2.T @ D @ U2.T; sc = float(np.trace(np.diag(S2) @ D) / (Ac**2).sum() * len(A))
    t2 = mb - sc * R2 @ ma
    Q = np.eye(3); Q[:2, :2] = R2; Q = Q @ Rw
    T = np.array([t2[0], t2[1], -sc * zplane])
    for n, im in imgs.items():
        Rn = im['R'] @ Q.T; tn = sc * im['t'] - Rn @ T
        im['R'], im['t'], im['C'] = Rn, tn, -Rn.T @ tn
    res = np.hypot(*(sc * (R2 @ A.T).T + t2 - B).T)
    print(f"  self-align: scale {sc:.4f}, heading {math.degrees(math.atan2(R2[1,0], R2[0,0])):+.2f} deg, "
          f"camera-vs-catalogue residual median {np.median(res):.0f} m (catalogue is approximate)", flush=True)
    # the 3D points in the aligned frame, for a ground SURFACE rather than a plane
    Pa = (sc * (Q @ P.T)).T + T
    self_align.points = Pa
    return 0.0


def refine_against(work, meta, cams, imgs, ref, mpp=5.0, search_m=600.0, log=print):
    """Refine the block's heading, scale and shift against a reference raster.

    The Procrustes to catalogue positions is only as good as the catalogue: 1949's
    scatter 361 m around the solved cameras, and the heading it fitted was a few
    degrees off, which is ~900 m at the block's ends -- far outside any fine
    search, and a mis-shape no later warp should be asked to absorb. So: render
    the block coarsely, measure its offset from the reference on the alias-proof
    arterial coarse field (6 km windows, wide-road ridge, +/-600 m search), fit a
    2D similarity to those offsets with robust rejection, and fold it into the
    camera alignment. Four parameters from ~30 windows: well determined, and it
    cannot bend anything."""
    import gridval
    from validate import reference_for
    E0, N0 = meta['E0'], meta['N0']
    C = np.array([imgs[n]['C'] for n in imgs]); half = 2600.0
    minE = C[:, 0].min() + E0 - half; maxE = C[:, 0].max() + E0 + half
    minN = C[:, 1].min() + N0 - half; maxN = C[:, 1].max() + N0 + half
    W = int((maxE - minE) / mpp); H = int((maxN - minN) / mpp)
    comp = np.zeros((H, W), np.uint8)
    for name, img in imgs.items():
        im = np.asarray(Image.open(f"{work}/images/{name}").convert('L'))
        r = render_frame(name, img, cams[img['cam']], im, 0.0, E0, N0, minE, maxN, W, H, mpp, crop=0.95)
        if r is None: continue
        a, x0, y0 = r; sub = comp[y0:y0 + a.shape[0], x0:x0 + a.shape[1]]
        np.copyto(sub, a, where=(a > 0) & (sub == 0))
    bbox = [dtmap.LAT0 + minN / dtmap.MLAT, dtmap.LON0 + minE / dtmap.MLON,
            dtmap.LAT0 + maxN / dtmap.MLAT, dtmap.LON0 + maxE / dtmap.MLON]
    refr = reference_for(bbox, W, H, ref)
    t = gridval.prepare(comp, refr, mpp, log=lambda *_: None)
    pts = gridval.coarse_field(t, search_m=search_m, log=log)
    if len(pts) < 8:
        log(f"  refine-against {ref}: only {len(pts)} windows locked; leaving the alignment as it is")
        return 0.0
    # window position (map m) -> where its content actually is: p + d
    Pw = np.array([[minE + p['x'] * mpp, maxN - p['y'] * mpp] for p in pts])
    D = np.array([[p['dE'], p['dN']] for p in pts])
    # content sits +d from the reference, so the block must move by -d
    keep = np.ones(len(D), bool)
    for _ in range(4):
        A_ = Pw[keep]; B_ = (Pw - D)[keep]
        ma, mb = A_.mean(0), B_.mean(0); Ac, Bc = A_ - ma, B_ - mb
        Hm = Ac.T @ Bc / len(A_); U2, S2, Vt2 = np.linalg.svd(Hm)
        d_ = np.sign(np.linalg.det(Vt2.T @ U2.T)); Dm = np.diag([1, d_])
        R2 = Vt2.T @ Dm @ U2.T; sc = float(np.trace(np.diag(S2) @ Dm) / (Ac**2).sum() * len(A_))
        t2 = mb - sc * R2 @ ma
        res = np.hypot(*((sc * (R2 @ Pw.T).T + t2) - (Pw - D)).T)
        s_ = max(np.median(res[keep]) * 1.4826, 5.0); keep = res < 3 * s_
    th = math.degrees(math.atan2(R2[1, 0], R2[0, 0]))
    log(f"  refine-against {ref}: {len(pts)} windows, {int(keep.sum())} kept; similarity shift "
        f"{t2[0]:+.0f},{t2[1]:+.0f} m (about origin), scale {(sc-1)*1e2:+.3f}%, rotation {th:+.3f} deg; "
        f"residual median {np.median(res[keep]):.1f} m", flush=True)
    # fold into the camera alignment: X' = sc*R2*X + t2 in the (E-E0, N-N0) frame
    Q = np.eye(3); Q[:2, :2] = R2
    T = np.array([t2[0] + sc * (R2 @ np.array([E0, N0]))[0] - E0, t2[1] + sc * (R2 @ np.array([E0, N0]))[1] - N0, 0.0])
    for n, im in imgs.items():
        Rn = im['R'] @ Q.T; tn = sc * im['t'] - Rn @ T
        im['R'], im['t'], im['C'] = Rn, tn, -Rn.T @ tn
    if hasattr(self_align, 'points') and self_align.points is not None:
        P_ = self_align.points; self_align.points = (sc * (Q @ P_.T)).T + T
    return 0.0


class Surface:
    """Ground height as a quadratic in (x, y), fitted to the aligned 3D points.

    A block bundled with a slightly wrong camera on flat terrain bows into a
    shallow dome -- measured on 1961: the point cloud's plane residual runs -51 to
    +36 m and camera heights follow a quadratic of -166 m over 10 km. The cameras
    are consistent with THAT surface, so projecting onto it closes the seams;
    projecting onto a plane through it left 12-16 m at every join. What the dome
    leaves is a smooth planimetric stretch, which the absolute warp removes."""

    def __init__(self, pts, order=2):
        P = np.asarray(pts); self.cx, self.cy = P[:, 0].mean(), P[:, 1].mean(); self.sc = 10000.0
        x = (P[:, 0] - self.cx) / self.sc; y = (P[:, 1] - self.cy) / self.sc
        cols = [np.ones_like(x), x, y] + ([x*x, x*y, y*y] if order >= 2 else [])
        A = np.stack(cols, 1)
        # robust: two passes trimming outliers (trees, roofs, mismatches)
        m = np.ones(len(P), bool)
        for _ in range(3):
            self.c = np.linalg.lstsq(A[m], P[m, 2], rcond=None)[0]
            r = A @ self.c - P[:, 2]; s = max(np.median(np.abs(r[m])) * 1.4826, 0.5)
            m = np.abs(r) < 3 * s
        self.order = order; self.resid = float(np.median(np.abs((A @ self.c - P[:, 2])[m])))

    def z(self, E, N):
        x = (np.asarray(E) - self.cx) / self.sc; y = (np.asarray(N) - self.cy) / self.sc
        cols = [np.ones_like(x), x, y] + ([x*x, x*y, y*y] if self.order >= 2 else [])
        return sum(c * v for c, v in zip(self.c, cols))


def read_model(d):
    cams = {}
    for l in open(f"{d}/cameras.txt"):
        if l.startswith('#') or not l.strip(): continue
        p = l.split(); cams[int(p[0])] = dict(model=p[1], w=int(p[2]), h=int(p[3]), params=[float(x) for x in p[4:]])
    imgs = {}
    lines = [l for l in open(f"{d}/images.txt") if not l.startswith('#') and l.strip()]
    for l in lines[0::2]:
        p = l.split()
        q = [float(x) for x in p[1:5]]; t = np.array([float(x) for x in p[5:8]])
        R = quat_to_R(*q); C = -R.T @ t
        imgs[p[9]] = dict(R=R, t=t, C=C, cam=int(p[8]))
    pts = []
    if os.path.exists(f"{d}/points3D.txt"):
        for l in open(f"{d}/points3D.txt"):
            if l.startswith('#') or not l.strip(): continue
            q = l.split(); pts.append((float(q[1]), float(q[2]), float(q[3])))
    zs = [p[2] for p in pts]
    return cams, imgs, (float(np.median(zs)) if zs else None), pts


def project(cam, R, t, X):
    """world points (N,3) -> pixel (u,v) for SIMPLE_RADIAL / PINHOLE / RADIAL."""
    x = (R @ X.T).T + t
    z = x[:, 2]; ok = z > 1e-6
    xn = x[:, 0] / np.where(ok, z, 1); yn = x[:, 1] / np.where(ok, z, 1)
    p = cam['params']
    if cam['model'] == 'SIMPLE_RADIAL':
        f, cx, cy, k = p; fx = fy = f
        r2 = xn*xn + yn*yn; d = 1 + k*r2; xn, yn = xn*d, yn*d
    elif cam['model'] == 'RADIAL':
        f, cx, cy, k1, k2 = p; fx = fy = f
        r2 = xn*xn + yn*yn; d = 1 + k1*r2 + k2*r2*r2; xn, yn = xn*d, yn*d
    elif cam['model'] in ('PINHOLE', 'SIMPLE_PINHOLE'):
        if cam['model'] == 'PINHOLE': fx, fy, cx, cy = p
        else: fx, cx, cy = p; fy = fx
    else:
        raise ValueError(cam['model'])
    return fx*xn + cx, fy*yn + cy, ok


def render_frame(name, img, cam, im, zg, E0, N0, minE, maxN, W, H, mpp, crop=0.95):
    """One negative on the common grid. `zg` is a constant ground height or a
    Surface; footprint estimated from the corners."""
    R, t = img['R'], img['t']
    surf = zg if isinstance(zg, Surface) else None
    z0 = float(surf.z(img['C'][0], img['C'][1])) if surf is not None else float(zg)
    w, h = cam['w'], cam['h']; f = cam['params'][0]; cx, cy = cam['params'][1], cam['params'][2]
    pts = []
    for (u, v) in ((w*(1-crop)/2, h*(1-crop)/2), (w*(1+crop)/2, h*(1-crop)/2), (w*(1+crop)/2, h*(1+crop)/2), (w*(1-crop)/2, h*(1+crop)/2)):
        d = R.T @ np.array([(u - cx)/f, (v - cy)/f, 1.0]); C = img['C']
        s = (z0 - C[2]) / d[2]; X = C + s*d; pts.append(X[:2])
    pts = np.array(pts)
    e0, e1 = pts[:, 0].min() + E0, pts[:, 0].max() + E0; n0, n1 = pts[:, 1].min() + N0, pts[:, 1].max() + N0
    x0 = max(0, int((e0 - minE)/mpp)); x1 = min(W, int((e1 - minE)/mpp) + 1)
    y0 = max(0, int((maxN - n1)/mpp)); y1 = min(H, int((maxN - n0)/mpp) + 1)
    if x1 - x0 < 8 or y1 - y0 < 8: return None
    Es = (minE + (np.arange(x0, x1) + 0.5)*mpp - E0).astype(np.float64)
    Ns = (maxN - (np.arange(y0, y1) + 0.5)*mpp - N0).astype(np.float64)
    EE, NN = np.meshgrid(Es, Ns)
    ZZ = surf.z(EE, NN) if surf is not None else np.full(EE.shape, z0)
    X = np.stack([EE.ravel(), NN.ravel(), np.asarray(ZZ).ravel()], 1)
    u, v, ok = project(cam, R, t, X)
    u = u.reshape(EE.shape); v = v.reshape(EE.shape); ok = ok.reshape(EE.shape)
    m = ok & (u >= w*(1-crop)/2) & (u < w*(1+crop)/2 - 1) & (v >= h*(1-crop)/2) & (v < h*(1+crop)/2 - 1)
    out = np.zeros(EE.shape, np.uint8)
    if not m.any(): return None
    ui = np.clip(u, 0, w - 2); vi = np.clip(v, 0, h - 2)
    ua = ui.astype(np.int32); va = vi.astype(np.int32); fu = (ui - ua).astype(np.float32); fv = (vi - va).astype(np.float32)
    val = (im[va, ua]*(1-fu)*(1-fv) + im[va, ua+1]*fu*(1-fv) + im[va+1, ua]*(1-fu)*fv + im[va+1, ua+1]*fu*fv)
    np.copyto(out, val.astype(np.uint8), where=m)
    return out, x0, y0


def load_block(work, model_dir=None, self_align_=True, surface=None, refine_ref=None):
    """Cameras, images (aligned into the local metric frame), ground z, offsets.
    surface: None for a plane, or the polynomial order of a fitted ground surface.
    refine_ref: a reference ('1961:placedC') to refine heading/scale/shift against."""
    if surface is None and '--surface' in sys.argv:
        surface = int(arg('--surface', 2))
    if refine_ref is None and '--refine-ref' in sys.argv:
        refine_ref = arg('--refine-ref')
    meta = json.load(open(f"{work}/meta.json"))
    model_dir = model_dir or (f"{work}/aligned" if os.path.exists(f"{work}/aligned/images.txt") else f"{work}/sparse_txt")
    cams, imgs, zg, pts = read_model(model_dir)
    if self_align_:
        priors = {}
        for l in open(f"{work}/priors.txt"):
            q = l.split(); priors[q[0]] = (float(q[1]), float(q[2]), float(q[3]))
        zg = self_align(cams, imgs, pts, priors)
        if refine_ref:
            zg = refine_against(work, meta, cams, imgs, refine_ref, log=print)
        if surface:
            order = int(surface)
            S_ = Surface(self_align.points, order=order)
            print(f"  ground surface: quadratic fitted to {len(self_align.points)} points, "
                  f"residual {S_.resid:.1f} m, height range over the block "
                  f"{S_.z(self_align.points[:,0].min(), self_align.points[:,1].min()) - S_.z(self_align.points[:,0].mean(), self_align.points[:,1].mean()):+.0f} m at a corner", flush=True)
            zg = S_
    return meta, cams, imgs, zg


def render_all_colmap(work, mpp, model_dir=None, crop=0.95, log=print, surface=None, refine_ref=None):
    """Every registered negative on one common grid at `mpp`, plus a placements-
    like dict (camera centres) so seam tools can label flight lines."""
    meta, cams, imgs, zg = load_block(work, model_dir, surface=surface, refine_ref=refine_ref)
    E0, N0 = meta['E0'], meta['N0']
    C = np.array([imgs[n]['C'] for n in imgs]); half = 2600.0
    minE = C[:, 0].min() + E0 - half; maxE = C[:, 0].max() + E0 + half
    minN = C[:, 1].min() + N0 - half; maxN = C[:, 1].max() + N0 + half
    W = int((maxE - minE) / mpp); H = int((maxN - minN) / mpp)
    rend = {}; sol = {}
    for name, img in imgs.items():
        rec = name.replace('.jpg', '')
        im = np.asarray(Image.open(f"{work}/images/{name}").convert('L'))
        r = render_frame(name, img, cams[img['cam']], im, zg, E0, N0, minE, maxN, W, H, mpp, crop=crop)
        if r is None: continue
        rend[rec] = r
        sol[rec] = dict(ok=True, e=img['C'][0] + E0, n=img['C'][1] + N0, dE=0.0, dN=0.0)
    return rend, sol, minE, maxN, W, H


def mosaic(work, tag, mpp=0.63, model_dir=None, crop=0.95, log=print):
    """Production composite: every pixel from the negative whose centre is
    nearest among those that cover it, at full resolution, streamed into a
    memmap so an 800 MP block never needs a float buffer of its own size."""
    from PIL import ImageDraw as _D
    meta, cams, imgs, zg = load_block(work, model_dir)
    E0, N0 = meta['E0'], meta['N0']
    names = sorted(imgs); C = np.array([imgs[n]['C'] for n in names]); half = 2600.0
    minE = C[:, 0].min() + E0 - half; maxE = C[:, 0].max() + E0 + half
    minN = C[:, 1].min() + N0 - half; maxN = C[:, 1].max() + N0 + half
    W = int((maxE - minE) / mpp); H = int((maxN - minN) / mpp)
    # coarse owner map: nearest camera centre among the frames whose footprint covers the pixel
    cm = 8.0; Wc = int(W * mpp / cm) + 1; Hc = int(H * mpp / cm) + 1
    owner = np.full((Hc, Wc), -1, np.int16); bestd = np.full((Hc, Wc), np.inf, np.float32)
    yy = (maxN - (np.arange(Hc) + 0.5) * cm)[:, None]; xx = (minE + (np.arange(Wc) + 0.5) * cm)[None, :]
    foot = {}
    for i, n in enumerate(names):
        img = imgs[n]; cam = cams[img['cam']]; w, h = cam['w'], cam['h']; f = cam['params'][0]; cx, cy = cam['params'][1], cam['params'][2]
        pts = []
        for (u, v) in ((w*(1-crop)/2, h*(1-crop)/2), (w*(1+crop)/2, h*(1-crop)/2), (w*(1+crop)/2, h*(1+crop)/2), (w*(1-crop)/2, h*(1+crop)/2)):
            d = img['R'].T @ np.array([(u - cx)/f, (v - cy)/f, 1.0]); Cc = img['C']
            zc = float(zg.z(Cc[0], Cc[1])) if isinstance(zg, Surface) else float(zg)
            s_ = (zc - Cc[2]) / d[2]; X = Cc + s_*d
            pts.append(((X[0] + E0 - minE) / cm, (maxN - (X[1] + N0)) / cm))
        m = Image.new('L', (Wc, Hc), 0); _D.Draw(m).polygon(pts, fill=255); m = np.asarray(m) > 0
        d2 = (xx - (img['C'][0] + E0))**2 + (yy - (img['C'][1] + N0))**2
        upd = m & (d2 < bestd); owner[upd] = i; bestd[upd] = d2[upd]
        foot[n] = pts
    log(f"  owner map {Wc}x{Hc} @ {cm} m; {len(names)} frames; output {W}x{H} @ {mpp} m/px")
    raw = f"/tmp/colmap_{tag}.raw"
    comp = np.memmap(raw, dtype=np.uint8, mode='w+', shape=(H, W)); comp[:] = 0
    k = cm / mpp
    for i, n in enumerate(names):
        img = imgs[n]; cam = cams[img['cam']]
        im = np.asarray(Image.open(f"{work}/images/{n}").convert('L'))
        # frame bbox in output px from its footprint
        P = np.array(foot[n]) * k; x0 = max(0, int(P[:, 0].min())); x1 = min(W, int(P[:, 0].max()) + 1)
        y0 = max(0, int(P[:, 1].min())); y1 = min(H, int(P[:, 1].max()) + 1)
        if x1 - x0 < 8 or y1 - y0 < 8: continue
        band = 1024
        for by in range(y0, y1, band):
            be = min(y1, by + band)
            r = render_frame(n, img, cam, im, zg, E0, N0, minE + x0 * mpp, maxN - by * mpp, x1 - x0, be - by, mpp, crop=crop)
            if r is None: continue
            arr, ax, ay = r
            oy = ((np.arange(by + ay, by + ay + arr.shape[0]) / k).astype(np.int32)).clip(0, Hc - 1)
            ox = ((np.arange(x0 + ax, x0 + ax + arr.shape[1]) / k).astype(np.int32)).clip(0, Wc - 1)
            own = owner[oy][:, ox] == i
            m = own & (arr > 0)
            sub = comp[by + ay:by + ay + arr.shape[0], x0 + ax:x0 + ax + arr.shape[1]]
            sub[m] = arr[m]
        if (i + 1) % 10 == 0: log(f"    {i+1}/{len(names)} frames")
    comp.flush()
    la_top = dtmap.LAT0 + maxN / dtmap.MLAT; lo_left = dtmap.LON0 + minE / dtmap.MLON
    tr = from_origin(lo_left, la_top, mpp / dtmap.MLON, mpp / dtmap.MLAT)
    from validate import P as _P
    out = _P('mosaics', f'detroit_{tag}_colmap.tif')
    with rasterio.open(out, 'w', driver='GTiff', height=H, width=W, count=1, dtype='uint8', crs='EPSG:4326',
                       transform=tr, tiled=True, blockxsize=512, blockysize=512, compress='DEFLATE', predictor=2,
                       num_threads='ALL_CPUS', BIGTIFF='YES') as ds:
        for y in range(0, H, 2048):
            yb = min(H, y + 2048)
            ds.write(np.asarray(comp[y:yb]), 1, window=rasterio.windows.Window(0, y, W, yb - y))
        ds.build_overviews([2, 4, 8, 16, 32, 64], rasterio.enums.Resampling.average)
    del comp
    try: os.remove(raw)
    except OSError: pass
    geo = dict(minE=float(minE), maxN=float(maxN), W=W, H=H, mpp=mpp,
               bbox=[dtmap.LAT0 + (maxN - H * mpp) / dtmap.MLAT, lo_left, la_top, dtmap.LON0 + (minE + W * mpp) / dtmap.MLON])
    json.dump(geo, open(_P('data', f'{tag}_colmap_geo.json'), 'w'))
    json.dump(dict(work=work, model_dir=model_dir,
                   surface=(int(arg('--surface', 2)) if '--surface' in sys.argv else None),
                   refine_ref=arg('--refine-ref')),
              open(_P('data', f'colmap_{tag}.json'), 'w'))
    log(f"  wrote {out}")
    return geo


def main():
    work = sys.argv[1]
    if '--mosaic' in sys.argv:
        tag = arg('--tag', 'colmap'); mpp = float(arg('--mpp', 0.63))
        mosaic(work, tag, mpp=mpp, model_dir=arg('--model'))
        return
    mpp = float(arg('--mpp', 2.0)); out_dir = arg('--out', '/tmp'); tag = arg('--tag', 'colmap')
    # one loader for every path, so --surface applies here as it does in mosaic()
    meta, cams, imgs, zg = load_block(work, arg('--model'), self_align_='--self-align' in sys.argv)
    E0, N0 = meta['E0'], meta['N0']
    if zg is None:
        raise SystemExit("no points3D.txt -- run model_converter with points")
    zdesc = 'fitted surface' if isinstance(zg, Surface) else f"z {float(zg):.1f}"
    print(f"{len(imgs)} registered frames, camera {cams[list(cams)[0]]['model']} {cams[list(cams)[0]]['params']}, ground {zdesc}")
    C = np.array([imgs[n]['C'] for n in imgs])
    zc = np.array([float(zg.z(c[0], c[1])) for c in C]) if isinstance(zg, Surface) else np.full(len(C), float(zg))
    H_fly = float(np.median(C[:, 2] - zc))
    print(f"  flying height above ground (solved) {H_fly:.0f} m", flush=True)
    # common grid
    half = 2600.0
    minE = C[:, 0].min() + E0 - half; maxE = C[:, 0].max() + E0 + half
    minN = C[:, 1].min() + N0 - half; maxN = C[:, 1].max() + N0 + half
    W = int((maxE - minE)/mpp); H = int((maxN - minN)/mpp)
    rend = {}; sol = {}
    for name, img in imgs.items():
        rec = name.replace('.jpg', '')
        im = np.asarray(Image.open(f"{work}/images/{name}").convert('L'))
        r = render_frame(name, img, cams[img['cam']], im, zg, E0, N0, minE, maxN, W, H, mpp)
        if r is None: continue
        rend[rec] = r
        sol[rec] = dict(ok=True, e=img['C'][0] + E0, n=img['C'][1] + N0, dE=0.0, dN=0.0)
    recs = sorted(rend); lab = SC.lines_of(sol, recs)
    print(f"  rendered {len(rend)} frames on {W}x{H} @ {mpp} m/px, {max(lab.values())+1} line(s)", flush=True)
    ties = close3.tie_windows(rend, lab, mpp, log=print)
    close3.seam_stats(ties, "COLMAP")
    # composite: nearest camera centre wins
    comp = np.zeros((H, W), np.uint8); best = np.full((H, W), np.inf, np.float32)
    for rec, (img, x0, y0) in rend.items():
        h, w = img.shape
        ys = (maxN - (np.arange(y0, y0+h) + 0.5)*mpp)[:, None]; xs = (minE + (np.arange(x0, x0+w) + 0.5)*mpp)[None, :]
        d = (xs - sol[rec]['e'])**2 + (ys - sol[rec]['n'])**2
        sub = best[y0:y0+h, x0:x0+w]; m = (img > 0) & (d < sub)
        sub[m] = d[m]; comp[y0:y0+h, x0:x0+w][m] = img[m]
    la_top = dtmap.LAT0 + maxN/dtmap.MLAT; lo_left = dtmap.LON0 + minE/dtmap.MLON
    tr = from_origin(lo_left, la_top, mpp/dtmap.MLON, mpp/dtmap.MLAT)
    p = os.path.join(out_dir, f"colmap_{tag}.tif")
    with rasterio.open(p, 'w', driver='GTiff', height=H, width=W, count=1, dtype='uint8', crs='EPSG:4326',
                       transform=tr, tiled=True, compress='DEFLATE') as ds:
        ds.write(comp, 1)
    print(f"  wrote {p}", flush=True)
    # a red/green seam picture: one along-track pair and one cross-line pair if present
    pairs = [(t['a'], t['b'], t['cross']) for t in ties]
    seen = {}
    for a, b, c in pairs:
        seen.setdefault(c, (a, b))
    panels = []
    for c in (False, True):
        if c not in seen: continue
        a, b = seen[c]; (ia, xa, ya), (ib, xb, yb) = rend[a], rend[b]
        x0 = max(xa, xb); y0 = max(ya, yb); x1 = min(xa+ia.shape[1], xb+ib.shape[1]); y1 = min(ya+ia.shape[0], yb+ib.shape[0])
        cx = (x0+x1)//2; cy = (y0+y1)//2; hh = int(700/mpp)
        def cut(im, xo, yo):
            o = np.zeros((2*hh, 2*hh), np.uint8); a0 = max(cy-hh, yo); a1 = min(cy+hh, yo+im.shape[0]); b0 = max(cx-hh, xo); b1 = min(cx+hh, xo+im.shape[1])
            if a1 > a0 and b1 > b0: o[a0-(cy-hh):a1-(cy-hh), b0-(cx-hh):b1-(cx-hh)] = im[a0-yo:a1-yo, b0-xo:b1-xo]
            return o
        def st(x):
            x = x.astype(np.float32); v = x > 0
            if v.sum() < 10: return np.zeros_like(x, np.uint8)
            lo, hi = np.percentile(x[v], [2, 98]); o = np.clip((x-lo)/max(hi-lo, 1)*255, 0, 255); o[~v] = 0; return o.astype(np.uint8)
        rgb = np.zeros((2*hh, 2*hh, 3), np.uint8); rgb[..., 0] = st(cut(ia, xa, ya)); rgb[..., 1] = st(cut(ib, xb, yb))
        im = Image.fromarray(rgb); d = ImageDraw.Draw(im); d.rectangle((0, 0, im.width, 20), fill=(0, 0, 0))
        d.text((6, 4), f"COLMAP {tag} {'CROSS-LINE' if c else 'ALONG-TRACK'}: {a} red / {b} green", fill=(255, 255, 255))
        panels.append(im)
    if panels:
        sheet = Image.new('RGB', (sum(p_.width for p_ in panels) + 10*(len(panels)-1), panels[0].height), (0, 0, 0))
        x = 0
        for p_ in panels: sheet.paste(p_, (x, 0)); x += p_.width + 10
        pp = os.path.join(out_dir, f"colmap_seams_{tag}.png"); sheet.save(pp); print(f"  wrote {pp}")


if __name__ == '__main__':
    main()
