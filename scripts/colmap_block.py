#!/usr/bin/env python
"""Solve a block's geometry with COLMAP instead of the hand-rolled bundle.

What this replaces: self-calibration, relative orientation, the per-frame
similarity bundle. What it keeps: the verified absolute measurement, the rigid
placement, and the seam metric as the acceptance test. COLMAP gives each negative
a full pose (position + 3-axis rotation, i.e. tilt) and a shared camera model
(focal, principal point, radial distortion) from SIFT ties -- the standard tool
doing the standard job -- and Detroit is flat enough that projecting each frame
through its pose onto one ground plane is orthorectification.

Steps (CLI only, no GUI):
  1. crop each scan's outer margin (film edge) into a work dir
  2. write catalogue positions as position priors, in a local metric frame
  3. feature_extractor  (one shared SIMPLE_RADIAL camera, focal prior)
  4. spatial_matcher    (match only neighbours, by the priors)
  5. mapper             (incremental SfM with bundle adjustment)
  6. model_aligner      (similarity onto the catalogue positions)
  7. model_converter    -> TXT  (cameras.txt, images.txt for rendering)

Usage:  ./.venv/bin/python scripts/colmap_block.py 1961 [--work /tmp/colmap_1961] [--margin 0.03]
        [--reuse-images DIR] [--focal 3602] [--k -0.00027]
"""
import sys, os, json, math, subprocess, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline')); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from PIL import Image
import dtmap
from validate import P
from rebuild import placements, SP
Image.MAX_IMAGE_PIXELS = None


def arg(n, d=None):
    return sys.argv[sys.argv.index(n) + 1] if n in sys.argv else d


def sh(cmd, log):
    print("$ " + " ".join(cmd), flush=True)
    with open(log, 'a') as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        raise SystemExit(f"failed: {cmd[0]} {cmd[1]} (see {log})")


def main():
    tag = sys.argv[1]
    work = arg('--work', f'/tmp/colmap_{tag}')
    margin = float(arg('--margin', 0.03))
    t0 = time.time()
    reuse = arg('--reuse-images')
    os.makedirs(f"{work}/images", exist_ok=True)
    log = f"{work}/colmap.log"
    sol, ppm = placements(tag)
    recs = sorted(r for r in sol if sol[r].get('ok'))
    # --plan FILE: solve an explicit set of negatives rather than whatever the old
    # catalogue happens to hold. {"recs": [...], "pos": {rec: [E, N]}, "pairs": [[a,b]]}.
    # A frame may appear with NO position -- the catalogue only ever covered the
    # frames someone had already located -- in which case it is still extracted and
    # matched, through `pairs`, and the reconstruction places it. E0/N0, the priors
    # file and the distance-plausible pair rule all use the positioned subset only.
    if arg('--plan'):
        pl = json.load(open(arg('--plan')))
        have = [r for r in sol if sol[r].get('ok')]
        gw = float(np.median([sol[r]['gw'] for r in have])) if have else 3200.0
        gh = float(np.median([sol[r]['gh'] for r in have])) if have else 3200.0
        pos = {str(k): v for k, v in (pl.get('pos') or {}).items()}
        newsol = {}
        for r in [str(x) for x in pl['recs']]:
            if r in sol and sol[r].get('ok'):
                newsol[r] = sol[r]
            else:
                p_ = pos.get(r)
                newsol[r] = dict(ok=True, e=(p_[0] if p_ else None), n=(p_[1] if p_ else None),
                                 rot=0.0, dE=0.0, dN=0.0, gw=gw, gh=gh, ratio=2.0)
        sol = newsol
        recs = [r for r in sorted(newsol, key=lambda x: int(x))
                if os.path.exists(f"{SP}/fullres/{r}.jpg")]
        miss = len(newsol) - len(recs)
        npos = sum(1 for r in recs if sol[r]['e'] is not None)
        print(f"  plan: {len(recs)} frames on disk ({miss} missing a scan), "
              f"{npos} with a position prior, {len(recs)-npos} without", flush=True)
    # --first N: a pilot on N consecutive frames of one flight line, to learn in
    # minutes whether COLMAP registers this film at all before spending an hour
    first = arg('--first')
    if first:
        import seamclass as SC
        lab = SC.lines_of(sol, recs)
        line0 = sorted([r for r in recs if lab[r] == 0], key=lambda r: sol[r]['n'])
        recs = line0[:int(first)]
    print(f"{tag}: {len(recs)} frames -> {work}", flush=True)

    # 1. crop the film edge; 2. position priors in a local metric frame (E, N, up)
    posrecs = [r for r in recs if sol[r].get('e') is not None]
    if not posrecs:
        raise SystemExit('no frame has a position prior; the block has nothing to anchor to')
    E0 = np.mean([sol[r]['e'] for r in posrecs]); N0 = np.mean([sol[r]['n'] for r in posrecs])
    # Crop the film edge, then PAD every crop to one common canvas, centred. The
    # scans come in five slightly different widths (5034-5088 px on 1961) at the
    # same resolution, and COLMAP with a single shared camera silently skips any
    # image whose size differs from the first -- 38 of 62 were dropped, and two
    # runs were misdiagnosed as a block that would not connect. Padding keeps the
    # scan's px/mm and the principal point exactly; resizing would not.
    crops = {}
    for r in recs:
        src = f"{reuse}/{r}.jpg" if reuse and os.path.exists(f"{reuse}/{r}.jpg") else None
        if src:
            crops[r] = Image.open(src).convert('L')
        else:
            im = Image.open(f"{SP}/fullres/{r}.jpg").convert('L')
            w, h = im.size; mx, my = int(w * margin), int(h * margin)
            crops[r] = im.crop((mx, my, w - mx, h - my))
    # Normalise scan RESOLUTION before padding. 1949 has a batch scanned square
    # (4998 x 4998) where the rest are 5076 x 4794: same 9 x 9 inch negative, ~4.5%
    # more pixels per millimetre. With one fixed focal for every negative those
    # frames are 4.5% wrong and the bundle bends the whole block to fit them --
    # measured as a -20% east-west affine against the old build. The film height
    # is the same physical dimension on every scan, so scale each crop to the
    # majority height; then pad widths, which are scan-window differences only.
    from collections import Counter
    hmaj = Counter(im.size[1] for im in crops.values()).most_common(1)[0][0]
    nres = 0
    for r, im in list(crops.items()):
        if abs(im.size[1] - hmaj) > 0.01 * hmaj:
            f_ = hmaj / im.size[1]
            crops[r] = im.resize((int(round(im.size[0] * f_)), hmaj), Image.LANCZOS); nres += 1
    if nres:
        print(f"  {nres} scans at a different resolution rescaled to height {hmaj}", flush=True)
    W_ = max(im.size[0] for im in crops.values()); H_ = max(im.size[1] for im in crops.values())
    sizes = sorted(set(im.size for im in crops.values()))
    print(f"  {len(sizes)} distinct crop sizes {sizes[0]}..{sizes[-1]}; padding all to {W_}x{H_}", flush=True)
    with open(f"{work}/priors.txt", 'w') as pf:
        for r in recs:
            dst = f"{work}/images/{r}.jpg"
            if not os.path.exists(dst):
                im = crops[r]
                if im.size != (W_, H_):
                    canvas = Image.new('L', (W_, H_), 0)
                    canvas.paste(im, ((W_ - im.size[0]) // 2, (H_ - im.size[1]) // 2)); im = canvas
                im.save(dst, quality=95)
            if sol[r].get('e') is not None:
                pf.write(f"{r}.jpg {sol[r]['e']-E0:.2f} {sol[r]['n']-N0:.2f} 0.0\n")
    del crops
    im = Image.open(f"{work}/images/{recs[0]}.jpg"); w, h = im.size
    # Camera: fixed. On flat ground focal and flying height trade off, so
    # self-calibrating focal is degenerate -- the full-block run produced 7198 px
    # for one component. The pilot on six frames solved 3602 px with the principal
    # point fixed, 0.6% from the 6-inch-lens prior; use that and refine nothing.
    focal_px = float(arg('--focal', 3602.0)); k = float(arg('--k', -0.00027))
    print(f"  images {w}x{h}, camera fixed: focal {focal_px:.0f} px, k {k:+.5f}", flush=True)
    # Pairs: only negatives that plausibly overlap by the catalogue. Exhaustive
    # matching on a street grid that repeats every 97.5 m matched non-overlapping
    # negatives confidently and split the block.
    rad = {r: math.hypot(sol[r]['gw'], sol[r]['gh']) / 2 for r in recs}
    seen, pairs = set(), []

    def add(a, b):
        k = (a, b) if a < b else (b, a)
        if a != b and k not in seen and a in rad and b in rad:
            seen.add(k); pairs.append(f"{k[0]}.jpg {k[1]}.jpg")

    for i, a in enumerate(posrecs):
        for b in posrecs[i + 1:]:
            d = math.hypot(sol[a]['e'] - sol[b]['e'], sol[a]['n'] - sol[b]['n'])
            if 1.0 < d < 0.75 * (rad[a] + rad[b]):
                add(a, b)
    n_dist = len(pairs)
    # Pairs the plan forces. For a roll of film, consecutive exposure numbers are
    # consecutive shots along a flight line and overlap about 60%, which is a
    # reliable pairing for frames no catalogue has ever placed.
    if arg('--plan'):
        for a, b in (json.load(open(arg('--plan'))).get('pairs') or []):
            add(str(a), str(b))
    open(f"{work}/pairs.txt", 'w').write("\n".join(pairs) + "\n")
    print(f"  {len(pairs)} pairs ({n_dist} by position, {len(pairs)-n_dist} forced) "
          f"of {len(recs)*(len(recs)-1)//2} possible", flush=True)
    json.dump(dict(E0=E0, N0=N0, focal_px=focal_px, w=w, h=h, margin=margin, recs=recs),
              open(f"{work}/meta.json", 'w'))

    db = f"{work}/db.db"
    # option names as of COLMAP 4.1 (FeatureExtraction.* / FeatureMatching.*).
    # 2800 px on a 3.2 km negative is 1.15 m/px for the tie features -- the pilot
    # closed to 2 m at 3600 -- and it is 1.65x cheaper per image on CPU, where
    # 62 images at 3600 took ~40 s each.
    if os.path.exists(db):
        os.remove(db)
    sh(['colmap', 'feature_extractor', '--database_path', db, '--image_path', f"{work}/images",
        '--ImageReader.single_camera', '1', '--ImageReader.camera_model', 'SIMPLE_RADIAL',
        '--ImageReader.camera_params', f"{focal_px:.1f},{w/2:.1f},{h/2:.1f},{k}",
        '--FeatureExtraction.max_image_size', '3600', '--SiftExtraction.max_num_features', '12000',
        '--FeatureExtraction.use_gpu', '0', '--FeatureExtraction.num_threads', '8'], log)
    import sqlite3
    n_db = sqlite3.connect(db).execute('select count(*) from images').fetchone()[0]
    print(f"  {n_db}/{len(recs)} images in the database after extraction", flush=True)
    if n_db < len(recs):
        raise SystemExit(f"failed: only {n_db} of {len(recs)} images extracted -- see {log}")
    sh(['colmap', 'matches_importer', '--database_path', db, '--match_list_path', f"{work}/pairs.txt",
        '--match_type', 'pairs', '--FeatureMatching.use_gpu', '0', '--FeatureMatching.num_threads', '8'], log)
    os.makedirs(f"{work}/sparse", exist_ok=True)
    sh(['colmap', 'mapper', '--database_path', db, '--image_path', f"{work}/images",
        '--output_path', f"{work}/sparse",
        '--Mapper.ba_refine_focal_length', '0', '--Mapper.ba_refine_principal_point', '0',
        '--Mapper.ba_refine_extra_params', '0'], log)
    # The mapper writes one model per connected component. Report them all, keep
    # the largest as sparse_txt. COLMAP's model_aligner is skipped: it needs a 3D
    # similarity from camera centres, which is degenerate for a single flight line
    # and fussy otherwise; colmap_render.py --self-align fits the ground plane and
    # a 2D similarity to the catalogue instead, and never fails on geometry.
    models = sorted(os.listdir(f"{work}/sparse"))
    counts = []
    for m in models:
        tmp = f"{work}/_m"; os.makedirs(tmp, exist_ok=True)
        sh(['colmap', 'model_converter', '--input_path', f"{work}/sparse/{m}",
            '--output_path', tmp, '--output_type', 'TXT'], log)
        n = sum(1 for l in open(f"{tmp}/images.txt") if l.strip() and not l.startswith('#')) // 2
        counts.append((n, m))
    counts.sort(reverse=True)
    print(f"  mapper produced {len(models)} model(s): " +
          ", ".join(f"model {m}: {n} frames" for n, m in counts), flush=True)
    best = counts[0][1]
    os.makedirs(f"{work}/sparse_txt", exist_ok=True)
    sh(['colmap', 'model_converter', '--input_path', f"{work}/sparse/{best}",
        '--output_path', f"{work}/sparse_txt", '--output_type', 'TXT'], log)
    n_img = counts[0][0]
    print(f"  registered {n_img}/{len(recs)} frames in the largest model; cameras.txt / "
          f"images.txt in {work}/sparse_txt  [{time.time()-t0:.0f}s]", flush=True)
    if n_img < len(recs) * 0.9:
        print(f"  WARNING: {len(recs)-n_img} frames not in the largest model -- the block "
              f"did not connect; see the other models", flush=True)


if __name__ == '__main__':
    main()
