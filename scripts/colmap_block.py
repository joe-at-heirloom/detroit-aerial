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
    if reuse and not os.path.exists(f"{work}/images"):
        os.makedirs(work, exist_ok=True); os.symlink(reuse, f"{work}/images")
    os.makedirs(f"{work}/images", exist_ok=True)
    log = f"{work}/colmap.log"
    sol, ppm = placements(tag)
    recs = sorted(r for r in sol if sol[r].get('ok'))
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
    E0 = np.mean([sol[r]['e'] for r in recs]); N0 = np.mean([sol[r]['n'] for r in recs])
    with open(f"{work}/priors.txt", 'w') as pf:
        for r in recs:
            dst = f"{work}/images/{r}.jpg"
            if not os.path.exists(dst):
                im = Image.open(f"{SP}/fullres/{r}.jpg").convert('L')
                w, h = im.size; mx, my = int(w * margin), int(h * margin)
                im.crop((mx, my, w - mx, h - my)).save(dst, quality=95)
            pf.write(f"{r}.jpg {sol[r]['e']-E0:.2f} {sol[r]['n']-N0:.2f} 0.0\n")
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
    pairs = []
    for i, a in enumerate(recs):
        for b in recs[i + 1:]:
            d = math.hypot(sol[a]['e'] - sol[b]['e'], sol[a]['n'] - sol[b]['n'])
            if 1.0 < d < 0.75 * (rad[a] + rad[b]):
                pairs.append(f"{a}.jpg {b}.jpg")
    open(f"{work}/pairs.txt", 'w').write("\n".join(pairs) + "\n")
    print(f"  {len(pairs)} plausible pairs (of {len(recs)*(len(recs)-1)//2})", flush=True)
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
