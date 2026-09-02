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
    # focal prior: a 9x9 inch negative scanned at this size, 6 inch lens (the usual)
    px_per_mm = w / (228.6 * (1 - 2 * margin)); focal_px = 152.4 * px_per_mm
    print(f"  images {w}x{h}, focal prior {focal_px:.0f} px (6-inch lens)", flush=True)
    json.dump(dict(E0=E0, N0=N0, focal_px=focal_px, w=w, h=h, margin=margin, recs=recs),
              open(f"{work}/meta.json", 'w'))

    db = f"{work}/db.db"
    # option names as of COLMAP 4.1 (FeatureExtraction.* / FeatureMatching.*)
    sh(['colmap', 'feature_extractor', '--database_path', db, '--image_path', f"{work}/images",
        '--ImageReader.single_camera', '1', '--ImageReader.camera_model', 'SIMPLE_RADIAL',
        '--ImageReader.camera_params', f"{focal_px:.1f},{w/2:.1f},{h/2:.1f},0.0",
        '--FeatureExtraction.max_image_size', '3600', '--SiftExtraction.max_num_features', '12000',
        '--FeatureExtraction.use_gpu', '0', '--FeatureExtraction.num_threads', '8'], log)
    # exhaustive matching over ~60 images is ~1800 pairs, fine on CPU
    sh(['colmap', 'exhaustive_matcher', '--database_path', db, '--FeatureMatching.use_gpu', '0',
        '--FeatureMatching.num_threads', '8'], log)
    os.makedirs(f"{work}/sparse", exist_ok=True)
    sh(['colmap', 'mapper', '--database_path', db, '--image_path', f"{work}/images",
        '--output_path', f"{work}/sparse",
        '--Mapper.ba_refine_focal_length', '1', '--Mapper.ba_refine_principal_point', '0',
        '--Mapper.ba_refine_extra_params', '1'], log)
    models = sorted(os.listdir(f"{work}/sparse"))
    print(f"  mapper produced {len(models)} model(s): {models}", flush=True)
    best = f"{work}/sparse/{models[0]}"
    os.makedirs(f"{work}/aligned", exist_ok=True)
    sh(['colmap', 'model_aligner', '--input_path', best, '--output_path', f"{work}/aligned",
        '--ref_images_path', f"{work}/priors.txt", '--ref_is_gps', '0',
        '--alignment_type', 'custom', '--alignment_max_error', '150'], log)
    sh(['colmap', 'model_converter', '--input_path', f"{work}/aligned",
        '--output_path', f"{work}/aligned", '--output_type', 'TXT'], log)
    n_img = sum(1 for l in open(f"{work}/aligned/images.txt") if l and not l.startswith('#')) // 2
    print(f"  registered {n_img}/{len(recs)} frames; cameras.txt / images.txt in {work}/aligned  "
          f"[{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
