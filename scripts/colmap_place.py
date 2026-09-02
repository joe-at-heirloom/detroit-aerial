#!/usr/bin/env python
"""COLMAP block -> placed mosaic, both numbers verified, in one command.

   render (2 m/px) + seam stats -> production mosaic -> absolute field -> fit ->
   continuous warp with seams verified on COLMAP-rendered frames

Usage:  ./.venv/bin/python scripts/colmap_place.py 1961 /tmp/colmap_1961 [--model auto]
"""
import sys, os, subprocess, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, '.venv', 'bin', 'python')


def run(args):
    print(f"\n$ {' '.join(args)}", flush=True)
    r = subprocess.run([PY, os.path.join(ROOT, 'scripts', args[0])] + args[1:], cwd=ROOT, text=True, capture_output=True)
    print(r.stdout[-6000:], flush=True)
    if r.returncode != 0:
        print(r.stderr[-3000:], flush=True); raise SystemExit(f"{args[0]} failed")


def main():
    tag, work = sys.argv[1], sys.argv[2]
    model = sys.argv[sys.argv.index('--model') + 1] if '--model' in sys.argv else None
    # --surface N : orthorectify onto the model's fitted ground surface (see
    # colmap_render.Surface); a domed block needs it. Off unless asked.
    surf = ['--surface', sys.argv[sys.argv.index('--surface') + 1]] if '--surface' in sys.argv else []
    cm = ['--model', sys.argv[sys.argv.index('--colmap-model') + 1]] if '--colmap-model' in sys.argv else []
    # --refine-ref REF : refine the block's heading/scale/shift against a reference
    # before rendering (see colmap_render.refine_against); the same reference is
    # then used for the absolute placement.
    rr = ['--refine-ref', sys.argv[sys.argv.index('--refine-ref') + 1]] if '--refine-ref' in sys.argv else []
    ref = ['--ref', rr[1]] if rr else []
    t0 = time.time()
    run(['colmap_render.py', work, '--self-align', '--tag', tag, '--out', '/tmp', '--mpp', '2.0'] + surf + cm + rr)
    run(['colmap_render.py', work, '--mosaic', '--self-align', '--tag', tag, '--mpp', '0.63'] + surf + cm + rr)
    run(['fieldmeasure.py', tag, 'colmap'] + ref)
    run(['fieldfit.py', tag, 'colmap', '--use', 'cells'])
    run(['fieldapply.py', tag, 'colmap', '--out', 'placedC'] + (['--model', model] if model else []) + ref)
    if ref:
        run(['fieldmeasure.py', tag, 'placedC'])      # and against modern, for the record
    print(f"\nDONE {tag}: mosaics/detroit_{tag}_placedC.tif   [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
