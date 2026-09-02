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
    t0 = time.time()
    run(['colmap_render.py', work, '--self-align', '--tag', tag, '--out', '/tmp', '--mpp', '2.0'])
    run(['colmap_render.py', work, '--mosaic', '--self-align', '--tag', tag, '--mpp', '0.63'])
    run(['fieldmeasure.py', tag, 'colmap'])
    run(['fieldfit.py', tag, 'colmap', '--use', 'cells'])
    run(['fieldapply.py', tag, 'colmap', '--out', 'placedC'] + (['--model', model] if model else []))
    print(f"\nDONE {tag}: mosaics/detroit_{tag}_placedC.tif   [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
