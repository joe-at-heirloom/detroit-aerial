#!/usr/bin/env python
"""One command from a closed placement to a placed mosaic, both numbers reported.

   composite stageA2 -> measure absolute field -> fit simplest model by held-out
   error -> apply as a continuous warp with seams verified -> seam pictures

Usage:  ./.venv/bin/python scripts/place2.py 1961 [--stage stageA2] [--model auto]
"""
import sys, os, subprocess, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, '.venv', 'bin', 'python')


def run(args, log):
    print(f"\n$ {' '.join(a for a in args)}", flush=True)
    r = subprocess.run([PY] + [os.path.join(ROOT, 'scripts', args[0])] + args[1:],
                       cwd=ROOT, text=True, capture_output=True)
    print(r.stdout[-6000:], flush=True)
    if r.returncode != 0:
        print(r.stderr[-3000:], flush=True); raise SystemExit(f"{args[0]} failed")


def main():
    tag = sys.argv[1]
    stage = sys.argv[sys.argv.index('--stage') + 1] if '--stage' in sys.argv else 'stageA3'
    model = sys.argv[sys.argv.index('--model') + 1] if '--model' in sys.argv else None
    label = {'stageA3': 'closedA3', 'stageA2': 'closedA2', 'stageA': 'closedA'}.get(stage, stage)
    out = {'stageA3': 'placed3', 'stageA2': 'placed2'}.get(stage, f'placed_{stage}')
    t0 = time.time()
    run(['applyclosed.py', tag, '--stage', stage, '--label', label], None)
    run(['fieldmeasure.py', tag, label], None)
    run(['fieldfit.py', tag, label, '--use', 'cells'], None)
    args = ['fieldapply.py', tag, label, '--out', out] + (['--model', model] if model else [])
    run(args, None)
    run(['seampic.py', tag, 'stageA', '--stage2', stage, '--out', '/tmp'], None)
    print(f"\nDONE {tag}: mosaics/detroit_{tag}_{out}.tif   [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
