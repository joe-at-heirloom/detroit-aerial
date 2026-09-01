#!/usr/bin/env python
"""Picture of two seams: one along-track, one cross-line, red over green.

The single most informative view of this block: two negatives rendered on the
same map grid, one in red, one in green. Where they agree the picture is yellow;
where they disagree every street appears twice. It needs no reference imagery and
no metric, which is exactly why it is worth trusting.

Usage:  ./.venv/bin/python scripts/seampic.py 1961 stageA [--stage2 stageA2] [--out /tmp]
  With --stage2, the same two pairs are rendered from both placements side by side.
"""
import sys, os, json, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'pipeline'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import frameadjust
from validate import P
from rebuild import placements, SP
import close as C, seamclass as SC
from PIL import Image, ImageDraw

MPP = 2.0


def load(tag, stage):
    sol, ppm = placements(tag)
    saved = json.load(open(P('data', f'{stage}_{tag}.json')))
    for r in list(sol):
        v = saved.get(str(r)) or saved.get(r)
        if v:
            sol[r].update(dE=v['dE'], dN=v['dN'], rot=v.get('rot', sol[r]['rot']),
                          gw=v.get('gw', sol[r]['gw']), gh=v.get('gh', sol[r]['gh']))
    return sol


def stretch(x):
    x = x.astype(np.float32); v = x > 0
    if v.sum() < 10:
        return np.zeros_like(x, np.uint8)
    lo, hi = np.percentile(x[v], [2, 98])
    o = np.clip((x - lo) / max(hi - lo, 1) * 255, 0, 255); o[~v] = 0
    return o.astype(np.uint8)


def overlay(sol, r1, r2, span=1400):
    minE, maxN, W, H = C.canvas(sol, MPP)
    rr = [frameadjust.render_frame(r, sol[r], f"{SP}/fullres", minE, maxN, W, H, MPP, crop=0.95)
          for r in (r1, r2)]
    if any(v is None for v in rr):
        return None
    (i1, x1, y1), (i2, x2, y2) = rr
    ox0 = max(x1, x2); oy0 = max(y1, y2)
    ox1 = min(x1 + i1.shape[1], x2 + i2.shape[1]); oy1 = min(y1 + i1.shape[0], y2 + i2.shape[0])
    cx = (ox0 + ox1) // 2; cy = (oy0 + oy1) // 2; h = int(span / 2 / MPP)
    X0, X1, Y0, Y1 = cx - h, cx + h, cy - h, cy + h
    def cut(im, x0, y0):
        out = np.zeros((2 * h, 2 * h), np.uint8)
        a0 = max(Y0, y0); a1 = min(Y1, y0 + im.shape[0]); b0 = max(X0, x0); b1 = min(X1, x0 + im.shape[1])
        if a1 > a0 and b1 > b0:
            out[a0 - Y0:a1 - Y0, b0 - X0:b1 - X0] = im[a0 - y0:a1 - y0, b0 - x0:b1 - x0]
        return out
    A = stretch(cut(i1, x1, y1)); B = stretch(cut(i2, x2, y2))
    rgb = np.zeros(A.shape + (3,), np.uint8); rgb[..., 0] = A; rgb[..., 1] = B
    return Image.fromarray(rgb)


def pick_pairs(tag, sol):
    """The cross-line pair with the biggest sidelap, and an along-track neighbour."""
    recs = sorted(r for r in sol if sol[r].get('ok')); lab = SC.lines_of(sol, recs)
    dp = P('data', f'crossdiag_{tag}_stageA.json')
    cross = None
    if os.path.exists(dp):
        ok = [d for d in json.load(open(dp)) if d['why'] == 'ok']
        if ok:
            cross = max(ok, key=lambda d: min(d['ow'], d['oh']) * (1 if 60 < d['mag'] < 200 else 0.1))
            cross = (cross['a'], cross['b'])
    if cross is None:
        best = None
        for i, a in enumerate(recs):
            for b in recs[i + 1:]:
                if lab[a] == lab[b]:
                    continue
                d = math.hypot(sol[a]['e'] - sol[b]['e'], sol[a]['n'] - sol[b]['n'])
                if best is None or d < best[0]:
                    best = (d, a, b)
        cross = (best[1], best[2])
    a = cross[0]; along = None
    for r in recs:
        if r != a and lab[r] == lab[a] and abs(sol[r]['n'] - sol[a]['n']) < 1600:
            along = r; break
    return (a, along), cross


def main():
    tag, stage = sys.argv[1], sys.argv[2]
    stage2 = sys.argv[sys.argv.index('--stage2') + 1] if '--stage2' in sys.argv else None
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else '/tmp'
    sol = load(tag, stage)
    along, cross = pick_pairs(tag, sol)
    stages = [(stage, sol)] + ([(stage2, load(tag, stage2))] if stage2 else [])
    panels = []
    for st, so in stages:
        for nm, pr in (('ALONG-TRACK', along), ('CROSS-LINE', cross)):
            im = overlay(so, *pr)
            if im is None:
                continue
            d = ImageDraw.Draw(im)
            d.rectangle((0, 0, im.width, 22), fill=(0, 0, 0))
            d.text((6, 5), f"{tag} {st}  {nm}: {pr[0]} red / {pr[1]} green", fill=(255, 255, 255))
            panels.append(im)
    n = len(panels); cols = 2; rows = (n + 1) // 2
    w, h = panels[0].size
    sheet = Image.new('RGB', (cols * w + 10 * (cols - 1), rows * h + 10 * (rows - 1)), (0, 0, 0))
    for i, im in enumerate(panels):
        sheet.paste(im, ((i % cols) * (w + 10), (i // cols) * (h + 10)))
    p = os.path.join(out, f"seams_{tag}_{stage}{'_vs_' + stage2 if stage2 else ''}.png")
    sheet.save(p); print(p)


if __name__ == '__main__':
    main()
