#!/usr/bin/env python
"""Turn a crossmatrix run into the two things you read: a table and a map.

  table  median / p90 / count for every pair, the pair matrix of medians, the
         instrument's own error on each pair (planted-field tracking), and
         every reliable cell over 25 m with its location.
  map    one panel per layer against a reference: the layer's own imagery with
         an arrow per locked cell, exaggerated, coloured by size. Constant
         arrows are a placement; arrows that swirl or grow toward an edge are
         geometry no single shift can remove.

Usage:  ./.venv/bin/python scripts/xmreport.py runs/crossmatrix_west_adj.json [--ref esrihi]
            [--exag 40] [--out runs/xm_west]
"""
import sys, os, json, math, argparse
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = lambda *a: os.path.join(ROOT, *a)


def locked(p, min_ratio):
    return [c for c in p.get('cells', []) if 'skip' not in c and not c['pegged']
            and c['ratio'] >= min_ratio and c.get('track', 0.0) <= 5.0]


def table(d):
    ids = [l['id'] for l in d['layers']]
    lab = {l['id']: f"{l['label']}" + (' (hand)' if l.get('adj') else '') for l in d['layers']}
    M = {}
    for p in d['pairs']:
        s = p.get('summary')
        M[(p['a'], p['b'])] = M[(p['b'], p['a'])] = s
    out = [f"# Pair matrix: {d['group']} ({d['when']}, {d['grid']['mpp']} m/px, "
           f"{d['cell_m']:.0f} m cells, ratio >= {d['min_ratio']}, hand alignment {d['adjust']})", '',
           'Median displacement between the two layers over cells both cover and the '
           'instrument locked, in metres. A blank is no overlap or nothing locked.', '']
    out.append('| | ' + ' | '.join(lab[i] for i in ids) + ' |')
    out.append('|---' * (len(ids) + 1) + '|')
    for a in ids:
        row = []
        for b in ids:
            s = M.get((a, b))
            row.append('·' if a == b else (f"{s['median']:.1f}" if s and s.get('n') else ''))
        out.append(f"| **{lab[a]}** | " + ' | '.join(row) + ' |')
    out += ['', '| pair | cells | median | p90 | max | >25 m | bias E,N | instrument (median/p90) |',
            '|---|---|---|---|---|---|---|---|']
    rows = []
    for p in d['pairs']:
        s = p.get('summary')
        if not s or not s.get('n'):
            why = p.get('skip') or (f"0 of {s['tried']} locked" if s else '')
            rows.append((9e9, f"| {lab[p['a']]} / {lab[p['b']]} | {why} | | | | | | |"))
            continue
        tk = f"{s['track_median']:.1f} / {s['track_p90']:.1f}" if s.get('track_median') is not None else ''
        rows.append((s['median'], f"| {lab[p['a']]} / {lab[p['b']]} | {s['n']}/{s['tried']} | "
                     f"{s['median']:.1f} | {s['p90']:.1f} | {s['max']:.1f} | {s['over25']} | "
                     f"{s['bias_dE']:+.1f}, {s['bias_dN']:+.1f} | {tk} |"))
    out += [r for _, r in sorted(rows, key=lambda x: x[0])]
    big = []
    for p in d['pairs']:
        for c in locked(p, d['min_ratio']):
            if c['mag'] > 25:
                big.append((c['mag'], lab[p['a']], lab[p['b']], c))
    if big:
        out += ['', f'## Reliable cells over 25 m ({len(big)})', '',
                '| pair | lat, lon | dE, dN | ratio |', '|---|---|---|---|']
        for m, a, b, c in sorted(big, key=lambda x: -x[0])[:60]:
            out.append(f"| {a} / {b} | {c['lat']:.4f}, {c['lon']:.4f} | "
                       f"{c['dE']:+.1f}, {c['dN']:+.1f} | {c['ratio']:.2f} |")
    return '\n'.join(out) + '\n'


def maps(d, ref, exag, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    g = d['grid']; s, w, n, e = g['bbox']
    cdir = P('runs', 'cache', f"xm_{d['group']}_{g['mpp']}")
    lab = {l['id']: l['label'] + (' (hand)' if l.get('adj') else '') for l in d['layers']}
    panels = []
    for p in d['pairs']:
        if ref in (p['a'], p['b']) and p.get('summary') and p['summary'].get('n'):
            other = p['b'] if p['a'] == ref else p['a']
            sign = 1 if p['a'] == other else -1       # displacement of `other` relative to ref
            panels.append((other, sign, p))
    if not panels:
        return None
    k = max(1, int(round(g['W'] / 700)))
    ncol = min(6, len(panels)); nrow = -(-len(panels) // ncol)
    asp = (n - s) * 111132.0 / ((e - w) * 82300.0)
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.2 * asp * nrow + 0.6), squeeze=False)
    for ax in axs.flat:
        ax.axis('off')
    for ax, (lid, sign, p) in zip(axs.flat, panels):
        img = os.path.join(cdir, f'{lid}.img.npy')
        if os.path.exists(img):
            I = np.load(img, mmap_mode='r')[::k, ::k]
            ax.imshow(I, cmap='gray', extent=(w, e, s, n), aspect='auto',
                      interpolation='bilinear', vmin=0, vmax=255)
        C = locked(p, d['min_ratio'])
        x = np.array([c['lon'] for c in C]); y = np.array([c['lat'] for c in C])
        u = np.array([sign * c['dE'] for c in C]); v = np.array([sign * c['dN'] for c in C])
        m = np.hypot(u, v)
        ax.quiver(x, y, u * exag / 82300.0, v * exag / 111132.0, m, angles='xy', scale_units='xy',
                  scale=1, cmap='plasma', clim=(0, 20), width=0.006, headwidth=3)
        sm = p['summary']
        ax.set_title(f"{lab[lid]}  {sm['median']:.1f} m (p90 {sm['p90']:.0f})", fontsize=9)
        ax.set_xlim(w, e); ax.set_ylim(s, n)
    fig.suptitle(f"{d['group']}: each layer against {lab[ref]}; arrows x{exag}, colour 0-20 m",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('--ref', default=None)
    ap.add_argument('--exag', type=float, default=40.0)
    ap.add_argument('--out')
    a = ap.parse_args()
    d = json.load(open(a.run))
    base = a.out or os.path.splitext(a.run)[0]
    open(base + '.md', 'w').write(table(d))
    print(f"wrote {base}.md")
    ref = a.ref or ('esrihi' if any(l['id'] == 'esrihi' for l in d['layers']) else d['layers'][-1]['id'])
    pth = maps(d, ref, a.exag, f"{base}_vs_{ref}.png")
    if pth:
        print(f"wrote {pth}")


if __name__ == '__main__':
    main()
