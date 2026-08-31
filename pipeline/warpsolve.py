"""Solve the local displacement field by iteration, never aliasing.

The constraint that shapes everything here: Detroit's residential streets repeat
every 97.5 m in this block, so any single correlation search wider than about half
that can lock onto the wrong block and return a confident, wrong answer. The
project's pass-one control stations used a +/-330 m search and that is where the
bad control came from.

You cannot simply narrow the search, because the mosaic starts further out than
one block. So: never search wide, search repeatedly.

  1. a regional coarse field, solved on 6 km windows where the mile grid is
     actually present, gives a starting displacement good to a few tens of metres.
  2. every station then refines by at most +/-40 m -- narrower than half a block,
     so it cannot alias -- around whatever the current field says.
  3. the refined stations are fitted into a new field, and step 2 repeats.

Each pass can move a station by up to 40 m, so three or four passes cover a couple
of hundred metres of local error while every individual search stays unambiguous.
The iteration is in station space only: the raster is resampled once, at the end.
"""
import math, json
import numpy as np
import gridval, stations, rbfwarp


class PixelField:
    """Adapts a world-coordinate warp function to the pixel coords the matcher
    addresses windows in."""

    def __init__(self, warp_at, minE, maxN, mpp):
        self.warp_at = warp_at; self.minE = minE; self.maxN = maxN; self.mpp = mpp

    def at(self, y, x):
        d = self.warp_at(self.minE + x * self.mpp, self.maxN - y * self.mpp)
        return (d[0], d[1])


def solve(t, bbox, MLAT, MLON, win_m=2000.0, overlap=0.4, iters=4,
          min_valid=0.50, min_ratio=1.20, tol=2.0, lengths=None, log=print):
    """Returns (stations, trend, RBFWarp, warp_at)."""
    mpp = t['mpp']
    minE = (bbox[1] - (-83.0450)) * MLON
    maxN = (bbox[2] - 42.3340) * MLAT
    field = gridval.Prior(gridval.coarse_field(t, log=log), mpp)
    kept = None; tr = R = warp_at = None
    best = None            # (held-out median, kept, trend, RBF, warp_at)
    for it in range(iters):
        st = stations.solve_grid_prep(t, bbox, MLAT, MLON, win_m=win_m,
                                      overlap=overlap, min_valid=min_valid,
                                      prior=field)
        good = [s for s in st if s['ratio'] >= min_ratio]
        if len(good) < 20:
            log(f"  iter {it}: only {len(good)} confident stations, stopping")
            break
        kept, dropped = stations.reject_outliers(good, length=3000.0, k=2.5)
        moved = np.array([math.hypot(s['fine_dE'], s['fine_dN']) for s in kept])
        pegged = int((moved >= gridval.FINE_SEARCH - 1e-6).sum())
        log(f"  iter {it}: {len(st)} solved -> {len(good)} confident -> {len(kept)} kept"
            f"   moved this pass: median {np.median(moved):5.1f}  p90 "
            f"{np.percentile(moved,90):5.1f}  at search limit {pegged}")
        L = None; ho = float('inf')
        if lengths:
            sel, table = rbfwarp.cv_length(kept, lengths=lengths)
            L = sel[0]; ho = sel[1]
            log(f"    RBF length by held-out station error: {L:.0f} m "
                f"(held-out median {ho:.1f} m)")
        tr, R, warp_at = rbfwarp.fit(kept, length=L or 1000.0, trim=True)
        # Keep the best pass, not the last one. On the harder blocks the control is
        # noisy enough that a later pass can fit it worse -- 1949 went 15.2 -> 19.7 m
        # held-out across five passes -- and returning whichever happened to be last
        # ships the worse field.
        if best is None or ho < best[0]:
            best = (ho, kept, tr, R, warp_at)
        field = PixelField(warp_at, minE, maxN, mpp)
        if np.median(moved) < tol and pegged == 0:
            log(f"  converged after {it+1} passes")
            break
    if best is not None:
        if best[0] != float('inf'):
            log(f"  keeping the pass with the lowest held-out error: {best[0]:.1f} m")
        return best[1], best[2], best[3], best[4]
    return kept, tr, R, warp_at
