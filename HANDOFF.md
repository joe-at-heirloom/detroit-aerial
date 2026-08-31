# What we are trying to do

## The goal, in one sentence

Place scanned aerial photographs of Detroit from **1949, 1956, 1961 and 1967**
onto the map so precisely that when you wipe between any two years — or between
any year and modern satellite imagery — **the streets do not move.**

## The test of success

Pick any spot in the covered area, at any zoom, and wipe between any two layers.
A street corner must stay on the same screen pixel. Currently it doesn't
everywhere, and that is the whole remaining problem.

**Target: every cell under ~10 m. No cell over ~25 m.**

## How we measure it

Cut the mosaic into a **16 × 3 grid** of cells. For each cell, correlate a ridge
(road-like linear feature) response from the historical imagery against the same
response from modern imagery, and record how far it had to shift to match. That
shift is the error for that cell.

Measuring on coarse full-width chunks is what fooled us before — it averages the
left, middle and right of the block together and hides real lateral error. Always
report **median, p90, max, and the count of cells over 25 m**, never a single
average.

## Where it stands

| block | frames | mosaic | control stations | median | p90 | worst cell | cells >25 m |
|---|---|---|---|---|---|---|---|
| 1961 | 62 | 820 MP | ~130 | 15.3 m | 55.5 m | 231.9 m | 14 / 44 |
| 1967 | 51 | 825 MP | ~110 | 11.2 m | 43.0 m | 162.3 m | 8 / 31 |
| 1949 | 50 | 834 MP | 61 of 121 | not yet measured | | | |
| 1956 | 70 | 1070 MP | 75 of 159 | not yet measured | | | |

All four `mosaics/detroit_*_rbf.tif` exist (438-547 MB each). 1949 and 1956 have
been through pass one only and their grids have never been measured — do that
first, they may be worse than 1961/1967.

Note the low confident-station yield on the new blocks (61/121 and 75/159, versus
~85% on 1961). Coverage is thinner (62% and 67%), so there is less overlap with
modern imagery to correlate against. That is a likely weak point.

Downtown (a separate 3-frame scene, not part of the blocks) is good: 3–5 m.

## Four things already ruled out — do not repeat these

1. **Correlating against the residential street grid.** Detroit's side streets
   repeat every ~100 m, so a mosaic shifted by exactly one block scores as
   *perfectly aligned*. Every metric built on this lied, including one reporting
   "0.28 m bias" on 31,202 samples while the mosaic was 150 m out. Use the
   **mile-grid arterials** (1609 m spacing) for absolute orientation — they
   cannot alias inside a ±400 m search.

2. **Validating against road centreline vectors.** Modern satellite imagery
   scores the same 5–6 m against those centrelines as our historical imagery
   does, so that is the *metric's* noise floor, not our error. Validate against
   modern **imagery**, not vectors.

3. **A single global transform per block.** Right on average, wrong everywhere in
   particular. Per-chunk error stayed at 87–217 m.

4. **A translation-only bundle adjustment.** Aerial film has a per-frame crab
   angle (measured spread 3.0° in 1961, 6.4° in 1956). Ignoring it is what made
   the stitching look crooked. Rotation must be solved per frame.

## The pipeline as it stands

1. **Self-calibrate** each mission — scale and rotation derived from the imagery,
   never hardcoded (`pipeline/calib.py`).
2. **Relative orientation with rotation** — windowed phase correlation swept over
   ±3° per overlapping pair (`pipeline/tiesim.py`). Tie residuals 0.5–1.7 m.
3. **Bundle adjust** rotations, then positions.
4. **Absolute orientation** against mile-grid arterials
   (`pipeline/absorient.py`).
5. **Local warp** — control stations on a grid, correlating historical ridge
   response against modern ridge response, then a linear trend plus a
   Gaussian-kernel RBF fitted to them (`pipeline/stations.py`,
   `pipeline/rbfwarp.py`, `pipeline/rbfapply.py`). Length scale chosen by
   cross-validation on held-out stations.

## The current lead

A **second RBF pass**. Pass one solved stations against an *unwarped* mosaic, so
the search had to be wide (±330 m) and was ambiguous — bad stations produce a bad
warp. Re-solving on the already-warped result allows a tight ±90 m search, which
is unambiguous, so the stations are far cleaner and a shorter-length RBF (800 m)
can remove what pass one missed. Script: `/tmp/refine.py`, reports before/after
on the same 16 × 3 grid.

**Status: pass two is running on all four blocks as of this handoff.** If its
output is lost, just re-run `/tmp/refine.py 1961 1967 1949 1956` from the project
root (copy it into `scripts/` first — it is only in /tmp). It writes
`mosaics/detroit_*_p2.tif` and prints before/after for each block. Nothing is
wired into `data/manifest.json` yet, so the viewer still serves the pass-one
mosaics; swap the manifest `file` and `bbox` entries to the `_p2` versions only
if the numbers actually improve.

If that is not enough, the next steps in order:
- denser control stations, with per-station outlier rejection against neighbours
- iterate passes until the grid stops improving
- per-frame rather than per-region correction (the frames are the physical unit)
- orthorectification, which is the real ceiling — buildings lean differently in
  every negative, so rooftops can never align without a DEM or dense stereo

## Practical notes

- Run the viewer: `./scripts/run.sh` then <http://localhost:8770>
- Source frames cached at
  `/private/tmp/claude-501/.../scratchpad/detroit/fullres/` (~4 GB); if that is
  gone, re-download from Wayne State's ContentDM IIIF endpoint.
- Pipeline modules resolve data paths against `data/` — run scripts from the
  project root.
- `mosaics/` is ~1.2 GB and gitignored, inside iCloud-synced Documents.
