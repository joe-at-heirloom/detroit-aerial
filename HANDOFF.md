# What we are trying to do

## The goal, in one sentence

Place scanned aerial photographs of Detroit from **1949, 1956, 1961 and 1967**
onto the map so precisely that when you wipe between any two years — or between
any year and modern satellite imagery — **the streets do not move.**

## The test of success

Pick any spot in the covered area, at any zoom, and wipe between any two layers.
A street corner must stay on the same screen pixel.

**Target: every cell under ~10 m. No cell over ~25 m.**

---

# How we measure it — and how we know the ruler is straight

Cut the mosaic into a **16 × 3 grid**. For each cell, correlate a ridge (road-like
linear feature) response from the historical imagery against the same response
from modern imagery, and record how far it had to shift. Report **median, p90,
max, and the count of cells over 25 m** — never a single average.

That much was already right. What was missing is that **the ruler was bent, and
every number the project has ever reported was measured with it.** So
`scripts/validate.py` now proves the metric before it reports anything:

1. **Planted-shift calibration.** Take the modern imagery, mask it to the
   historical footprint, move it by a known amount, and measure. Whatever it
   reports is the metric's own error. It must recover 0 m, and it must recover a
   planted 120 m too.
2. **Self-consistency.** Plant a shift on the *real* historical mosaic. Every
   cell's measurement must move by exactly that much. A metric that locks onto
   the wrong street lattice passes step 1 and fails this one.
3. **Only then, the measurement.**

Run it: `./.venv/bin/python scripts/validate.py 1961 1967`

---

# What was actually wrong

## 1. The validator itself was aliasing — the project's own trap, re-entering by the back door

Detroit's residential streets in this block repeat every **97.5 m** (measured off
the modern imagery: ridge autocorrelation r = 0.46–0.67 on the E–W axis). A
single-stage ±200 m search over a ridge map that contains side streets therefore
holds six alias positions per axis, and it will confidently return one of the
wrong ones.

The evidence was sitting in the old numbers. Cell (5,1) reported **dE = +97.5 m** —
exactly one block, to the decimetre. The "231.9 m worst cell" and the "98 m" cell
were the metric wandering, not the mosaic moving.

The project already knew not to do this for *absolute orientation*. It came back
through the validator.

**The fix — never search wide, search repeatedly.** No single correlation is
allowed to be ambiguous:

- a **regional coarse field**, solved on 6 km windows where the mile grid is
  actually present (1609 m spacing, cannot alias inside ±250 m), gives every cell
  a starting displacement;
- each cell then refines by at most **±40 m** — narrower than half a block, so
  aliasing is impossible — and that refinement is *repeated*, re-centring each
  time, so the reach extends to ~120 m while every individual search stays
  unambiguous.

Self-consistency went from *"up to 405 m out, 3 cells inconsistent"* to **0.0 m,
none**.

## 2. The modern reference was Web Mercator with a lat/lon box stamped on it

Comparing the old reference against a correctly reprojected one, band by band up
the raster, the north–south offset went

    -7.4   -13.8   -17.3   -17.2   -13.8   -6.4  m

— a **symmetric parabola**: zero at both ends, ~17 m at mid-span. That shape, at
that amplitude, is produced by exactly one thing: Mercator pixels addressed as if
they were linear in latitude, matched at the endpoints. The closed form for that
error over this span predicts **15.5 m**.

So every control point in the project was being fitted to a map bowed by up to
17 m through the middle. The warp dutifully absorbed the bow.

`scripts/fetchmodern.py` refetches Esri World Imagery at z16 (**1.77 m/px**,
against 3.54 before) and resamples **each output row from its own true Mercator
latitude**, onto this project's linear-latitude grid.

## 3. The reference did not cover 1949 or 1956 — and the code hid it

1949 runs ~3 km further south than the old reference, 1956 ~3 km further south
and west. **15.5 km² of 1949's imagery and 65.4 km² of 1956's had no modern
reference at all.** `modern_for` cropped to the overlap and then *resized the crop
to the full block extent* — silently stretching the reference across the whole
block. Anything measured for those two epochs against it was meaningless.

It now samples properly and leaves the uncovered margin as nodata, which the
masked correlation already knows to ignore. The z16 refetch covers all four
blocks.

## 4. A smooth field cannot fix a per-frame step

The mosaic is a Voronoi composite of individual negatives, and each negative
carries its own residual position and crab angle. Where two frames meet the error
can step discontinuously, and no continuous RBF can represent a step — it splits
the difference and leaves half the disagreement on both sides. Spatially-blocked
cross-validation put the floor on a smooth field at **~9 m held-out**, and frame
centres sit ~1.4 km apart, which is exactly the scale at which the field stopped
being predictable.

So the frame is corrected as a frame: each negative is rendered alone into map
space (middle 78% only — the outer margin is where vignetting and relief
displacement are worst), matched against modern imagery with the same alias-proof
search, and fed back as its own (dE, dN) before re-compositing.

A frame that does not lock **must not be left at zero** while its neighbours move
by 70 m — that tears the mosaic exactly where it used to be continuous. Frames
that fail borrow their neighbourhood's consensus; frames that disagree with their
neighbourhood are replaced by it.

## 5. Smaller things that were quietly wrong

- The "second peak" exclusion radius was fixed at 140 m, wider than the fine
  search window, so it blanked the whole surface and **every confidence ratio came
  back as ~0**.
- Peak position was integer-pixel, quantising every measurement to ±1.25 m for no
  reason. Now parabolic sub-pixel.
- The displacement lattice was sampled every 200 m regardless of the RBF length
  scale, which cross-validation keeps choosing at 400–600 m. It now scales with it.

---

# The pipeline

1. **Self-calibrate** each mission — scale and rotation from the imagery, never
   hardcoded (`pipeline/calib.py`).
2. **Relative orientation with rotation** — windowed phase correlation swept over
   ±3° per overlapping pair (`pipeline/tiesim.py`). Tie residuals 0.5–1.7 m.
3. **Bundle adjust** rotations, then positions.
4. **Absolute orientation** against the mile-grid arterials (`pipeline/absorient.py`).
5. **Composite** the block (`pipeline/mosaic_sim.py`).
6. **Per-frame adjustment** against modern imagery (`pipeline/frameadjust.py`),
   then re-composite.
7. **Residual local field** — control stations solved coarse-to-fine and iterated
   to convergence (`pipeline/warpsolve.py`), length scale chosen by
   *spatially-blocked* cross-validation (`pipeline/rbfwarp.py`), applied once
   (`pipeline/rbfapply.py`).

Steps 5–7 are one command: `./.venv/bin/python scripts/rebuild.py 1961`

---

# Ruled out — do not repeat these

1. **Correlating against the residential street grid.** It repeats every 97.5 m
   here, so a mosaic shifted exactly one block scores as perfectly aligned. One
   metric reported "0.28 m bias" on 31,202 samples while the mosaic was 150 m out.
2. **Validating against road centreline vectors.** Modern imagery scores the same
   5–6 m against them as our imagery does — that is the metric's noise floor.
3. **Per-cell arterial validation.** Modern imagery, which is zero-error by
   definition, scores 90–212 m against the arterials in the eastern column with
   ratio ≈ 1.0. One 3.6 × 2.0 km cell contains about two arterials per axis, which
   is not enough to lock. Useful only where bearing diversity is high, and only to
   tell "a few metres" from "a whole block".
4. **A single global transform per block.** Right on average, wrong everywhere in
   particular.
5. **Translation-only bundle adjustment.** Per-frame crab is 3.0° in 1961, 6.4° in
   1956; 1° over a 3.4 km frame throws the corners 30 m.
6. **Any correlation search wider than half a block that is not on the arterials.**
   This includes the previous session's `/tmp/refine.py`, which tightened the
   search to ±90 m — still wider than the 97.5 m block, so it still aliased. Its
   `mosaics/detroit_*_p2.tif` outputs are not used.
7. **Asking an ArcGIS ImageServer for a pixel grid whose aspect does not match
   the bbox.** `exportImage` does not simply honour the bbox: if the requested
   pixel grid has a different aspect ratio, it silently *expands the bbox* to
   match and returns imagery covering more ground than you asked for. Our check
   windows are square in metres, which is 1.35:1 in degrees at this latitude, so
   requesting a square pixel grid stretched the latitude extent by **243 m** --
   and the independent verification then compared two different footprints and
   reported ~74 m of error that did not exist. Three separate registration
   libraries all confirmed that phantom error, because all three were being handed
   the same mis-framed image. Request a grid whose aspect already matches the bbox
   in degrees, and assert the returned extent (`naipcheck.verify_extent`).

8. **Denser control stations.** The obvious lever for the ~1 km-scale residual is
   more control, so 1961 was re-solved with 1200 m station windows at 50% overlap
   instead of 2000 m at 40% -- 480 stations against 260. Held-out error got
   **worse**: 12.3 m against 5.9 m. A smaller window is a worse measurement, and
   more bad measurements do not make a better field. The cross-validation caught
   it; without that it would have looked like progress, because the in-sample grid
   would have improved.

9. **Chaining a block through an earlier epoch when the epochs do not share an
   extent.** Matching 1956 against the already-solved 1961 instead of modern
   imagery is the right idea -- five years apart rather than sixty-eight -- and the
   per-frame lock rate did improve, 53 of 70 frames. But 1956 runs several km
   further west and south than 1961, so outside 1961's coverage the reference falls
   back to modern, and frames straddling that boundary were matched against a
   patchwork of 1961 film and modern satellite. Held-out error got worse, 36.5 m
   against 27.8 m. Chaining needs the reference epoch to actually cover the block,
   or per-frame choice of which reference to use -- not a blended raster.

10. **Random-fold cross-validation on overlapping control windows.** A held-out
   station almost always has a near-duplicate left in the training set, so it
   scores interpolation and rewards ever-shorter length scales. Random folds
   claimed 5.1 m held-out where spatial folds said 9.1 m.

---

# Independent verification

Everything above is measured against a reference this project built itself. That
is a single point of failure, and it has been wrong once already. So the whole
chain is checked against imagery and tooling that share nothing with it:

- **USGS NAIP** (`scripts/naipcheck.py`) -- public-domain orthoimagery from
  USDA/USGS, a different organisation and a different processing chain than Esri.
  Its ImageServer returns EPSG:4326 on request, so *the server does the
  reprojection* and this project's own Mercator handling is taken out of the loop.
- **Four registration implementations** (`scripts/crosscheck.py`) -- ours, plus
  scikit-image's masked phase correlation, OpenCV's ECC (gradient-based, not
  correlation-peak-based) and SimpleITK's Mattes mutual information (an
  information-theoretic criterion that assumes nothing about the two images having
  similar brightness). Each is pinned by a planted-shift self-test before it is
  trusted -- that test caught scikit-image and SimpleITK returning the opposite
  sign to what their documentation implied.

Measured, our modern reference agrees with NAIP to **1.4-6.2 m** at widely
separated locations, all four backends agreeing; over 24 windows, median 2.16 m,
max 7.56 m.

**1961 is verified end to end.** Against NAIP, on 1.5 km windows, on the road
response, by three implementations that agree:

| backend | n | median | p90 | max | bias |
|---|---|---|---|---|---|
| ours | 10 | 3.84 m | 9.24 | 14.04 | +1.5, +1.7 |
| scikit-image | 8 | 3.16 m | 4.92 | 8.00 | 0.0, +2.0 |
| OpenCV ECC | 12 | 3.79 m | 12.23 | 16.09 | +0.2, +2.3 |

which matches this project's own 32 x 6 grid figure of 3.5 m. The residual ~2 m
northward bias is inside the reference's own agreement with NAIP, so it is not
distinguishable from the reference.

## Verify the verifier, every time

Two independent checks were themselves broken, and both looked convincing first:

- **The NAIP request was mis-framed** (see trap 7 below). Three separate
  registration libraries all agreed on ~74 m of error that did not exist, because
  all three were being handed the same wrongly-framed image. Agreement between
  tools is not evidence when they share an input.
- **AROSICS cannot measure this data.** It is a published, peer-reviewed
  co-registration package and it reported a beautifully tidy "median 0.51 m, bias
  0.00" across 349 tie points on the 1967 block. Then a planted 30 m shift was put
  on the target and it reported the same thing: its global mode returns `None` for
  two of three planted shifts and reliability 0 for the third, its local mode has
  median reliability 0. It is built for multi-sensor satellite imagery of the same
  era, and 1967 panchromatic film against modern orthoimagery is outside what it
  can match. **Its number was meaningless and would have been reassuring.**

So no backend counts as verification until it has recovered a *planted shift on the
real data* -- not on a synthetic pair, which only pins its sign convention.
`scripts/crosscheck.py --realcontrol` is that test, and it is the gate that decides
which tools are allowed to have an opinion.

What passed it, on 1967 film against NAIP (median error recovering the plant):

| backend | median error | within 5 m | verdict |
|---|---|---|---|
| ours (masked NCC on a ridge response) | 0.43 m | 71% | can measure this |
| scikit-image masked phase correlation | 1.00 m | 79% | can measure this |
| OpenCV ECC | 0.27 m | 67% | can measure this |
| SimpleITK Mattes mutual information | 43.69 m | 0% | **cannot** |
| AROSICS | did not respond to a 30 m plant | -- | **cannot** |

Two independent, validated verifiers therefore remain: **scikit-image** and
**OpenCV**. Of the two, only OpenCV's ECC is a *local* method, so for absolute
offsets it is the one to read -- scikit-image's phase correlation searches the
whole image and will happily lock a block off on this street grid. It was fine in
the control only because that measures a *change*, which cancels a consistent
alias.

# Tools

| | |
|---|---|
| `scripts/validate.py` | calibrate the metric, then measure a block |
| `scripts/rebuild.py` | composite → per-frame → residual field → final mosaic |
| `scripts/rewarp.py` | field-only re-solve, on the raw composite or on an existing warp |
| `scripts/montage.py` | alignment sheet: historical in red, modern in green, spread over the block |
| `scripts/inspectcells.py` | the same, for the worst cells specifically, as-built vs corrected |
| `scripts/fetchmodern.py` | refetch the modern reference, correctly reprojected |
| `scripts/manifest.py` | point the viewer at whatever the best build now is |
| `scripts/naipcheck.py` | verify against USGS NAIP -- independent imagery |
| `scripts/crosscheck.py` | measure the same windows four ways, three of them third-party |
| `pipeline/glue.py` | LightGlue frame matching (measured worse than ridge here; see the module) |
| `scripts/arosicscheck.py` | AROSICS run (kept for the record; it cannot match this data — see above) |

Where both epochs have visible roads, aligned streets render **yellow** and a
misalignment splits into a red ghost beside a green one. Areas that are red-only
or green-only are land-use change, not error — a freeway that did not exist in
1961 has no 1961 counterpart to align to.

---

# Practical notes

- Run the viewer: `./scripts/run.sh` then <http://localhost:8770>
- Source frames cached at
  `/private/tmp/claude-501/-Users-joelint-Documents-Apps/b1877f49-.../scratchpad/detroit/fullres/`
  (~4 GB, 238 frames); if that is gone, re-download from Wayne State's ContentDM
  IIIF endpoint.
- Ridge maps are memoised under `/tmp/das_ridge`, keyed on the reference imagery
  as well as the raster, so changing the reference invalidates them.
- Pipeline modules resolve data paths against `data/` — run scripts from the
  project root.
- `mosaics/` is gitignored and lives inside iCloud-synced Documents. **Two
  mosaics disappeared from disk mid-session** (`detroit_1961_rbf.tif`,
  `detroit_1967_rbf.tif`) — not iCloud placeholders, simply gone. Everything is
  regenerable from `data/` plus the cached frames, but do not treat `mosaics/` as
  durable storage.
