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

# READ THIS FIRST — the project was measuring the wrong thing

Every accuracy figure this project ever produced answers one question: does this
AREA sit in the right place against modern imagery. None of them answered the
question a viewer notices first: **does a road stay straight across the join
between two negatives.**

Those are different properties, and the project optimised the first while
destroying the second. Measured pair by pair, straight out of the bundle
adjustment, before any correction:

| block | overlapping frames disagree with each other by |
|---|---|
| 1961 | median **16.5 m**, p90 37.1, max 65.7 — 31 of 44 pairs over 10 m |
| 1949 | median **22.6 m**, p90 60.5, max 73.5 — 36 of 46 pairs over 10 m |

1961 is the block that measures 1.3 m against modern imagery and is independently
confirmed at 3.2-3.8 m by three registration libraries. It is also torn. **The
bundle does not close, for any block.**

The per-frame correction stage then made it much worse — it moved frames
independently by up to 300 m chasing modern imagery, taking 1949 from 22.6 m to
**56.4 m** of frame-to-frame disagreement while every number being watched improved.

`scripts/seamcheck.py` is the metric that catches this. It uses no external
reference, no CRS and no third-party library: it renders each overlapping pair with
its final placement and correlates the two in the region they share. Perfectly
consistent frames give zero. **It is the acceptance test, and it must never feed the
pipeline.**

## What the disagreement is made of

`scripts/characterise.py` splits each overlap into sub-windows instead of measuring
it once, because the two candidate causes want opposite responses. A blunder gives a
roughly CONSTANT offset across the overlap and is fixable by re-solving. Geometry —
per-frame tilt, film deformation — VARIES across the overlap and no per-frame shift
can represent it.

Measured on 1961, across 14 pairs: constant per-pair offset median **15.8 m**,
variation within each overlap median **8.7 m** (p90 17.3). Within one overlap the
error ran 6.9 → 21.3 → 25.0 m across the band.

So it is roughly two-thirds blunder, one-third geometry. Re-solving the block can
remove the first. The second sets a floor near **8-9 m** until there is a camera
model. **There are no fiducials** — the scans are cropped to the image area, with
annotation burned into the frame and no film margin — so the classic archival
interior-orientation route is closed.

## Stage A works. Placing the closed block does not — and here is exactly why

`scripts/close.py` Stage A closes a block using only frame-to-frame observations,
and it works, verified end to end:

| | frame-to-frame disagreement |
|---|---|
| bundle output | median **20.6 m**, p90 71.7 |
| after closing | median **0.19 m**, p90 28.4 |

The solver was checked against planted offsets (exact recovery) and against planted
97.5 m blunders (12 of 12 caught) BEFORE it was run on real data. On the first pass
it rejected 10 of 64 observations as blunders. The closure survives compositing: the
composited closed block measures 0.19 m, and a pure translation leaves it at 0.18 m,
as a translation must.

**Placing that closed block absolutely is unsolved.** The closed block sits ~72 m
off, and the error is not uniform: it is a **bend of 94 m east and 180 m north end
to end over 26 km**. That is the textbook behaviour of relative-only adjustment on a
long strip — local pairwise constraints fix local shape exactly and then accumulate
drift down a 50-link chain, giving a chain of perfectly-fitted links curving away
from reality.

Three placement attempts failed, each producing a plausible small number rather than
a visible failure:

1. **Similarity fit, wrong coordinate frame.** Frames rendered into a canvas from the
   current solution, matched against a modern ridge map built on a different grid.
   Reported a confident 9.4 m shift that meant nothing.
2. **Similarity applied to frame centres only.** A similarity rotates and scales the
   frames too; moving only centres shears the block. Consistency 0.2 -> 4.4 m.
3. **Similarity with an alias-proof prior.** Fitted a 0.78% scale that was not real
   scale but the four-parameter fit absorbing residual shape error. Consistency
   0.2 -> 21.2 m, worse than the bundle it started from.
4. **Translation plus a measured degree-1 bend.** The bend is real and was measured
   before being fitted. But removing 180 m of drift over ~50 frames should cost about
   3.6 m of seam disagreement per adjacent pair, and it cost **12.9 m** while
   improving absolute only 72 -> 45 m. Over three times the cost the arithmetic
   allows means it is absorbing error, not removing it. Reverted.

### The circularity that makes this hard

Measuring the drift requires matching the block against modern imagery. That
matching aliases on Detroit's 97.5 m residential lattice unless it has a good prior.
A good prior requires a well-placed block. **A badly bent block cannot measure its
own bend.**

### The way out, not yet tried

Use absolute control that *cannot* alias. The mile-grid arterials repeat every
1609 m, so correlation against them is unambiguous inside a +/-400 m search — this is
the one absolute signal in this scene immune to the trap that has caused most of the
project's failures. `pipeline/absorient.py` already implements it and this session
demoted it in favour of ridge-vs-modern-imagery, which aliases.

The directive: constrain the closed strip's drift at several stations along its
length using arterial correlation, not imagery correlation. Few parameters, smooth
along the strip, and every candidate correction judged by BOTH numbers — a
correction that costs more seam disagreement than the drift arithmetic predicts is
absorbing error and must be rejected.

## Stage A closed each line, not the block (found 2026-09-01)

Split by pair class, the "0.19 m" closure was 47 along-track pairs hiding 16
cross-line pairs:

| pair class | measured / geometric | median | p90 | max |
|---|---|---|---|---|
| along-line | 47 / 147 | 0.1 m | 0.4 | 28 |
| cross-line | 16 / 241 | 16.5 m | 77 | 109 |

1961 is four parallel N-S flight lines, 15-16 frames each, 2.3 km apart with 1.1 km
sidelap. Each line is closed perfectly along its length and the lines disagree with
each other by 50-109 m. Along-track ties cannot observe a line's scale -- every
frame in the line shares it -- only sidelap ties can, and 225 of 241 sidelap pairs
failed the confidence filter and never entered the solve. Never again report one
seam median: `scripts/seamclass.py` reports the two classes.

Measured against modern, the lines' along-track scale errors alternate:
-1.01%, -0.59%, -1.04%, -0.57% -- the signature of alternating flight direction.
Measured properly (`scripts/crossdiag.py`, every sidelap pair regardless of
threshold, search +/-150 m): 151 of 241 pairs hit the search EDGE, 40 fail on
coverage, 9 are weak, and the 41 that lock disagree by **median 114 m, max 195 m**.
The lines sit 100-200 m apart from each other. The same disease on every block:
1949 cross-line 65 m median (9 of 172 pairs measurable), 1967 54-90 m. And the
iteration-1 affine warp of 1961 confirms it from the other side: it took absolute
error from 73 to 43 m and no further, because the 43 m that remains IS the lines
disagreeing, which no global model can express.

The fix is `scripts/close2.py`: Stage A with a scale and a rotation per flight
line as unknowns alongside the per-frame translations, solved from along-track and
cross-line observations together, cross-line searched at +/-250 m. Verified on a
synthetic two-line block with a planted 0.5% differential scale: recovered exactly,
residual 0.000 m. A per-line similarity (4 lines x 4 params) is the physical model.

**The sidelap needs a different matcher than the along-track overlap.** With a
+/-250 m normalised-correlation search, 1967's 61 cross-line "matches" came back at
a median of 239 m -- the search boundary -- because the normalised correlation of
the few pixels still overlapping at a large shift is high by chance. The solver
then fitted those blunders (per-line scales of 5%, frames moved 700 m); the runs
were killed before anything was applied, and close2.py now refuses any solution
asking for more than 1.5% of line scale or moving frames further than twice the
disagreement it is removing. Cross-line offsets are now measured the way the
original tie-point stage measured along-track ones: **multi-window phase
correlation** on the raw film (same epoch, so a delta peak), accepting a pair only
when several 800 m windows independently agree. Verified on a planted 190 m shift:
36 of 36 windows agree.

On real sidelap the limiting factor is coverage, not peak strength (pcdiag.py, 30
sidelaps): every window that reached the peak test passed it at median sharpness
15-20 against a threshold of 7, but 78% never got there because one frame had
under 60% content in them. Each frame's bounding box is sized by its diagonal, so
two sidelapping boxes intersect ~2.2 km wide where the film overlaps ~0.9 km, and
half of every window sits in empty box. 800 m windows accepted 3 of 30 sidelaps,
512 m accepted 10; 384 m with a 75% coverage floor is the setting -- tiling the
film's own extent rather than the bounding box, with a boundary rule at 60% of the
window so a real 185 m offset is not rejected.

**close2's per-line model is the wrong shape, measured.** On 1961 it took the
cross-line disagreement from 63 to 21 m by pushing along-track from 0.2 to 9 m:
least squares sacrificed 50 good observations to serve 121 the model could not fit.
The lines are not rigidly scaled and rotated relative to each other.

The structural mistake was one offset per pair. A single offset determines a
translation and nothing else -- scale and rotation are visible only as the offset
VARYING across the overlap, and collapsing the phase-correlation windows to a median
throws exactly that away. `scripts/close3.py` keeps every window as its own tie
point, along-track and cross-line alike (thousands per block), and solves a
translation, scale and rotation per FRAME: along-track windows pin neighbouring
frames' scales together (0.2 m over a 2 km overlap allows 0.03% of difference),
cross-line windows pin the lines to each other. This is the textbook tie-point
bundle adjustment the panel recommended. `scripts/test_close3.py` plants per-frame
scale and rotation on a synthetic block and recovers them to 0.02% with zero
residual; run it before trusting any change.

**The bundle's per-frame rotations may be ~1 deg wrong, and nothing before close3
could see it.** `tiesim.match` estimated each frame's crab from one 384 px window
at the overlap centre, swept in 0.5 deg steps, keeping the sharpest peak -- and
the overlap centre is exactly where rotation has no effect. The "0.5 m tie
residual" was measured there. Dense tie windows across the whole overlap show
along-track disagreement of 8 m median, 27 m p90, inside overlaps whose central
offset is 0.2 m: that is the signature of ~1 deg of rotation between adjacent
frames, and close3's solve asks for a median of 1.0 deg per frame. Whether those
rotations are real is decided by re-measuring after applying them (round 1), not
by a threshold, which is why close3's guard bounds scale but not rotation.

**close3 works, measured on 1961.** After one round of the per-frame tie-point
bundle, along-track windows went 8.2 -> 4.0 m median (p90 27 -> 10) and cross-line
windows 67.2 -> 5.7 m (p90 159 -> 14), both improving together. Converged after
three rounds: along-track 4.0 m (p90 8.2, max 45), cross-line 4.5 m (p90 11.3,
max 44), over ~9000 tie windows. That 4 m floor is per-window variation inside an
overlap -- tilt, relief, film -- which a similarity per frame cannot remove; it is
the published floor, not hidden.

Same solve on the other blocks (Stage A -> close3, no absolute placement yet):

| block | along-track | cross-line | note |
|---|---|---|---|
| 1961 | 4.0 m (p90 8) | 4.5 m (p90 11) | converged |
| 1967 | 2.0 m (p90 6) | 4.0 m (p90 38, max 132) | a tail of bad sidelaps; round 2 refused by the guard |
| 1949 | 2.0 m (p90 6) | 20.0 m (p90 36) | cross-line STUCK: only ~46 of 172 sidelap pairs yield ties on this 3-line block |
| 1956 | 4.0 m (p90 9) | 5.7 m (p90 26, max 184) | round 2 refused (3.8% scale asked); closed on round 1's state |

1949's remaining 20 m, diagnosed (`scripts/pairdiag.py`): the worst cross-line
pairs yielded 2-3 tie windows each -- no constraint -- and the pairs that yielded
many had windows scattering by +/-45 m, which is noise. Grainy 1949 film over open
land gives raw-texture phase correlation nothing to lock to. `close3.py --feature
ridge` phase-correlates the road-ridge response instead; a window with no road is
flat and skipped rather than wrong. A few along-track pairs (1398/481, 631/708)
also sit at 20-38 m with high spread: likely tilt across the overlap, which a
similarity per frame cannot express. The ridge-feature pass took 1949's cross-line
from 20 to 14 m (95 windows) -- modest; the block's sidelap simply carries little
that phase-correlates. This is where COLMAP's real camera model should do better
than any per-frame similarity, and 1949 is queued behind 1961 on that path.

Two things to look at next, both about the guard-refused rounds: on a block that is
already closed the solver still proposes 20-36 m moves with 2-4% scales, which
means some frames' normal equations are near-singular -- block-edge frames with
ties on one side only, where translation trades off against scale. A stronger
prior for frames with few ties, or dropping their scale/rotation to the
neighbours' consensus, would let round 2 refine instead of being refused. And
1949's sidelap yields few ties; a per-pair breakdown of tie counts will say
whether the 20 m is a handful of bad pairs or the block's real floor. The rotations were
real: the block's per-frame crab was ~1 deg off out of the bundle. Nothing was
composited from any earlier Stage A variant.

**Run `scripts/test_close2.py` before trusting any change to close2.py.** It plants
four sidelap shifts from 45 to 225 m and a 0.5% differential line scale, and
derives the expected sign from how it moved the content. The matcher was right and
the hand-written expectation wrong about north three times in one session.

Known inefficiency, deliberately not fixed mid-run: `frameadjust.render_frame`
sizes every frame's canvas box by its half-diagonal so any rotation fits, which at
the actual crab angles (under 3 deg) renders about twice the pixels the footprint
needs and inflates every overlap box. Tightening it to the rotated extent would
roughly halve Stage A's per-round cost and raise sidelap coverage for free. It
touches every seam measurement, so do it between runs, and re-verify seams after. Its held-out
residual is still ~40 m, but that is measurement noise, not model error: the
per-frame absolute observations jump 26-38 m between ADJACENT frames that Stage A
closed to 0.1 m, which no geometry can produce. The noise has a cause -- a 0.7%
scale error across a 4 km frame smears the correlation peak by ~29 m -- so it
shrinks as the block converges. Place iteratively: fit the smooth part, apply it as
a continuous warp of the closed composite (never per-frame; a continuous warp
cannot open a seam), re-measure on sharper peaks, refit.

The earlier similarity placement failed for a two-character reason, not a deep
one: rotation and scale were applied with the wrong sign (`rot += th`, `gw *= s`
where the correction is the inverse), doubling the error. Fixed in close.py. The
"scale absorbs error" diagnosis in the section above was wrong.

The arterial cross-check cannot referee outside the city: the named-centreline data
stops at the Detroit limit and most cells have no arterials at all.

## The architecture the evidence supports

Two stages, coupled only at block level, with nothing in between
(`scripts/close.py`):

**Stage A — close the block.** Solve where each frame sits relative to its
neighbours from frame-to-frame observations alone, measured in map space on the
frames themselves. Modern imagery is forbidden to touch an individual frame.
Same-epoch frames were often taken seconds apart, so this is an easy, trustworthy
match — unlike film against satellite across seventy years, which is where every bad
control point in this project came from. Nothing proceeds until the block closes.

**Stage B — place the closed block.** Move the whole thing as one near-rigid body:
a handful of parameters against many well-spread observations, so no single bad
match can bend anything. A rigid transform cannot tear a seam. That is the point.

**Deleted permanently: the per-frame correction and the RBF warp.** Both are
error-concealment devices — they bought excellent absolute numbers by tearing the
mosaic. A high-degree-of-freedom warp converts measurement error into permanent
geometry. Also deleted: the mile grid as datum, since it is periodic and that is
what aliased in the first place; keep it for bias removal only.

Ship with per-frame error published. **Seams first, absolute second.**

---

# The simpler path is also the better one: COLMAP (2026-09-01, late)

The user called the hand-rolled pipeline over-engineered for the material, and
the pilot proved them right. COLMAP (brew install colmap, CLI only) on six
consecutive 1961 negatives, default settings, principal point fixed at centre,
focal prior from a 6-inch lens: 6 of 6 registered, focal solved 0.6% off the
prior, flying height 2285 m, and on the seam metric that closed the block --
close3's dense tie windows -- along-track median **2.0 m, p90 4.5, max 10**,
against 4.0 / 8.2 / 45 from the hand-rolled bundle. It models tilt and a real
camera; the hand-rolled similarity per frame cannot. Four minutes of compute.

`scripts/colmap_block.py` runs it; `scripts/colmap_render.py` projects each
negative through its camera onto the ground plane (Detroit is flat: that is
orthorectification), seam-checks the result, and composites it. COLMAP's own
model_aligner is degenerate for cameras along one line, so `--self-align` fits
the ground plane to the 3D points and a 2D similarity to the catalogue instead.
Two things to keep in mind: on flat ground a free principal point trades off
against tilt (refine_principal_point must stay 0), and the pilot had no sidelap --
the full-block run is the real test of cross-line seams.

**The full-block run failed for a reason neither diagnosis guessed.** The scans
come in five slightly different widths (5034-5088 px after cropping, same height,
so same resolution but different scan windows). COLMAP with a single shared camera
silently skips any image whose size differs from the first: 24 of 62 entered the
database, and the 8- and 10-frame "split" models were all 24 images could build.
Two runs were misdiagnosed -- first as false matches on the repeating grid, then
as a killed extractor -- before COLMAP's own log line `CAMERA_SINGLE_DIM_ERROR`
was read. The driver now pads every crop to one common canvas, centred (padding
keeps px/mm and the principal point; resizing would not), and refuses to continue
if the database is short. It also matches only catalogue-plausible pairs and fixes
the camera at the pilot's values (focal 3602 px, k -0.00027) -- both sensible,
neither was the fault.

What is kept from the hand-rolled work: the absolute measurement, the rigid /
low-order placement, the seam metric, the viewer. 1961's hand-rolled result stays
in the viewer until COLMAP's beats it on both numbers. 1967 and 1949 placed at
30 and 52 m absolute from their hand-rolled closures -- their blocks were still
partly inconsistent -- so the hand-rolled path is the fallback from here.

# 1961 is done, both axes verified (2026-09-01, late)

Closed by the per-frame tie-point bundle (`close3.py`), placed by one continuous
quadratic warp (`place2.py`), and the viewer serves it (`detroit_1961_placed3.tif`):

| | along-track seams | cross-line seams | absolute vs modern |
|---|---|---|---|
| before warp | 4.0 m (p90 8.2) | 4.5 m (p90 11.3) | 67.5 m |
| after warp | 4.0 m (p90 8.2) | 4.5 m (p90 11.7) | **4.7 m** (p90 12.2, max 25) |

The seam numbers are from close3's dense tie windows, the instrument that closed
the block. fieldapply's verifier had still been using the whole-overlap normalised
correlation that returns the search boundary on sidelaps and printed 105 m on this
same block; that is fixed and it now uses the tie windows.

# Where it stands (earlier in the day)

Every number below comes from a metric that was proved first, per block: it must
recover a planted shift on modern-vs-modern imagery (0.0 m, all four blocks), and a
shift planted on the real mosaic must move every cell's measurement by exactly that
much (0.0 m, all four blocks).

| block | frames | 3.6 x 2.0 km cells | 1.8 x 1.0 km cells |
|---|---|---|---|
| 1961 | 62 | **1.3 m** (p90 4.8, max 6.8, none over 10 m) | **3.5 m** |
| 1967 | 51 | 3.0 m (p90 9.8, max 11.6, none over 25 m) | 8.7 m |
| 1949 | 50 | 3.7 m (p90 27.8) | 18.5 m |
| 1956 | 70 | 5.0 m (p90 10.4, max 38.0) | 12.9 m |

Read the fine column. The residual varies at roughly the kilometre scale -- the
spacing of the frame centres -- so a 3.6 km cell averages most of it away, and
somebody looking at a street corner sees the local number.

**All four are independently verified** against USGS NAIP, on 1.5 km windows, on
the road response:

| block | ours | scikit-image | OpenCV ECC |
|---|---|---|---|
| 1961 | 3.84 m | 3.16 m | 3.79 m |
| 1949 | 3.47 m | 3.50 m | 5.47 m |
| 1956 | 7.74 m | -- | 12.71 m |
| 1967 | 16.22 m | 8.64 m | 11.57 m |

Note 1949, where independent measurement says 3.5-5.5 m and our own fine grid says
18.5 m. The grid uses fixed cells and many of 1949's are only partly covered
(coverage 0.63), which measures badly; the independent windows require 97%
coverage, so they test the imagery where it exists rather than where it does not.
Both are honest, and they answer different questions -- treat the grid as the
pessimistic bound.

The whole serving path is verified too, not just the files on disk: tiles pulled
from the running viewer at z17 and overlaid on modern show single yellow streets
for all four blocks, which exercises the Mercator conversion `serve.py` performs.

Downtown is a separate, weaker story. It had never been checked against anything
and was 20-42 m out, not the 3-5 m claimed. It is now corrected (1961 20.7 -> 4.5,
1949 37.4 -> 5.3) except 1956, which stays around 32 m. Downtown measurement is
unreliable: a large part of the frame is the Detroit River, which has nothing to
correlate, and the core is high-rise, so relief displacement leans every building.

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
