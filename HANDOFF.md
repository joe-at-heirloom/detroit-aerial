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

# 1961 done on the COLMAP path -- both axes, best yet (2026-09-01, 23:20)

`detroit_1961_placedC.tif` is what the viewer serves. COLMAP cameras (fixed at
the pilot's values), rendered onto the block's fitted quadratic surface,
composited, placed by one quadratic warp (held-out 3.2 m), seams verified on the
COLMAP-rendered frames before and after:

| | along-track | cross-line | absolute 3.6 km | absolute 1.8 km |
|---|---|---|---|---|
| hand-rolled tie bundle | 4.0 (p90 8) | 4.5 (p90 11) | 4.7 (p90 12) | 5.7 (p90 15) |
| **COLMAP + surface** | **2.0 (p90 4.0)** | **2.0 (p90 4.5)** | **2.6 (p90 5.8)** | **3.0 (p90 8.0)** |

Seams identical before and after the warp; 41/41 coarse cells locked.
Independently, against USGS NAIP on 1.5 km windows: scikit-image **2.24 m**
(p90 5.4, max 5.8), OpenCV **3.10 m**, ours **3.53 m** -- the previous build's
were 3.16 / 3.79 / 3.84. One command per block: `colmap_block.py TAG` then
`colmap_place.py TAG WORK --surface 2`.

**1967 on the same path -- done, in the viewer**: 51/51 frames in one model,
surface residual 3.1 m (205 m of bow at a corner), seams along-track **2.0 m (p90
4.5)**, cross-line **4.0 m (p90 7.2)**, prior 30/30 windows locked, quadratic 6.1 m
held-out, seams unchanged through the warp, absolute **2.6 m at 3.6 km (p90 7.7)**
and **4.1 m at 1.8 km (p90 11.6)**.

**1949 -- done, in the viewer** (`detroit_1949_placedC.tif`): re-aligned against
1961 (heading +1.27 deg, scale +0.26%), prior 27/27 windows locked, quadratic
6.8 m held-out, seams unchanged through the warp at 2.0 / 2.0 m; absolute against
1961 5.1 m; against modern imagery, which it was NOT placed against, **4.2 m at
3.6 km (p90 9.9)** and **4.7 m at 1.8 km (p90 11.9)**, nothing unlockable.

**1949 on the same path**: 50/50 frames in one model (seven scan sizes, padded),
surface residual 1.7 m with only 12 m of bow; seams along-track **2.0 m (p90
4.0, max 24)**, cross-line **2.0 m (p90 4.5, max 18)** over 7,629 windows. That
sidelap was 65 m out of the original bundle and 20 m at the hand-rolled best.
Its catalogue positions scatter 361 m around the solved cameras, so the absolute
placement's prior is the thing to watch there -- and it failed: the alias-proof
coarse field locked 12 of 23 windows, fell back to a zero prior, and 17 of 17
coarse cells then could not lock at all (the block sits 40-100 m off and a +/-40 m
fine search cannot reach it). The placement against 2024 imagery is therefore
unmeasured, not bad, and is NOT in the viewer. 1949 is being placed against
1961's verified build instead (`fieldmeasure.py 1949 colmap --ref 1961:placedC
--label-out colmap61`, then fit, then `fieldapply.py ... --fits colmap61 --ref
1961:placedC --out placedC2`): film to film twelve years apart is a far easier
match than film to satellite across seventy-five, and it inherits 1961's 2.6 m.
The same will apply to 1956.

**Why 1949 would not place against anything, including its own old build:** the
COLMAP block is rotated a few degrees and shifted up to ~900 m (self-aligned to
catalogue positions scattering 361 m), AND smoothly warped -- against the old build
a similarity does not fit, an affine wants -20% scale east-west only, a quadratic
leaves 93 m. Cause found in the scans: 1949 has a batch scanned square
(4998 x 4998) where the rest are 5076 x 4794 -- same negative, ~4.5% more pixels
per mm. One fixed focal for every negative makes those frames 4.5% wrong and the
bundle bends the block to fit them. `colmap_block.py` now rescales any scan of a
different resolution to the majority film height before padding (heights are the
physical film edge; widths are scan-window differences). `coarsealign.py` measures
hundreds-of-metres offsets between same-film mosaics at 20 m/px for the cases where
a block is seeded far off. Two routes ran in parallel: coarse-to-fine placement of
the warped block (fallback), and a re-solve with normalised scans as tag `1949n`.

The fallback failed and said why: the old 1949 build is itself torn 100+ m across
its lines, so window-by-window offsets against it are not a similarity, not a
quadratic, and warping to them bent the block (36 -> 80 m, seams opening). Only
two scans were square, so the resolution batch was not the main story either. The
main story is the self-alignment: seeded from a catalogue scattering 361 m, its
heading was **1.27 deg** off -- 550 m at the block's ends. `colmap_render.py
--refine-ref 1961:placedC` now refines heading, scale and shift against a
reference on the alias-proof arterial coarse field (24 of 28 windows locked, all
kept, scale +0.26%, residual 43 m -- the remaining smooth distortion the placement
warp takes) and folds it into the camera alignment before rendering; seams after
it are unchanged at 2.0 / 2.0 m, as a similarity guarantees. The whole chain takes
`--refine-ref` and persists it for the seam verification.

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

With padded crops, catalogue-plausible pairs and the fixed camera, 1961 solves as
**one model, 62 of 62 frames**, all 268 pairs verified at ~1,240 inliers each,
in 17 minutes of CPU. `colmap_place.py` then renders it, measures both seam
classes with the tie-window instrument, builds the production mosaic and places
it absolutely.

**The connected block domes.** Rendered onto one plane, its seams are along-track
16.0 m (p90 32) and cross-line 12.2 m (p90 35) -- worse than the pilot's 2.0 and
the hand-rolled 4.0/4.5. The 3D points' plane residual runs -51 to +36 m and the
camera heights follow a quadratic of -166 m across 10 km: a flat block bundled
with a slightly wrong fixed camera bows into a shallow bowl, the textbook failure
of planar aerial SfM. Two remedies, both in flight: refine focal and radial
distortion now that four lines with sidelap constrain them (mapper rerun), and
render onto the model's own fitted quadratic surface (`colmap_render.py
--surface 2`) rather than a plane -- the cameras are consistent with that
surface, so seams close, and the smooth planimetric stretch it leaves is what the
absolute warp already removes.

**Rendered onto its own surface, the COLMAP block closes to 2 m in both classes.**
`colmap_render.py --surface 2`: quadratic fitted to the 46,564 aligned points at
2.2 m residual, 144 m of bow at a corner. Full 62-frame 1961: along-track
**2.0 m median, p90 4.0, max 16** over 5,166 windows; cross-line **2.0 m, p90
4.5, max 20** over 2,193 windows; 11 windows rejected by consensus (835 on the
plane). The hand-rolled bundle's best was 4.0 / 4.5 with p90 8 / 11. Refining
focal instead was degenerate (6574 px). The whole COLMAP chain now takes
`--surface 2`; the mosaic writes the choice into data/colmap_<tag>.json and the
seam verification after the absolute warp reads it, so every step renders the
same way.

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
1949 37.4 -> 5.3, 1956 168.4 -> 2.2).

**1956 was 168 m out, and read 32 m.** That is the project's own alias trap
arriving a third time, in the one place nobody had pointed the whole-raster test
at. The fine search is capped at +/-40 m and walks ~120 m over three re-centrings;
the street lattice repeats every 97.5 m; so an epoch beyond that reach locks onto
the wrong lattice member and reports a confident small number rather than
failing. The tell was in the control all along and was read as ordinary
difficulty: 23 of 39 stations sitting at the search limit. Solved on the whole
raster with a 250 m search the peak is unambiguous (ratio 1.90) and every cell of
a 4x4 grid agrees with it. `fixdowntown.global_prior` now runs that test first on
every epoch and reports what it finds, so a scene out of reach says so.

Lesson, again: *a measurement that cannot reach the answer will return a wrong one
that looks fine.* Bound the search, then check the bound. Downtown measurement is
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

## Which build the viewer serves (2026-09-02)

The manifest once served a failed 1949 placement because its label (`placedC3`)
sorted "newer" than the verified one. Builds are now chosen explicitly, per
block, and only after seams and absolute placement are verified:

    ./.venv/bin/python scripts/manifest.py --serve 1961:placedC,1967:placedC,1956:placedC,1949:placedE

Nothing is picked by label order. Failed 1949 placements (placedC2..C5) were
deleted.

**Never regenerate the manifest without that flag.** A bare run is not a no-op and
it is not merely "unpinned": in `best()`, every `placedC*`/`placed3` candidate is
skipped unless the tag appears in the `--placed3` allow list, so with no arguments
it falls through `PREF` to `final` -- the earlier per-frame builds whose seams
were torn -- and silently swaps all four blocks for them. It changes the bboxes
too, which is the visible tell: 1956 gets noticeably wider to the west.

This happened on 2026-09-04. The manifest was regenerated bare for no reason
except to bump `ver` for cache-busting after a downtown rebuild, and every West
Detroit block was quietly downgraded off its verified COLMAP placement. A second
session spotted the bbox change and it was restored by re-running with the full
`--serve` list. Nothing downstream had been measured in between, so no numbers
are contaminated -- but nothing would have said so if they had been.

The lesson is the project's usual one wearing different clothes: *the safe default
was the one that required remembering a flag.*

**Closed the same evening** by the session adding the modern layers. `main()` now
derives the served build per tag from the files it actually chose and writes it
as `"serve": {"1949": "placedE", "1956": "placedC", ...}` on every run;
`best(tag, recorded)` takes that recorded choice first, `--serve` overrides it,
and `PREF` is consulted only when neither names a build that exists -- and then
it says so. Verified independently: the key is on disk and a bare run reprints
all four builds unchanged. A chosen or recorded build whose files are missing
now **drops the block** with `chosen build X MISSING -- block left out`; the old
walk down `PREF` to the torn `final` builds is reachable only with `--fallback`.
A hole is honest, a torn block is not. The `serve` record survives a dropped
block, so when the files return a bare run restores the right build. Pass
`--serve` when you mean to *change* a build; for anything else, a bare run is
now safe. `scripts/wipepic.py TAG LABEL` draws a build wiped against today at
mile-road crossings, the quickest visual check that the arterials line up.

## Downtown: what "warped" turned out to be (2026-09-04)

Reported: line up Campus Martius exactly on downtown 1949, to scale, and Grand
Circus (660 m away) is nowhere near. That reads as an internal warp, which no
similarity handle can fix. It is not one. Measured:

- **Ground.** 1949 vs 1961 at window sizes the metric can be trusted on downtown
  (600-1000 m, ratio >= 1.3, 65-82 windows): median 2.0 m, p90 8.5, max 17.5.
  Campus Martius +8.9 E, Grand Circus +4.6 E -- 4.5 m apart over 660 m. The
  uncorrected `_src` raster was a uniform ~40 m east (a placement, not a warp)
  and the fix removed it cleanly. City street-centreline vectors sit on the
  streets in 1949, 1961 and Esri Today at both landmarks to a few metres.
- **Roofs.** Image-to-image correlation in the high-rise core is not measuring
  ground. 1961 vs Esri Today at Campus Martius: -24.7 E at ratio 2.2-3.8, the
  strongest lock in the exercise -- and false, because the vectors say the ground
  agrees. Esri vs NAIP there: +13 to +24 E; 1961 vs NAIP: 0-3 m. The
  1949/1961/Esri triangle does not close (0.5, 7.9, -24.7), which is what pairs
  locking on different leaning features look like. The lean is ~25-35 m at the
  core and ~5 m by the river.
- **The metric below ~500 m downtown is unreliable.** The same 400 m window at
  Campus Martius returned +33 E on one run and -33 E on the next, ratio ~1.1; a
  window-size sweep flipped between ~40 m and ~10 m. Irregular blocks plus
  high-rise. Use >= 600 m and demand ratio >= 1.3, or use the vectors.

The concrete case: Joe had lined 1949 up on the First National Bank Building
(26 storeys, 104 m, 1921 -- the stepped "Z" footprint south-east of the park).
Crop it at 0.42 m/px with the street vectors over 1949, 1961 and Esri: the
streets sit on the vectors in all three, Esri's roof sits over the OSM footprint,
and the *same roof* in 1949 and 1961 sits ~35-40 m west/north-west of it. Align
that roof onto Esri's and the whole 1949 layer moves ~35-40 m -- which is what
Grand Circus then showed.

So: lining Campus Martius up "perfectly" by eye means matching building outlines,
the dominant feature in a high-rise core, and those lean 25+ m in Esri and
differently again in 1949. That puts the ground 25 m off, and the park at Grand
Circus shows it. The hand-alignment handle is the wrong tool for the core; the
right reference for aligning by eye downtown is street geometry, not buildings.

Two dead ends, recorded so nobody repeats them: correlating with the same street
mask on both sides returns exactly zero at ratio 1e9 (the mask correlates with
itself); and correlating ridge-filtered imagery against a thin rasterised
vector line read ~30 m east for every image including NAIP, while the overlay
of that same raster visibly sat on the roads -- a bias of the method, not the
map. The vector check that worked was the picture, not the number.

## 1949 against USGS NAIP (2026-09-02)

`crosscheck.py --block 1949 --suffix placedC --n 12 --span 1500 --mpp 1.0 --bound 60 --ridge`:
ours 8.39 m (p90 9.85), scikit-image 8.03 (p90 10.66), OpenCV 8.20 (p90 9.70);
bias under 1 m both axes. Random scatter, not a shift; the leave-one-out figure
of the placement fit (6.8 m) predicted most of it. The candidates from the
normalised scans (1949n) are the only open route to improving it; switch only
if seams and absolute both measure better.

## 1949 from the normalised scans (1949n): not better, discarded (2026-09-02)

Same 50 negatives after tone normalisation, solved and placed by the same
chain (refine-ref 1961:placedC): seams 2.0 / 2.0 m (identical), leave-one-out
7.2 m (vs 6.8), absolute vs 1961 5.2 m at 3.6 km, 6.4 m at 1.8 km (vs 5.1).
Normalising the scans changes nothing COLMAP cares about. The served 1949 stays.
Solve kept at /tmp/colmap_1949n for the session; mosaics deleted.

## 1956: a 33 km block needs a cubic placement (2026-09-02)

COLMAP registered all 70 negatives (four lines, 33 km long); seams 2.0 / 4.0 m
(p90 8.0 / 8.9). The catalogue placement was 815 m off, so `refine_against`
now iterates and widens its search when fewer than half the windows lock
(600 m -> 1500 m); three passes converged, rotation -0.69 deg.

Placement against 1961 with the quadratic model: 12.5 m median, p90 46 m. The
raw field is a parabola along the flight lines (-50 m in the middle, +80 to
+110 m at both ends) plus a cross-track scale gradient. The solved camera
heights bow by 250 m over each line -- the usual fixed-focal dome -- and the
rescaled square scans sit exactly on that curve, so the rescale is right. A
`cubic` model (20 params) added to fieldfit wins by leave-one-out: 7.8 m.

Residuals by region: inside 1961's coverage 5.6 m median (p90 14); the western
column and the northern 3 km lie outside 1961 and are measured against modern
imagery, where 1956 matches as badly as 1949 does: 18 m median with outliers
over 100 m that are failed matches, not block errors. No epoch covers those
strips, so they stay unverified beyond that.
Applied (cubic): seams unchanged 2.0 / 4.0 m; against modern 5.6 m median at
3.6 km (p90 21.6), 6.2 m at 1.8 km (p90 38.3). Served as `1956:placedC`.
Against USGS NAIP (`crosscheck.py --block 1956 --suffix placedC --n 12 --span 1500
--mpp 1.0 --bound 60 --ridge`): ours 5.90 m (p90 31.2), scikit-image 4.47 (p90 29.3),
OpenCV 4.59 (p90 21.2); bias dE -2.6 dN +1.0. Agrees with our own 5.6 m / p90 21.6.

## 1949 re-placed: direct to modern, cubic (2026-09-02)

"1949 is slightly off kilter from modern map." It was not rotated -- five
independent instruments put the block rotation at 0.01-0.09 deg against a 1961
control in the same range. What was real was a LOCAL bearing warp varying with
latitude (about -0.24 deg at Joy and Plymouth, nil at Fenkell and McNichols),
which no rigid measurement sees.

Fixed by two changes: fit the placement with the `cubic` model (added for 1956,
after 1949 had already been placed), and measure the field against modern
imagery directly instead of chaining through 1961:

    fieldmeasure.py 1949 colmap --label-out colmapM     # vs modern, not 1961
    fieldfit.py 1949 colmapM --use cells                # cubic wins, LOO 4.9 m
    fieldapply.py 1949 colmap --fits colmapM --out placedE

Independent arbiter (USGS NAIP, which nothing is fitted to), ours/skimage/opencv:

    placedC (was served)  8.39 / 8.03 / 8.20 m   p90 9.7-10.7   max 13.2-14.9
    placedD (cubic vs 1961) 4.46 / 4.04 / 4.30   p90 8.3-9.2    max 10.3-12.4
    placedE (cubic vs modern) 4.44 / 5.05 / 5.32 p90 6.6-7.1    max 7.3-7.5

Seams unchanged at 2.0/2.0 m through the warp. Served as `1949:placedE`;
placedD deleted. The same "fit direct to modern with a cubic" question is open
for 1956, which is still chained through 1961.

## Extending a block: what worked and what did not (2026-09-04)

Wayne County holds exactly four years, so "more years" for Detroit is not
available in this collection; more coverage in those four years is. The whole
1949 roll ha-17 was fetched (118 frames, against 50 in the served build) and
solved.

**What worked.**
- `fetchscans.py` at native resolution. A fixed width fails with HTTP 501 on any
  negative scanned narrower than it (IIIF level1 will not upscale) and these run
  5352-5400 wide, so the original store was silently missing frames.
- `rollplan.py`: consecutive exposure numbers are consecutive shots along a line
  and overlap ~60%, which pairs them with no position at all; within a run,
  position is near-linear in number (1200-1400 m a frame), so frames a few
  numbers off a solved run extrapolate to a few hundred metres.
- Bootstrapping. The first solve registered 84 frames but the new ones had only
  along-track pairs, leaving the four lines wanting -24/-11/+3/+24 m in E -- a
  0.9% cross-block scale error. Re-pairing from the solved positions took
  cross-line pairs 250 -> 352, cross seams p90 6.3 -> 4.5 m, and placement
  7.1 -> 5.1 m.

**What did not.** The 92-frame build is still not better than the served 50-frame
one on the same ground: median 4.3 m against 4.2, but p90 13.4 against 7.1. The
new ground it reaches (south to lat 42.06) is only 2-3 flight lines wide, so it
carries almost no measurable control -- 2 usable cells at 3.6 km, none at 1.8 km.
Lengthening a narrow block does not pay; widening it would.

**Widening is blocked on locating rolls.** ha-16 and ha-18 were fetched (172
frames) and locate nowhere on the block: rolls cover different regions, not
adjacent strips. `locate.py` places a raw negative by FFT cross-correlation
against a placed mosaic and works well there (3 of 4 known frames within 45 m,
the fourth correctly refused), but against modern imagery over the whole county
it fails outright -- 75 years of change plus a grid that repeats every 97.5 m.
The route that should work is matching to the ROAD NETWORK rather than to
imagery, which is how the original catalogue was built; `data/overlay_segs.npy`
covers Detroit city (30 x 49 km) but not Downriver, so it would need extending
from OSM first.

Served builds are unchanged. `1949x` (84 frames) and `1949y` (92 frames) are kept
as evidence and are not in the manifest.

## More years: where they can come from (2026-09-04, evening)

The DTE catalogue online holds Wayne County in exactly 1949 / 1956 / 1961 / 1967
(547 / 432 / 533 / 389 frames). The physical collection also has 1952, 1981 and
1997, but online those are Oakland 1952 (241 frames, rolls de-21..de-32), Monroe
1981, St Clair 1985 and Oakland + Macomb 1997 -- none over Detroit.

Public sources verified today over the West Detroit bbox (lat 42.19-42.49, lon
-83.38 to -83.12), all already orthorectified, none yet measured by validate.py:

| year | source | res | endpoint |
|---|---|---|---|
| 1951 | Michigan Tech tile layer of USGS EarthExplorer frames, spline-warped (seams visible) | z18 | tiles.arcgis.com/tiles/RPhrOu9XQzI31xTa/arcgis/rest/services/Aerial_Imagery_of_the_City_of_Detroit_1951/MapServer |
| 1998 | City of Detroit DOQ | 1 m | egis.detroitmi.gov/image/rest/services/Imagery/1998_Aerial_Imagery/ImageServer |
| 2005, 2010 | City of Detroit | 2 ft | egis.detroitmi.gov/image/rest/services/Imagery/{2005,2010}_Aerial_Imagery/ImageServer |
| 2012, 2014 / 2016-2022 | NAIP via Planetary Computer STAC (34 tiles a year over the bbox) | 1 m / 0.6 m | planetarycomputer.microsoft.com/api/stac/v1, collection naip |
| 2014-2026 | Esri World Imagery Wayback, 196 releases | 0.3 m | wayback.maptiles.arcgis.com |
| 2020, 2024 | City of Detroit MiSAIL tile caches | 6 in | tiles.arcgis.com/tiles/qvkbeam7Wirps6zC/arcgis/rest/services/{MiSAIL_2020_6in_Clip_webMerc,2024Sp_Wayne_6in_MiSAIL_tileCache}/MapServer |

The 1951 layer matters for what it proves, not what it is: USGS holds a 1951 flight
over the city, and the raw frames are on EarthExplorer (login required) -- the
right input for colmap_block, where the tile layer is a rubber-sheet.

North of 8 Mile (the blocks reach lat 42.475-42.486, up to 4 km into Oakland)
Oakland County serves public ImageServers for 1940 (1.4 m), 1949 (2 m), 1963,
1974, 1980, 1990, 1997, 2000, 2002, 2005, 2006, 2008, 2010, 2012, 2014, 2015,
2017, 2020 (0.25 m), 2023 and 2025 at
gisservices.oakgov.com/arcgis/rest/services/ImageServices/EnterpriseOrtho{BW,TC}<year>ImageService/ImageServer;
exportImage verified over Southfield for 1940-2000.

Not online: SEMCOG's flights (1966, then every five years 1970-2020; request
only), MSU RS&GIS archive (1930s-2000s statewide, $30 a frame, "roughly
referenced"), U-M Clark Library prints (Wayne 1963/64/69/80/90/2000), Wayne
County's own GIS servers (not publicly reachable). The IU "Historic Wayne County
Images" set (1936-1994) is Wayne County, Indiana.

Oakland 1952 negatives, fetched at 1000 px and run through locate.py against
1956:placedC: 9 of 241 locked. Four are credible -- de-31-30/31 at lat 42.449/42.462 and
de-31-105/106 at 42.454/42.443, consecutive numbers landing a frame apart -- so a
few 1952 frames straddle 8 Mile at the 1956 block's north end; the other five sit
deep inside Detroit at ratios of 1.27-1.45 and are the grid aliasing. Oakland
County already serves 1949 and 1963 orthos over that strip, so 1952 raw film
buys one year over ~2 frames' width. Not worth a solve on its own.

## Modern layers in the viewer: 1998, 2005, 2010, NAIP 2012-2022 (2026-09-04, night)

`scripts/fetchlayer.py` pulls an already-orthorectified layer onto the project
grid as a 3-band GeoTIFF (`mosaics/layer_<name>.tif` + `data/layer_<name>_geo.json`);
`manifest.py` slots every such layer between the film and Today in year order;
`serve.py` reads three bands when a GeoTIFF has them. Sources, all public:

| layer | source | native | how |
|---|---|---|---|
| 1998, 2005, 2010 | City of Detroit ImageServers (egis.detroitmi.gov, Imagery/<year>_Aerial_Imagery) | 1 m, 2 ft, 2 ft | exportImage in EPSG:4326 |
| 2012, 2014, 2016, 2018, 2020 | State of Michigan NAIP ImageServers (imagery.michigan.gov/server/rest/services/Michigan_NAIP_<year>) | 1 m, 1 m, 0.6 m x3 | exportImage in EPSG:4326 |
| 2022 | USGS NAIPPlus ImageServer (imagery.nationalmap.gov) with `--where Year=2022` | 0.6 m | exportImage in EPSG:4326 |

Two things learned on the way. An ImageServer asked for EPSG:4326 output honours
the bbox only if the requested pixel grid has the bbox's aspect in DEGREES, so
chunks are requested at square-degree pixels and warped onto the square-metre
grid afterwards; the returned extent is asserted against the request every time
(ruled-out #7 again). And Microsoft Planetary Computer's anonymous NAIP access is
throttled: 12 MB/s for the first half-gigabyte, then 0.2-1 MB/s for good, which
put a 0.6 m year at 17 GB and most of a day. GDAL range reads through it were
ten times slower per byte again. The State of Michigan serves the same NAIP at
30 MB/s, so `--naip` remains in the script but is not the route.

The city's 1998 ortho, measured by validate.py against the Esri reference exactly
as our own builds are: **2.3 m** median at 3.6 km (p90 3.2, max 4.8, 23 of 23
locked) and **3.0 m** at 1.8 km (p90 5.6, max 8.5). Every layer, measured the same way (median at 3.6 km / at 1.8 km, p90 at 1.8 km):

| layer | 3.6 km | 1.8 km | p90 | note |
|---|---|---|---|---|
| 1998 | 2.3 | 3.0 | 5.6 | city only; 23 of 48 coarse cells have coverage |
| 2005 | 3.7 | 3.6 | 4.7 | one 35 m cell at the layer's north edge is partial coverage, windows 1 km away lock at 2 m |
| 2010 | 2.2 | 2.6 | 3.5 | city only |
| 2012 | 3.7 | 3.9 | 4.9 | |
| 2014 | 3.8 | 4.1 | 5.3 | |
| 2016 | 1.5 | 1.6 | 2.5 | best of the set |
| 2018 | 3.0 | 3.4 | 9.6 | 18 fine cells at 10-14 m, all in the eastern two columns, all +10 to +14 m east: a shifted quarter-quad column in the NAIP 2018 product |
| 2020 | 3.6 | 3.7 | 6.6 | |
| 2022 | 2.0 | 2.1 | 3.9 | |

Nothing was corrected: these are served as published, and the table says which
wipes will show a few metres of movement that belong to the source, not to the
film. The 2018 column offset is the one a viewer will notice, east of Greenfield.

One more guard, added the same night: manifest.py records the build it served
per block (`serve` in data/manifest.json) and a run without `--serve` keeps
that choice, carried across runs even when a block's files are missing. A
missing build now leaves the block OUT of the viewer with a loud line rather
than falling through PREF to the torn `final` build; `--fallback` restores the
old order on purpose. A bare run at 20:19 had silently downgraded all four
blocks to `final`; it was caught only because a layer fetch printed a
different block extent.

## Alignment plan, steps 1-4: baseline, the lens question, and the seam metric's floor (2026-09-27)

Working from ALIGNMENT_PLAN.md. Nothing served has changed.

**Baseline.** `runs/baseline_20260923/` holds sha256 of every served raster and
the Esri reference, the manifest, and both versions of `data/adjust.json`. Every
west GeoTIFF's own transform equals its manifest bbox exactly; the three downtown
TIFFs differ by 1.2 m (south) and 0.75 m (east), which the viewer stretches away.

**Downtown 1949's hand alignment makes it ten times worse.** `dtcross.py` had
crashed on every run since hand alignment was added (the loop that prints the
offsets rebound the argparse `a`), so neither saved alignment was ever measured.
Fixed and measured (1.25 m/px, 5x5), 1949 against 1961 / 1956:

    unadjusted                              3.2 / 4.0 m
    committed  +41.2 E  +2.2 N  x0.95       35.2 / 43.3
    worktree   +20.6 E +17.0 N  x0.95       44.0 / 38.7

Roof lean again (see "Downtown: what 'warped' turned out to be"). Removed on
Joe's say-so the same day: `data/adjust.json` is now empty, `dist/tiles/dt1949`
was deleted and rebuilt (build_static skips tiles already on disk, so the old
ones must go first; `--layers dt1949` would also rewrite dist's manifest with
that one layer, so the default film build was run). 1077 tiles, against 1006
before -- the x0.95 had left edge tiles empty. Measured on the SERVED z17 tiles
against dt1961's: 1.6 m median (p90 5.9), where the old tiles read 89 m.

**`scripts/crossmatrix.py`**: every layer of a manifest group against every
other, on one lat/lon grid. Each layer is read once by windowed GDAL reads with
area averaging, hand alignment applied through `handadjust`, ridge-filtered at
gridval's two bands and cached in `runs/cache/` by source signature. Each pair
gets gridval's 6 km coarse field as prior, then fine cells; each cell is
re-measured with a planted field (30 m + 1.5 m/km scale + 1 mrad rotation) and a
cell that misses it by more than 5 m is dropped. `scripts/xmreport.py` makes the
table and an arrow map per layer. Smoke test: 1961 vs Esri 2.8 m median over 1.3
km cells (HANDOFF's 2.6 at 3.6 km), planted field tracked to 0.8 m (p90 1.1).
Downtown film-to-film: 1949/1956 3.7, 1949/1961 2.9, 1956/1961 2.4 m. Downtown
against Today does NOT measure: most cells fail their planted field, and four
roof locks at 30-160 m pass it -- a wrong lock can be self-consistent. Use the
street vectors there.

**The seam metric is whole-pixel.** `close2._pc` returns the integer peak, so at
2 m/px every tie window reads 0, 2, 2.8, 4... m, and "2.0 m (p90 4.0)" has meant
"within one pixel". With a parabolic peak (`distortion_pilot.py --subpixel`)
1961's rebuilt solve is **along-track 1.7 m (p90 4.0), cross-line 2.1 m (p90
5.0)**. Any acceptance test on seams below 2 m needs the sub-pixel correlator.

**1961 rebuilt** in `runs/colmap/1961` (not /tmp): 62/62 frames, one model,
seams identical to the served build's (integer metric 2.0/2.0, p90 4.0/4.5).
The 1956, 1961 and 1967 scans had been lost with /tmp; 1961's were re-fetched
(four of its frames are missing from dte_catalogue.json and were fetched by
pointer). 1956 and 1967 scans are still absent; 1949 has 48 of 50.

**Is there lens or scan distortion? Only a metre of it, and fixing it buys
nothing.** `scripts/distortion_pilot.py` renders every frame through its solved
camera and asks two questions:

- seams, no reference: a radial error dk shows as disagreement that varies
  across an overlap as dk(|ra|^2 ra - |rb|^2 rb)/H^2; per-pair translation is a
  nuisance. 5-fold by pair, held-out within-pair scatter, sub-pixel ties:

      fixed camera                        1.15 m (p90 3.30)
      CONTROL: r^3 about wrong centres    1.14   (3.31)
      radial r^3                          1.06   (3.16)   dk -1.8e-3
      radial r^3 + r^5                    1.02   (3.17)   coefficients collinear
      scan affinity per batch             1.14   (3.29)   nothing

  It beats the control, but it is not stable: dk is -2.1e-3 north, -1.6e-3
  south, -2.2e-3 from along-track ties and -0.3e-3 from cross-line ties, and
  leave-one-pair-out the gain is +0.013 m per pair (95% CI -0.008..+0.033; 118
  of 201 pairs improve).
- absolute, against NAIP 2016, per-frame similarity removed, 4,201 windows
  stacked in image coordinates: residual 4.89 m median (p90 11.96); with the
  shared cubic 4.87 (11.92), dk -4.2e-3. The radial profile has the cubic's
  shape (+0.4..+0.6 m inside 1 km, -1.3 m at 1.65 km), tangential under 0.15 m.

So a real radial term of about a metre at the frame edge exists, but its
estimate disagrees between tests by a factor of two and more, it is invisible
in the held-out seams and against modern imagery, and the plan's targets are
3 m median / 10 m p90. By the plan's rule the fixed SIMPLE_RADIAL camera stays;
raw-frame recalibration is not justified, which also means 1956 and 1967 do not
need re-fetching for it. What the 4.9 m within-frame scatter against NAIP IS
remains open: fifty-five years of change, 520 m window noise, and the fitted
ground surface are the candidates, and the pair matrix is the next evidence.

Running long jobs on this machine. CPU SIFT peaks at ~3.7 GB per thread on these
negatives, so the old fixed 8 threads is ~30 GB and pages a 32 GB machine to a
standstill; `colmap_block.py --threads 2` (added today) stays in RAM at ~40 s per
image per thread. zsh runs `cmd &` at nice 5, which starved COLMAP to half a core;
launch with `setopt NO_BG_NICE` and check `ps -o ni`. Killing crossmatrix.py's
parent leaves its pool workers running (`pgrep -fl multiprocessing`).

## The west pair matrix (2026-09-27)

`crossmatrix.py west`: 15 layers, 105 pairs, 1.3 km cells at 2 m/px, every
pair's planted field tracked to 0.8 m (p90 1.1). Report `runs/crossmatrix_west.md`,
arrow maps `runs/crossmatrix_west_vs_l2016.png` and `_vs_esrihi.png`. Medians
against NAIP 2016 / Esri hi:

    1961 3.4 / 2.8   1967 4.3 / 4.0   1949 5.5 / 4.4   1956 7.2 / 7.2
    1998 3.8   2005 3.1   2010 2.0   2012 4.0   2014 3.7   2018 3.0 (p90 17)
    2020 5.1   2022 2.1   Esri hi 1.9   Today (modern_west.png) 12.2

Film against film, what a wipe shows: 1949/1956 4.6, 1956/1961 4.8,
1961/1967 4.7, 1949/1961 6.0, 1949/1967 7.1, 1956/1967 8.4 m.

- **The local viewer's west Today is 12 m out.** Every layer, film and modern,
  sits 11-15 m north of `modern_west.png` with the same vector in every cell;
  the same layers agree with Esri hi to 2-4 m with no bias. That PNG is the old
  Web Mercator image (section 2 above) and serve.py still offers it as Today.
  The static site is unaffected: its Today is live Esri tiles.
- **2016 holds as the reference**: 1.9 m from Esri hi, 2.0 from 2010, 2.1 from
  2022. 2018's eastern column is its documented +10-14 m E quarter-quad; the
  rest of 2018 agrees.
- **1956's north-west is unverified, not 100 m out.** Its cells at -95 to -125 m
  E repeat identically against seven modern years and pass the planted field --
  but at full resolution the ground is 1956 farmland under 2016 subdivisions,
  and the one shared feature, the main east-west road, differs by ~20 m N-S.
  A consistent false lock on changed ground passes every self-check this
  instrument has; there, only roads that existed in both years can arbitrate.

## Step 5: joint placement of the film years (2026-09-27)

`jointfit.py` fits a correction field per film year to the pair matrix at once --
film vs NAIP 2016 (held fixed) AND film vs film -- and scores each model order on
cells withheld by 4 km block. Esri hi is never fitted and never used to choose.
`jointapply.py` composes the chosen correction with the served placement and
re-warps the COLMAP mosaic once (the saved placement reproduces the served raster:
patch correlation 0.93-0.997), after checking the field's size, gradient and
foldover; beyond the measured ground its coordinates are clamped, not extrapolated.

What the evidence supports:

- **1967: quadratic, a clear win.** Candidate `1967:placedJ` (NOT served), field
  under 16 m, gradient 0.35%/km (seam change < 1 cm by construction; the frames
  cannot be re-rendered, the solve and scans are gone). Cell-matched against
  every layer it was never fitted to, served -> candidate median (p90): 1998
  5.2 -> 3.9 (11.0 -> 8.7), 2005 5.7 -> 4.4, 2010 5.0 -> 4.0, 2012 6.3 -> 5.2,
  2014 6.1 -> 4.3, 2020 5.2 -> 4.2, 2022 5.4 -> 3.6, Esri hi 4.0 -> 3.2 (8.7 ->
  6.7). Wipes: vs 1961 4.7 -> 3.8, vs 1949 6.8 -> 6.4. wipepic.py shows no
  artefact at the six crossings.
- **1956: nothing earns it.** Unconstrained, a cubic scored 6.0 -> 3.7 on Esri
  hi -- by swinging ~200 m through the north-west, where 47 cells were gated as
  changed ground and nothing holds it. With a weak 'no change' prior on
  unmeasured ground (`--damp`) every order does WORSE than none on the withheld
  blocks (6.1 -> 6.4); the cubic also fails the gradient check (1.03%/km).
  1956's error lives where the ground changed; a smooth field cannot reach it.
- **1949 and 1961: no correction.** 1949's narrow block cannot carry more than a
  shift and a shift buys nothing (4.38 -> 4.31 on Esri hi); 1961's cubic is 2.79
  -> 2.61, inside the instrument's 0.8 m.
- The joint fit is worth doing: film-to-film agreement in-sample is 4.2 m joint
  against 4.8 m fitting each year to the reference alone.

Also 2026-09-27: the local viewer's west Today is now `modern_west_hi.tif`
(manifest.py prefers it); served Today vs 2016 at four spots 0.3-4.5 m in
scattered directions, where the old PNG sat a uniform 12 m south.

## Every year we could find, and a view that stacks them (2026-10-03)

Joe asked for more years and for a view with a year slider: the newest picture
taken in or before the chosen year on top, older ones wherever they reach
further, each outlined thinly with its year. Both are done. The viewer's first
tab, **Through the years**, now has 41 steps from 1931 to Today.

### The view

`scripts/footprints.py` reads every layer at ~15 m cells (any band non-zero, which
is exactly what serve.py draws as opaque), cleans and vectorises the coverage, and
for each year computes newest-first which part of which layer no newer layer
covers. It writes `data/timeline.json` (per layer: year, label, credit,
footprint; per year: the stack with tile bounds, and one outline + caption per
picture). Run it after manifest.py; restart serve.py after both. The viewer loads
only tiles inside each visible region, merges a layer with its `part_of` strips
into one outline, places captions so they do not collide, and names every
picture's publisher in the Source cell. Clicking the map opens the probe on
every layer whose outline holds the point.

### What was added, and how well it sits

validate.py against the Esri reference, exactly as every other layer: median
error on 3.6 km cells / on 1.8 km cells (p90), metres. The small partial layers
(Taylor, Oakland) get the same cell counts over a much smaller box, so their
"1.8 km" cells are a few hundred metres and noisier; read the first number.

| year | layer | source | coverage | 3.6 km | 1.8 km (p90) |
|---|---|---|---|---|---|
| 1931 | erca1931 | Essex Region Conservation Authority, Detroit River 1:10,000, rubber-sheeted | river corridor, downtown to Grosse Ile, both banks | not measurable: the reference stops short of the river (1 cell) | |
| 1940 | oak1940, tay1940 | Oakland County; City of Taylor | north of 8 Mile; Taylor | 8.0 (Taylor) | |
| 1949 | oak1949 | Oakland County (DTE prints, EDCA rubber-sheet) | north of 8 Mile, under our film | | |
| 1951 | 1951 | Michigan Tech tile layer of USGS frames (4-20-'51), spline | City of Detroit + downtown | 3.0 | 3.5 (7.9) |
| 1957 | tay1957 | City of Taylor, published as **1964**: frame stamps read 5-16-57 (USDA XU) | Taylor | 6.1 | |
| 1963-1997 | oak1963/74/80/90/97 | Oakland County | north of 8 Mile | 4.4 / 3.8 / - / - / 3.7 | |
| 1972, 1985 | tay1972, tay1985 | City of Taylor | Taylor | 5.7 (1972) | |
| 1983 | nhap1983 | USDA FPAC orthorectified NHAP CIR, 5/10 May 1983 | whole block + downtown | 4.2 | 4.9 (8.0) |
| 1987 | nhap1987 | USDA FPAC orthorectified NHAP2 CIR, 14/17 Jun 1987 | whole block + downtown | 2.8 | 3.1 (6.4) |
| 1999 | 1998 (relabelled), napp1998 | City + State NAPP, **published as 1998**: Hudson's is a cleared lot (imploded Oct 1998), EE dates the DOQs 1999-03-28 | city; whole block (CIR) | 3.1 (NAPP) | 3.6 (6.3) |
| 2000 | oak2000, tay2000 | Oakland County; Taylor (frame label MAR 26, 2000) | strips | | |
| 2002 | 2002 | USGS HRO 1 ft, 10 Apr-4 May 2002, served by Oakland County | whole block | 2.6 | 2.6 (6.6) |
| 2004 | tay2004 | City of Taylor (needs `--native`: blank when asked for EPSG:4326) | Taylor | 3.4 | |
| 2005 | naip2005 | NAIP (State mirror), summer; under the City's spring 2005 | whole block | 4.7 | 4.8 (6.5) |
| 2006 | naip2006, oak2006 | NAIP 2 m, gappy; Oakland B&W | ~40 %; strip | 6.1 | 6.6 (9.3) |
| 2008 | 2008, oak2008 | NOAA Great Lakes border ortho 0.3 m (31 Jul-1 Aug 2008); Oakland (16-bit, server stretch) | whole block + downtown; strip | 3.5; 4.8 | 3.6 (5.0) |
| 2009 | 2009 | NAIP (State mirror) | whole block | 2.8 | 3.1 (5.2) |
| 2010 | naip2010 | NAIP, summer; under the City's spring 2010 | whole block | 2.7 | 3.0 (5.0) |
| 2015 | 2015, oak2015 | Wayne County / SEMCOG spring 2015 tiles (white outside Wayne, masked); Oakland | Wayne; strip | 3.9 | 4.1 (5.7) |
| 2017 | oak2017 | Oakland County | strip | 2.6 | |
| 2019, 2023 | wb2019, wb2023 | Esri Wayback (Maxar satellite), **remote tiles, never downloaded** | everywhere, but only from z14 | | |
| 2021 | tay2021 | City of Taylor | Taylor | 1.4 | |
| 2023, 2025 | oak2023, oak2025 | Oakland County (white nodata, masked) | strip | -; 1.6 | |
| 2024 | 2024 | State MiSAIL spring 2024 (11 Mar Wayne), public tile cache, identical to the 6 in ImageServer | whole block | 2.4 | 2.7 (4.2) |
| 2025, 2026 | noaa2025, noaa2026 | NOAA NGS DSS 0.25 m (10 May 2025; 28-30 May 2026) | SE corner + downtown; river strip | 3.0 (2025) | 3.0 (3.8) |

A dash is a layer not run through validate.py (Oakland 1949/1980/1990/2000/2006/
2015/2023, Taylor 1985/2000, NOAA 2026); the rest of each family measured
1.4-8 m, and Taylor's film years (1940, 1957, 1972: 5.7-8.0 m with outliers) are
the weakest -- municipal rubber-sheets of single frames, a few metres worse than
our own solves. NAIP 2005 and 2009 were measured before their band order was
corrected; validate.py works on the mean of the three bands, which for 2009 is
the same and for 2005 nearly so.

Plus `<name>e` strips (`--part-of`) for every full-coverage year from 1999 to 2024:
the West block's layers stop at lon -83.115, so without them the year view showed
1951 film between Woodward and downtown in every year up to 2006.

Fixes made on the way, each one a thing the year view exposed:
- The City's 2005 and 2010 orthos fill outside coverage with **white**, which
  hid every older year (17 % of 2010 was white, all of Dearborn). `maskfill.py`
  turns large white regions into nodata; 2015, Taylor 2004/2021, the City east
  strips and Oakland 2023/2025 needed it too.
- Taylor's "1964" is 1957 and the City's "1998" is 1999 (above). Labels now say
  what the ground says; `published_as` in the geo JSON keeps the publisher's year.
- Wayback releases show their named capture only from zoom 16; below that Esri
  serves other dates. The remote layers draw z16 tiles from z14 and nothing
  further out, and the year view says "exists only close up" instead of showing
  an empty outline.
- Oakland 2008 is 16-bit; pixelType=U8 saturates it white. fetchlayer.py asks the
  server for a 0.5 % percent-clip stretch on any non-U8 service.
- Taylor's small server rate-limits tiles (403 after ~1000); get_tile backs off.
- The State's NAIP services do not share a band order: 2005 is NIR,R,G,B and
  2009 is B,G,R,NIR (2010-2020 and USGS 2022 are R,G,B,NIR). Reading bands 1-3
  blind made 2005 false-colour and swapped 2009's red and blue; the yearsheet
  showed it. fetchlayer.py now correlates a raw export against the service's own
  PNG rendering at the extent's centre and requests the matching bandIds
  (`--bands` overrides). Every other layer was checked the same way and is right.

A pre-commit review (four reviewers, each finding then attacked by an independent
skeptic; 13 of 18 held) found, and this commit fixes:
- maskfill's first version left the white fill's JPEG noise as a speckled sheet
  (2010 over Dearborn) and white slivers at its edges; it now takes near-white
  and neutral (>= 235, bands within 10), counts nodata as fill so it can be
  rerun, opens away anything under 5 px inside the fill, and never touches a
  white region that does not reach nodata or the raster edge (a big roof). Rerun
  on every layer that had fill, plus NAPP 1999 and Oakland 2015.
- The 16-bit stretch mapped the darkest picture to 0 = nodata (holes in Oakland
  2008); it now outputs 1..255. Oakland 2008 refetched.
- Wayback layers draw everywhere inside their tile bounds, but the timeline only
  counted two cover boxes, so the gap between them showed 2019 under a 2008
  caption. `cover` is now one rectangle, snapped to the z16 tile edges they draw
  from, and the viewer stops their tiles at that edge. They are no longer wipe
  chips (black at block zoom); the probe marks them "close up only" below z14.
- Probe links now carry their tab (#lat,lon,z,tab); a link without one opens the
  scene whose box holds it, never the year view, so old West links still open West.
- Smaller: centrelines ticked in the year view drew under its layers; the
  close-up notice judged the unrounded zoom; yearsheet crashed downtown (the
  Today JPEG has no georeferencing).

`scripts/yearsheet.py LAT LON` is the quick look: one spot, every year in the
timeline, composited exactly as the year view draws it, captioned with year and
publisher (logs/yearsheet_northland.jpg: 8 Mile & Greenfield, farmland in 1940,
Northland Center from 1954, cleared by 2023).

fetchlayer.py grew three source kinds: `--tiles` (any Web Mercator cache),
`--cogs` (a NOAA tile index or URL list, read at the overview nearest --mpp), and
MapServer `/export` (with `--native` for services that cannot reproject), plus
`--clip`, `--label`, `--prio`, `--credit` and `--part-of`.

### Found, not fetched

- **Duplicates of years already served** (skipped): SEMCOG 2005 regional (same
  flight as the City's 2005), SEMCOG 2010 Wayne (pixel-identical to the City's
  2010), Wayne County 2020 6 in and 2024 LERC, Detroit MiSAIL 2020/2024 caches,
  NAIP 2024 (FPAC), Ontario SWOOP 2015/2020/2025 (US bank only), Livonia/Ferndale/
  Allen Park/Taylor clips of county flights, Wayback 2012/2015/2020/2022. Wayne
  2020 6 in (leaf-off) would be a sharper 2020 than NAIP if wanted.
- **Licensed, not ours to republish**: Dearborn Nearmap Apr 2026, DTW airport
  Nearmap Jul 2025.
- **Too small or not usable**: ERCA 1947 (1.4 km2 of US bank), ERCA 1988 (0.7 km2),
  SEMCOG "Then & Now 1950" frames (stamped Apr-May **1949**, collars on), a 1981
  Woodbridge clip (no provenance), NARA 1937/1940 USDA photo-index sheets (index
  mosaics, 2-8 m, would need hand georeferencing), The Henry Ford's 1925 Ford
  Airport photomosaic (uncontrolled, Dearborn only).
- **Needs Joe** (login, request or purchase, all raw frames for colmap_block):
  EarthExplorer single frames 1951, 1956, 1966-68, 1973, 1978, 1980, NHAP 1983/87
  frames, NAPP 1993/94/99, HRO 2004/2008/2010/2011 (login); USAF/USAAF film at NARA
  RG 373 for 1942, 1944 (whole block), 1952 (whole block), 1953 (whole block),
  1957, 1960, 1962 (request, not digitised); Abrams 1935 SE Michigan survey
  (Archives of Michigan, on site); USDA FSA 1957/1964/1972 and NRCS 1976 film
  (purchase); SEMCOG 1966-1995 (request); MiSAIL secured 2010/2012/2014/2015/2020
  (partner request); CORONA/HEXAGON 1964-84 (EE, some paid); Canadian NAPL rolls
  1931-1995 over the riverfront ($25/scan, account).

Disk: the new layers are ~60 GB of mosaics (115 GB total, gitignored). The NOAA
tile indexes used are in `data/sources/`.
