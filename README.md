# Detroit Air Survey

Georeferenced historical aerial photography of Detroit, served locally at native
resolution. Wipe between **1949, 1956, 1961, 1967** and today.

Source imagery: Wayne State University, **DTE Aerial Photo Collection** —
© DTE Energy, administered by the Walter P. Reuther Library. Modern imagery: Esri.
Street centrelines: City of Detroit open data (+ OpenStreetMap outside city limits).

---

## Run it

```bash
cd "Detroit Air Survey"
python3 -m venv .venv && ./.venv/bin/pip install -r requirements.txt
./.venv/bin/python scripts/serve.py --port 8770
```

Then open <http://localhost:8770>.

There is **no tile pyramid on disk**. `scripts/serve.py` reads windows straight out
of the GeoTIFFs on demand, so you can zoom to the native **0.63 m/px** of the
negatives without pre-rendering anything.

---

## What's here

| Path | |
|---|---|
| `scripts/serve.py` | local XYZ tile server + static host |
| `viewer/index.html` | Leaflet viewer: tabs, epoch rails, swipe divider |
| `mosaics/` | the actual rasters (gitignored) |
| `data/` | solutions, calibrations, control stations, manifest |
| `pipeline/` | the photogrammetry modules |
| `scripts/validate.py` | calibrate the metric, then measure a block |
| `scripts/rebuild.py` | composite -> per-frame -> residual field -> final mosaic |
| `scripts/montage.py` | alignment sheets: historical in red, modern in green |

### Layers

- **Downtown** — 1949 / 1956 / 1961 / Today, on one 0.60 m/px grid
- **West Detroit** — the 1949 (50 frame), 1956 (70), 1961 (62) and 1967 (51)
  blocks, plus modern

---

## How the georeferencing works

Seven stages. Each exists because the previous one was measurably wrong.

**1. Relative orientation.** Overlapping frames are matched by windowed phase
correlation swept over ±3° of rotation. Aerial film carries a per-frame *crab
angle*; measured spread across the 1961 block is **3.0°**, and 1° over a 3.4 km
frame throws the corners 30 m. An earlier translation-only bundle is why the
stitching looked crooked.

**2. Bundle adjustment.** Rotations solve first (θⱼ − θᵢ = Δθᵢⱼ, robust least
squares), then positions with rotations fixed.

| | translation-only | with rotation |
|---|---|---|
| 1961 tie residual | 2.9 m | **0.5 m** |
| 1967 tie residual | 2.2 m | **0.8 m** |

**3. Absolute orientation.** Against Detroit's **mile-grid arterials**, which repeat
every 1609 m. This matters: the residential grid repeats every **97.5 m** here
(measured, ridge autocorrelation r = 0.46–0.67), so every correlation method
aliases against it — a block shifted by one lattice vector scores as perfectly
aligned. Matching arterials only cannot alias inside a ±400 m search.

**4. Composite** the block from those placements, nearest-centre (Voronoi) so seams
fall midway between exposures.

**5. Per-frame adjustment.** This is what fixes the worst of it. A continuous warp
cannot represent the *step* in error where two negatives meet, so it splits the
difference and leaves half the disagreement on either side. Each negative is
instead rendered alone into map space — its middle 78%, avoiding the vignetted
edge and the worst relief displacement — matched against modern imagery, and fed
back as its own (dE, dN). The solve iterates: a frame that fails to lock on the
first pass did so with a prior tens of metres off. Frames that still never lock
borrow their neighbourhood's consensus rather than staying put, which would tear
the mosaic exactly where it used to be continuous.

On 1961 this alone takes the block from **66.7 m to 6.0 m** median.

**6. Residual local field.** ~260 control stations solved coarse-to-fine and
iterated to convergence, then a linear trend plus a Gaussian-kernel local
regression. Length scale by **spatially-blocked** cross-validation — random folds
are optimistic here because the control windows overlap, so a held-out station
almost always has a near-duplicate left in training (random folds claimed 5.1 m
where spatial folds said 9.1 m).

**7. Apply once**, on a displacement lattice that scales with the chosen length
scale.

---

## Measuring it honestly

Everything above is only as good as the ruler, and the ruler was bent twice.

**The metric aliased.** Correlating a ridge response over a ±200 m window includes
the 97.5 m residential lattice, which gives six alias positions per axis. It
reported cell (5,1) of 1961 as `dE = +97.5 m` — exactly one block. So no single
search is allowed to be ambiguous: a regional coarse field on 6 km windows (mile
grid, cannot alias), then per-cell refinements of at most **±40 m**, repeated.
Self-consistency went from *405 m out* to **0.0 m**.

**The modern reference was bent.** It was Web Mercator with a lat/lon bounding box
stamped on it. Against a correctly reprojected reference it bows by −7.4, −13.8,
−17.3, −17.2, −13.8, −6.4 m up the raster — a symmetric parabola, ~17 m at
mid-span, where theory for that error predicts 15.5 m. `scripts/fetchmodern.py`
refetches at **1.77 m/px** and resamples each row from its own true Mercator
latitude.

`scripts/validate.py` now proves the metric before reporting anything: it must
recover a planted shift on modern-vs-modern, and a planted shift on the real
mosaic must move every cell's measurement by exactly that much.

**Confidence is reported, not hidden.** A cell whose correlation has no unique peak
does not have a large error — it has an unknown one. Quoting it as error is how
this project once got a "231.9 m worst cell" out of a mosaic that was roughly
right. Cell (5,1) is the worked example: it reports 39.9 m at ratio 1.07, and
applying that correction visibly splits an arterial that was already aligned.

---

## Accuracy

Every figure below comes from a metric that is proved before it is believed. Per
block, it must recover a shift planted on modern-vs-modern imagery (**0.0 m**, all
four), and a shift planted on the real mosaic must move every cell's measurement by
exactly that much (**0.0 m**, all four).

Residual against modern imagery, over cells where the correlation locked:

| block | frames | 3.6 x 2.0 km cells | 1.8 x 1.0 km cells |
|---|---|---|---|
| 1961 | 62 | **1.3 m** — p90 4.8, max 6.8, none over 10 m | **3.5 m** |
| 1967 | 51 | 3.0 m — p90 9.8, max 11.6, none over 25 m | 8.7 m |
| 1949 | 50 | 3.7 m | 18.5 m |
| 1956 | 70 | 5.0 m — p90 10.4 | 12.9 m |

Against the un-warped composites: 58-116 m median.

### Verified independently

None of that is worth much on its own, because it is measured against a reference
this project built with code this project wrote. So it is checked against **USGS
NAIP** — public-domain orthoimagery from a different organisation and a different
processing chain, delivered in EPSG:4326 so *their* server does the reprojection
and ours leaves the loop — using **three registration implementations**, each
admitted only after recovering a planted shift on the real data.

The reference itself agrees with NAIP over 24 windows: **median 2.16 m, max 7.56 m**.

The blocks, on 1.5 km windows, on the road response:

| block | ours | scikit-image | OpenCV ECC |
|---|---|---|---|
| 1961 | 3.84 m | 3.16 m | 3.79 m |
| 1949 | 3.47 m | 3.50 m | 5.47 m |
| 1956 | 7.74 m | — | 12.71 m |
| 1967 | 16.22 m | 8.64 m | 11.57 m |

Note 1949: independent measurement puts it at 3.5-5.5 m where our own fine grid
says 18.5 m. The grid divides the block into fixed cells and many of 1949's are
only partly covered (coverage 0.63), which measures badly; the independent windows
require 97% coverage, so they test the imagery where it actually exists. Read the
two together — the grid is the pessimistic bound.

And the whole serving path is verified, not just the files: tiles pulled from the
running viewer at zoom 17 and overlaid on modern show single yellow streets for all
four blocks, which exercises the Mercator conversion `serve.py` performs.

### Downtown

Downtown is a separate, weaker story, and it had never been checked against
anything independent. The claimed 3-5 m came from road-vector validation, which
this project had already established has a 5-6 m noise floor. Measured against
NAIP it was **20-42 m** out. Nobody caught it because downtown sits east of the
modern reference raster entirely, so the modern side was all zeros and every cell
silently reported "no peak".

| | before | after |
|---|---|---|
| downtown 1961 | 20.7 m | **4.5 m** |
| downtown 1949 | 37.4 m | **5.3 m** |
| downtown 1956 | 42.3 m | 32.5 m |

Downtown measurement is unreliable and should be read with care: a large part of
the frame is the Detroit River, which has nothing to correlate, and the core is
high-rise, so relief displacement leans every building differently in every
negative.

### Known limits

- **Not orthorectified.** Tall buildings lean differently in every negative, so
  rooftops jump between epochs even where the ground is aligned. Fixing that needs
  a DEM or dense stereo from the 60%-overlapping frames. Relief alone accounts for
  several metres at frame edges.
- **Land-use change is not misalignment.** Where a freeway or subdivision did not
  exist in the historical epoch there is nothing to align to, and the correlation
  there is weak by construction. Those cells are reported as unlocked.
- **The reference is 1.77 m/px**, so measurements below ~2 m are not meaningful.
- **Coverage** is 61–76% of each block's extent — gaps between flight lines, not
  holes.

---

## Note on iCloud

`~/Documents/Apps` is iCloud-synced and `mosaics/` is several GB. It is gitignored,
but iCloud will still upload it, and two mosaics disappeared from disk during a
working session. Everything there is regenerable from `data/` plus the cached
frames — do not treat it as durable storage. To keep it local only, move `mosaics/`
outside Documents and symlink it, or append `.nosync` to the folder name.
