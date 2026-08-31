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

Residual against modern imagery on a 16 × 3 grid, over cells where the correlation
locked, after the full pipeline:

| block | n | median | p90 | max | cells >25 m | unlocked |
|---|---|---|---|---|---|---|
| 1961 | 34 | **2.7 m** | 11.9 m | 53.2 m | 1 | 8 |
| 1961 (ratio ≥ 1.25) | 28 | **2.6 m** | 5.1 m | **12.4 m** | **0** | — |

Spatially-blocked held-out station error — the honest out-of-sample figure, since
the grid above is interpolation between control — is **5.7 m** for 1961.

Against the un-warped composite: median 66.7 m, p90 109.0 m, max 161.7 m.

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
