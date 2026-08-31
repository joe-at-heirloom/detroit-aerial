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
python3 -m venv .venv && ./.venv/bin/pip install rasterio pillow numpy
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
| `mosaics/` | the actual rasters (~1.2 GB, gitignored) |
| `data/` | solutions, calibrations, control stations, manifest |
| `pipeline/` | the photogrammetry modules |

### Layers

- **Downtown** — 1949 / 1956 / 1961 / Today, on one 0.60 m/px grid
- **West Detroit** — the 62-frame 1961 block and 51-frame 1967 block, plus modern

---

## How the georeferencing works

Four stages. Each exists because the previous one was measurably wrong.

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
every 1609 m. This matters: the residential grid repeats every ~100 m, so every
correlation method aliases against it — a block shifted by one lattice vector
scores as perfectly aligned. Matching arterials only cannot alias inside a ±400 m
search.

**4. Local warp.** One global transform is right on average and wrong everywhere
in particular. So ~110 control stations per block are solved on a 2.6 km grid by
correlating the historical *ridge response* against the modern *ridge response*
(imagery, not road vectors -- 1.8% vector coverage was the accuracy bottleneck).

A global polynomial was fitted to those stations first and it failed at the block
ends, where it extrapolates past the last station and diverges -- 217 m on 1967's
southern chunk. It is now a **linear trend plus a Gaussian-kernel local
regression** (1000 m length scale, ridge-damped, with neighbourhood-based outlier
trimming). Outside station coverage the RBF term decays to zero and the linear
trend takes over, so extrapolation degrades gracefully instead of exploding.

Length scale and trend degree were chosen by **3-seed x 6-fold cross-validation on
held-out stations**, scored separately on the end bands:

| model | held-out overall | end bands |
|---|---|---|
| polynomial deg 3 (1967) | 44.0 m | 37.2 m |
| polynomial deg 4 (1961) | 31.8 m | 39.5 m |
| **linear + RBF 1000 m (1967)** | **35.9 m** | **30.2 m** |
| **linear + RBF 1000 m (1961)** | **25.8 m** | **38.0 m** |

## Accuracy, measured per chunk

Residual against modern imagery, in eight chunks down each block:

| block | median | max |
|---|---|---|
| 1961 | **8.5 m** | 18.2 m |
| 1967 | **7.3 m** | 14.1 m |

Against the earlier global-polynomial build (median 71.4 / 30.3 m, max 118.8 / 216.8 m).

Downtown epochs register to **3.0-5.0 m** against road centrelines.

### Known limits

- **End-of-block accuracy is now the RBF's**, not a polynomial's: outside station
  coverage the local term decays and the linear trend carries, so error grows
  slowly rather than diverging.
- **Not orthorectified.** Tall buildings lean differently in every negative, so
  rooftops jump between epochs even where the ground is aligned. Fixing that needs
  a DEM or dense stereo from the 60%-overlapping frames.
- **Coverage 75% (1961) / 62% (1967).** Gaps between flight lines, not holes.
- **Road-based validation bottoms out around 5 m** — modern satellite imagery
  scores the same against the same centrelines, so that is the metric's noise
  floor, not the imagery's error.
- `mosaics/modern_1961.png` is greyscale; recolour if you want the West Detroit
  "Today" side in colour.

---

## Note on iCloud

`~/Documents/Apps` is iCloud-synced and `mosaics/` is ~1.2 GB. It is gitignored,
but iCloud will still upload it. To keep it local only, either move `mosaics/`
outside Documents and symlink it, or append `.nosync` to the folder name and
update `data/manifest.json`.
