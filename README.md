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

The negatives are solved as a photogrammetric block with **COLMAP**, then placed on
the map with one low-order warp, and every step is judged by an instrument that
touches no external reference: how far two overlapping negatives disagree about
the same ground, measured on thousands of small phase-correlation windows across
every overlap, along each flight line and across neighbouring lines separately.

1. **Prepare the scans.** Crop the film edge, then pad every negative to one common
   canvas -- the scans come in several slightly different sizes, and COLMAP with
   one shared camera silently drops any image whose size differs.
2. **Solve the block.** SIFT features, matched only between negatives that plausibly
   overlap by the catalogue positions (a street grid that repeats every ~100 m
   gives confident matches between negatives that do not overlap at all), then
   COLMAP's incremental bundle adjustment with the camera fixed at values solved
   on a short pilot strip -- on flat terrain a free focal length is degenerate.
3. **Orthorectify onto the block's own surface.** A block bundled on flat ground
   bows into a shallow dome (144 m at a corner on 1961). Projecting each negative
   through its camera onto a quadratic surface fitted to the block's 3D points
   reproduces the geometry the cameras agree on; projecting onto a plane left
   12-16 m at every seam.
4. **Place the block.** Measure its displacement from modern imagery on a grid,
   with an alias-proof prior from the mile-grid arterials; fit the simplest model
   that explains it by held-out error (a quadratic, 3.2 m on 1961); apply it as one
   continuous warp of the composite, with the seams re-measured afterwards.

Scripts: `colmap_block.py` (2), `colmap_render.py` (3), `colmap_place.py` (3-4).

### Why not the earlier pipeline

The pipeline this replaced solved each negative's rotation from one window at the
overlap centre -- where rotation is invisible -- then moved each negative
independently to match modern imagery and fitted a flexible warp on top. That
produced excellent absolute numbers on a mosaic whose flight lines disagreed with
each other by 100-200 m: a road crossed a sidelap and appeared twice. A hand-rolled
tie-point bundle (`close3.py`) fixed that to ~4 m before COLMAP's real camera model
took it to 2 m. `HANDOFF.md` records the whole sequence, including the eight ways
a measurement lied before being trusted.

## Accuracy

Every number below comes from a metric proved first: it must recover a shift
planted on modern-vs-modern imagery (0.0 m), and a shift planted on the real mosaic
must move every cell by exactly that much (0.0 m).

**1961** (62 negatives, four flight lines) -- complete and in the viewer:

| | along-track seams | cross-line seams | absolute, 3.6 km cells | absolute, 1.8 km cells |
|---|---|---|---|---|
| COLMAP + surface | **2.0 m** (p90 4.0, max 16) | **2.0 m** (p90 4.5, max 20) | **2.6 m** (p90 5.8) | **3.0 m** (p90 8.0) |

Seams are identical before and after the placement warp. Independently, against
USGS NAIP on 1.5 km windows, three registration implementations put it at
**2.2 / 3.1 / 3.5 m** (scikit-image / OpenCV / ours).

**1967** (51 negatives, four flight lines) -- complete and in the viewer:

| | along-track seams | cross-line seams | absolute, 3.6 km cells | absolute, 1.8 km cells |
|---|---|---|---|---|
| COLMAP + surface | **2.0 m** (p90 4.5) | **4.0 m** (p90 7.2) | **2.6 m** (p90 7.7) | **4.1 m** (p90 11.6) |

**1949** (50 negatives, three flight lines) -- complete and in the viewer. Its
catalogue positions were poor, so the block's heading was refined against 1961's
verified build before placement; the absolute figures below are against modern
imagery, which it was not placed against:

| | along-track seams | cross-line seams | absolute, 3.6 km cells | absolute, 1.8 km cells |
|---|---|---|---|---|
| COLMAP + surface | **2.0 m** (p90 4.0) | **2.0 m** (p90 4.0) | **4.2 m** (p90 9.9) | **4.7 m** (p90 11.9) |

Independently, against USGS NAIP on 1.5 km windows, the three implementations
agree at **8.0 / 8.2 / 8.4 m** (scikit-image / OpenCV / ours) with no systematic
shift (bias under 1 m in both axes). That is worse than our own figure against
Esri imagery and worse than 1961's 2-3.5 m. Seventy-five years separate the
negatives from either reference, the farmland of 1949 is subdivisions now, and
the placement's own leave-one-out estimate was 6.8 m -- the remaining error is
the difficulty of measuring 1949 against anything modern, not the block's shape.

**1956** (70 negatives, four flight lines, 33 km long) -- complete and in the
viewer. Its catalogue placement was 815 m off; the heading refinement against
1961 now iterates with a wider search. The block is long enough that the
placement needed a cubic model (chosen by leave-one-out, 7.8 m):

| | along-track seams | cross-line seams | absolute, 3.6 km cells | absolute, 1.8 km cells |
|---|---|---|---|---|
| COLMAP + surface | **2.0 m** (p90 8.0) | **4.0 m** (p90 8.9) | **5.6 m** (p90 21.6) | **6.2 m** (p90 38.3) |

Where 1961 covers the ground the placement residual is 5.6 m median; the
western column and the northern 3 km lie outside every other epoch and were
measured against modern imagery only, where 1956 is as hard to match as 1949,
so the p90 there is the measurement as much as the map.

Downtown is a separate three-frame scene: 1961 and 1949 corrected to ~5 m, 1956
still ~32 m; much of the frame is the river and the core is high-rise, so no 2D
correction can do better there.

### Known limits

- The surface is a quadratic, not terrain. Detroit is flat enough that this is the
  ground to a couple of metres; it would not be elsewhere.
- Tall buildings lean differently in every negative, so rooftops jump between
  epochs even where the ground is aligned. Fixing that needs a dense surface from
  the 60%-overlapping frames.
- Coverage follows the flight lines; gaps between them are gaps, not errors.

## Note on iCloud

`~/Documents/Apps` is iCloud-synced and `mosaics/` is several GB. It is gitignored,
but iCloud will still upload it, and two mosaics disappeared from disk during a
working session. Everything there is regenerable from `data/` plus the cached
frames — do not treat it as durable storage. To keep it local only, move `mosaics/`
outside Documents and symlink it, or append `.nosync` to the folder name.
