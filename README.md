# Detroit Air Survey

Georeferenced historical aerial photography of Detroit, served locally at native
resolution. Wipe between **1949, 1956, 1961, 1967**, the City of Detroit's
**1998, 2005 and 2010** orthophotos, USDA NAIP for **2012 to 2022**, and today.

Source imagery: Wayne State University, **DTE Aerial Photo Collection** —
© DTE Energy, administered by the Walter P. Reuther Library. 1998, 2005 and 2010
orthophotos: City of Detroit (egis.detroitmi.gov). 2012-2020: USDA NAIP from the
State of Michigan's image services (imagery.michigan.gov); 2022: USDA NAIP from
the USGS National Map. Modern imagery: Esri.
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
  blocks; then the City of Detroit's 1998 (1 m), 2005 and 2010 (2 ft)
  orthophotos and USDA NAIP for 2012, 2014 (1 m), 2016, 2018, 2020 and 2022
  (0.6 m), fetched onto the same grid by `scripts/fetchlayer.py`; plus modern

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
verified build before placement. The placement itself is fitted directly against
modern imagery with a cubic model:

| | along-track seams | cross-line seams | absolute, 3.6 km cells | absolute, 1.8 km cells |
|---|---|---|---|---|
| COLMAP + surface | **2.0 m** (p90 4.0) | **2.0 m** (p90 4.0) | **4.2 m** (p90 7.1) | **4.1 m** (p90 8.4) |

It was first placed by chaining through 1961, on the reasoning that a 12-year gap
matches far more easily than a 75-year one. That was wrong here: chaining
inherits 1961's own 2.6 m, and once the heading refinement had removed the gross
error there was enough signal to fit against modern imagery directly. Fitting
direct, and with a cubic rather than a quadratic, halved the error against an
independent reference -- USGS NAIP, which no build is fitted to -- from
**8.0-8.4 m** to **4.4-5.3 m**, with the worst of twelve windows going from
14.9 m to 7.5 m.

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
Independently, against USGS NAIP on 1.5 km windows, the three implementations
put it at **4.5 / 4.6 / 5.9 m** (scikit-image / OpenCV / ours), p90 21-31 m,
bias under 3 m: the same picture, a good centre and uncertain edges.

Downtown is a separate three-frame scene: 1961 and 1949 corrected to ~5 m, 1956
still ~32 m; much of the frame is the river and the core is high-rise, so no 2D
correction can do better there.

Those were *absolute* figures -- each epoch measured against a modern reference on
its own. The wipe shows a **pair**, and two epochs can each be 15 m out in
opposite directions and be 30 m apart on screen. `scripts/dtcross.py` measures the
pair directly, on a 6x6 grid at 1.25 m/px:

| | before | after | p90 | max |
|---|---|---|---|---|
| 1949 vs 1956 | 49.8 m | **3.9 m** | 10.8 | 31.1 |
| 1956 vs 1961 | 37.9 m | **3.2 m** | 14.1 | 30.7 |
| 1956 vs today | 41.5 m | **15.2 m** | 64.6 | 86.2 |
| 1949 vs 1961 | 3.9 m | 3.9 m | 31.8 | 56.0 |
| 1949 vs today | 20.2 m | 20.2 m | 32.9 | 42.7 |
| 1961 vs today | 23.1 m | 23.1 m | 52.1 | 56.3 |

### Downtown 1956 was 168 m out, not 32

Every measurement of it had been aliasing, and the alias is the one this project
already documents for the west blocks: Detroit's street lattice repeats every
~97.5 m, and the fine search is capped at +/-40 m, walking about 120 m over three
re-centrings. An epoch further out than that cannot be reached -- and it does not
fail loudly. It locks onto the wrong member of the lattice and reports a
confident small number. 1956 read ~32 m for years for exactly that reason.

Solved on the whole raster with a 250 m search it is 170 m out, at a coarse peak
ratio of 1.90, and all sixteen cells of a 4x4 grid independently agree to within
+/-20 m. `fixdowntown.py` now solves that one rigid shift before the residual
field (`global_prior`), and applies it only when the peak is convincing and the
scene is beyond the fine search's reach -- on 1949 and 1961 it finds 14-15 m at
ratios of 1.16-1.26 and correctly declines to use it.

With the shift in front of it, the residual field's control goes from 39 stations
(23 at the search limit, 15.1 m held-out) to 62 stations (none at the limit,
**2.0 m** held-out), and the grid goes 168.4 -> 2.2 m median, no cell over 10 m.
Rebuild it with:

    ./.venv/bin/python scripts/fixdowntown.py 1956 --ref 1961 --src

`--src` re-runs from the uncorrected raster; without it a second run reads the
manifest, which by then points at the *corrected* mosaic, and stacks a second
warp on the first.

### What is left downtown

The historical epochs now agree with each other at 3-4 m. What remains is the
"versus today" column, and it is not a placement error. Against USGS NAIP -- an
independent reference no downtown build is fitted to -- today's layer measures
4.0 m and 1961 measures 6.2 m, so neither is 23 m out of place. 1949 and 1956
measure 19-29 m against NAIP with only ~40% of cells locking at all, which is the
difficulty of matching seventy-year-old imagery to modern downtown, where the
whole core was rebuilt. The rest is relief displacement: a 560 m cell downtown is
mostly high-rise roofs, and they lean differently in every negative. Separating
that from the ground needs a DEM, not another 2D fit.

## Hosting it (GitHub Pages, Vercel, anything static)

The local viewer cuts tiles out of 54 GB of GeoTIFFs on demand. A static host
serves files, so the tiles have to exist first:

    ./.venv/bin/python scripts/build_static.py

That renders every film layer (four blocks and downtown) to a WebP pyramid at
its native z18, the road centrelines to z17, and lays out `dist/` -- the viewer,
the manifest, the places -- with everything relative, so it works under a repo
subpath. **Today** is Esri World Imagery fetched straight from Esri, so it costs
nothing to host. Hand alignments in `data/adjust.json` are baked into the tiles.
It is resumable: a killed run picks up where it stopped, and a rerun only renders
what is missing.

Measured, not estimated: **686 MB and 114,000 tiles** for the film at native
resolution (WebP q80 is 4-12 KB a tile where PNG is 30-140). That is inside what
a GitHub Pages repo will hold, and Vercel's 100 MB-per-file cap is nowhere near.
It takes about **106 minutes** on nine workers -- nine processes reading the same
GeoTIFFs at z18 contend on I/O, so more cores do not help much.

The ten modern layers (1998-2022) are not built by default; `--layers all
--zmax modern:17` adds ~2.7 GB, which wants object storage rather than a repo.
`--tiles-base https://your-bucket` writes the tiles locally as usual but points
`dist/index.html` at the bucket, so `dist/` on Pages stays a few hundred KB and
the tiles live somewhere built for it (Cloudflare R2's free tier is 10 GB with
no egress charge).

Then either:

- **GitHub Pages** -- commit `dist/`, set Settings -> Pages -> Source to "GitHub
  Actions", and `.github/workflows/pages.yml` publishes it on every push that
  touches `dist/`. The first push of 114k files is slow; after that only changed
  tiles move.
- **Vercel** -- `vercel.json` already points at `dist/` with no build step and
  a year-long cache header on tiles. Import the repo and deploy.

What the static build cannot do: no ALIGN handle (nothing to write to; it hides
itself) and no "+ PLACE" -- both are local-server tools. Everything else --
the wipe, the probe, the places, the filter -- is identical, and it was verified
by serving `dist/` with `python3 -m http.server` and driving it, not by reading
the code.

## Aligning by hand

Some things the correlation cannot be talked into, and some you just want to nudge.
Tick **ALIGN** in the viewer: the layer you pick is drawn semi-transparent over
whatever is on the other rail, and you move it.

| | |
|---|---|
| drag | move |
| shift-drag | rotate about the layer's centre |
| alt-drag | scale about the same point |
| arrows | nudge one pixel, ten with shift |
| space-drag | pan the map instead |
| hold **F** | blink the layer off to check the fit |

**SAVE** writes four numbers to `data/adjust.json` and `scripts/serve.py` applies
them from then on, so the alignment is what everything sees. Nothing is baked into
the rasters: it stays reversible, and it stays readable as

    "dt1956": { "dE": 28.8, "dN": -167.8, "deg": 0.0, "scale": 1.0 }

Dragging cannot wait for a round trip, so the live picture is a CSS transform on a
pane of its own and the server is asked only on save. Those two have to agree
exactly or the map would jump the moment you saved, so the preview draws the
*difference* between what your hands are doing and what the server already holds --
save, the difference becomes the identity, and the transform falls away by itself.

The measurement follows it. `scripts/dtcross.py` reads the same file through the
same geometry (`pipeline/handadjust.py`, which the viewer, the tile server and the
metric all come through, because two implementations of one transform is how this
project has gone wrong before). So a hand alignment can be checked rather than
trusted: plant 50 m by hand and the metric reports 49.3 m.

    ./.venv/bin/python scripts/dtcross.py       # prints the hand offsets it found

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
