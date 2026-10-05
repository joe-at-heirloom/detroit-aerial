<img src="viewer/logo.svg" width="72" alt="">

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
| `viewer/index.html` | Leaflet viewer: year slider; per-scene tabs with epoch rails and a swipe divider |
| `scripts/fetchlayer.py` | pull a published ortho (ImageServer, tile cache, NAIP) onto the project grid |
| `scripts/footprints.py` | each layer's real coverage outline, and what the year view shows per year |
| `scripts/maskfill.py` | turn a layer's white out-of-coverage fill into nodata |
| `mosaics/` | the actual rasters (gitignored) |
| `data/` | solutions, calibrations, control stations, manifest |
| `pipeline/` | the photogrammetry modules |
| `scripts/validate.py` | calibrate the metric, then measure a block |
| `scripts/rebuild.py` | composite -> per-frame -> residual field -> final mosaic |
| `scripts/montage.py` | alignment sheets: historical in red, modern in green |

### Layers

- **Downtown** — 1949 / 1956 / 1961 / Today, on one 0.60 m/px grid
- **West Detroit** — the 1949 (50 frame), 1956 (70), 1961 (62) and 1967 (51)
  blocks; then published orthophotos fetched onto the same grid by
  `scripts/fetchlayer.py`: the City of Detroit's 1999 (published as 1998), 2005
  and 2010, USGS 2002, NOAA 2008, State/USDA NAIP 2005-2022, NHAP 1983 and 1987,
  Wayne County 2015, the State's 2024, Michigan Tech's 1951; Oakland County's 15
  years north of 8 Mile (1940-2025), the City of Taylor's (1940-2021), Essex
  County's 1931 Detroit River, NOAA 2025/2026 along the river, and Esri Wayback
  2019/2023 as remote tiles. HANDOFF.md "Every year we could find" lists them
  with their measured placement. Plus modern imagery as **Today** -- locally a
  1.8 m Esri fetch on the same grid (`mosaics/modern_west_hi.tif`), on the
  static site live Esri World Imagery

### Through the years

The first tab is one map with a year slider instead of a wipe. For the chosen
year, the newest picture taken in or before it is on top, and wherever an older
one reaches further -- past the edge of a film block, outside the City's
orthophoto at the city limits, south of Oakland County's at 8 Mile -- the older
one shows there instead of nothing. Each picture's visible ground is outlined
thinly and captioned with its year: amber for the year picked, white for older
ground filling in. Clicking the map opens the probe on every layer whose
outline holds the point.

The outlines are each raster's real coverage, not its bbox:
`scripts/footprints.py` reads every layer at ~15 m cells (any band non-zero =
ground, exactly what `serve.py` draws as opaque), cleans and vectorises it, and
for every year works out which part of which layer no newer one covers. It
writes `data/timeline.json`; run it after `manifest.py`, then restart
`serve.py`:

```bash
./.venv/bin/python scripts/manifest.py && ./.venv/bin/python scripts/footprints.py
```

A source that fills outside its coverage with white instead of black (the
City's 2005 and 2010 orthos did) hides every older year there; `scripts/maskfill.py
<name>` turns large white regions of `mosaics/layer_<name>.tif` into nodata.
Two layers of the same year are ordered by `prio` in their `layer_<name>_geo.json`
(`fetchlayer.py --prio`), higher on top.

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

Every number below comes from an instrument proved first: it must recover a shift
planted on modern-vs-modern imagery, and each cell is re-measured with a known
field (30 m + 1.5 m/km scale + 1 mrad rotation) planted on one side -- a cell that
does not follow it by 5 m is dropped. Across every pair measured it tracks to
0.8 m (p90 1.1), so differences under a metre are not evidence. Figures are
medians over 1.3 km cells (p90 in brackets), measured 2026-09-27 by
`scripts/crossmatrix.py`; the full pair table is `runs/crossmatrix_west.md`.

**West Detroit, the film against modern imagery:**

| | vs Esri (never fitted to) | vs NAIP 2016 | seams: along / across the flight lines |
|---|---|---|---|
| **1961** (62 negatives) | **2.8 m** (6.6) | 3.4 m (7.0) | 1.7 / 2.1 m (p90 4.0 / 5.0) |
| **1967** (51) | **3.3 m** (7.1) | 3.0 m (6.1) | within a 2 m pixel |
| **1949** (50) | **4.4 m** (10.0) | 5.5 m (11.7) | within a 2 m pixel |
| **1956** (70) | **7.2 m** (41) | 7.2 m (47) | within a 2 m pixel along, 4 m across |

Seams are how far two overlapping negatives disagree about the same ground. The
tie windows that measure them report whole pixels, 2 m at the resolution they run
at, so "within a 2 m pixel" is the honest reading of the "2.0 m" this project
quoted for months; 1961 re-measured with a sub-pixel peak
(`distortion_pilot.py --subpixel`) gives the figures in its row.

**What a wipe shows** is a pair, not two absolute numbers: 1961/1967 **3.8 m**
(8.0), 1949/1956 4.6 m (13.7), 1956/1961 4.8 m (13.7), 1949/1961 6.0 m (11.1),
1949/1967 6.6 m (14.2).

**The modern layers** are served as published. Against NAIP 2016: 2022 2.1 m,
2010 2.0, 2005 3.1, 2014 3.7, 1998 3.8, 2012 4.0, 2020 5.1, Esri 1.9. NAIP 2018 is
3.0 m except for a strip east of Greenfield that sits 10-14 m east (a shifted
quarter-quad column in the product; p90 17 m).

Notes by block:

- **1961** is the reference build: the camera for every block was solved on a
  six-frame 1961 pilot. A pilot on the rebuilt 1961 solve found about a metre of
  radial distortion at the frame edges that the fixed camera does not model; fitting
  it changed nothing measurable (HANDOFF, "Alignment plan, steps 1-4").
- **1967** is placed by its original quadratic plus a second quadratic from the
  joint placement (`jointfit.py`, `jointapply.py`): every film year fitted at once
  to NAIP 2016 and to each other, the model chosen on 4 km blocks withheld. It is
  1-1.8 m closer to every modern year it was never fitted to than the build before.
- **1949**'s catalogue positions were poor, so its heading was refined against
  1961 before placement, and the placement is fitted straight to modern imagery
  with a cubic. Chaining it through 1961 first -- a 12-year gap should match more
  easily than a 75-year one -- inherited 1961's own error; fitting direct halved the
  error against USGS NAIP (8.0-8.4 m to 4.4-5.3 m on three independent
  registration implementations). The farmland of 1949 is subdivisions now, and
  that, not the block's shape, is what limits it.
- **1956** runs 33 km and needed a cubic placement. Where 1961 covers the ground
  it is good; its western column and northern edge were farmland that is suburbs
  now, where image matching locks onto the wrong streets -- 100 m "errors" that
  repeat identically against every modern year, while the one road present in both
  years says about 20 m. Its p90 is that edge, and it is unverified rather than
  measured. No smooth correction survives withheld validation there.

**Downtown** is a separate three-frame scene. On the ground the historical years
agree to **2.4-3.7 m** (1956/1961 2.4, 1949/1961 2.9, 1949/1956 3.7). Against
today image matching does not measure the ground at all: the core is high-rise,
its roofs lean 25-35 m and differently in every photograph, and correlation locks
onto the roofs. City street-centreline vectors sit on the streets in 1949, 1961
and today to a few metres. A hand alignment made by lining up roofs put 1949
35-89 m off the other years and was removed (2026-09-27). A true downtown ortho
needs a surface model; no 2D fit will make the roofs line up.

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
ratios of 1.16-1.26 and correctly declines to use it. Rebuild it with:

    ./.venv/bin/python scripts/fixdowntown.py 1956 --ref 1961 --src

`--src` re-runs from the uncorrected raster; without it a second run reads the
manifest, which by then points at the *corrected* mosaic, and stacks a second
warp on the first.

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

The year view on the static site covers only the layers that build carries:
`build_static.py` cuts its own `dist/data/timeline.json` from exactly those (a
year's regions assume every layer in its stack is there). The default film build
gets seven years -- 1949, 1956, 1961, 1967, and Esri Wayback 2019/2023 and Today
as remote tiles. All 41 need `--layers all` and somewhere to put ~8 GB of tiles
(z17), i.e. `--tiles-base`.

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
