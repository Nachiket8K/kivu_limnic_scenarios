# Lake Kivu Scenario Player — Web Viewer Plan (Staged)

This plan turns the current offline simulation into a browser-based scenario “video player” with MapLibre, precomputed scenario bundles, PNG overlays per timestep, and OSM-derived vector tiles generated from OSMnx features.  
Primary goals:
- **No physics in the browser** (fast playback).
- **Scenario selection via dropdown** (limited runs initially).
- **AOI-bounded view** (no wasted rendering).
- Save **GeoTIFFs for raw fields** (for analysis/3D later), but use **PNG overlays** for UI playback.
- Use **FastAPI** for local/dev first (debug-friendly), then **GitHub Pages** for static hosting.

---

## Stage 0 — Freeze the POC contract

**Decisions (POC):**
- Reporting timestep: **10 minutes** (fixed, not user-configurable).
- Scenario selection: **dropdown** (precomputed runs only).
- Viewer renders: **MapLibre** basemap (custom OSM vector tiles) + **PNG overlays** per timestep.
- Output storage: keep both:
  - **metrics CSV** (time series)
  - **manifest JSON** (scenario metadata)
  - **PNG frames** (one per timestep per layer)
  - optional **GeoTIFFs / Zarr** for raw arrays (for 3D and analysis)

**Deliverable:** A documented “scenario bundle” spec (see Stage 1).

---

## Stage 1 — Define a scenario bundle spec (files + metadata)

Create one folder per scenario:

```
scenarios/
  scenario_<id>/
    manifest.json
    metrics_10min.csv
    frames/
      ppm_z/
        t0000.png
        t0001.png
        ...
      lethal_mask/
        t0000.png
        ...
      pdeath/
        t0000.png
        ...
      base_h/            # optional (visual thickness)
        t0000.png
    raw/
      h_base.tif         # optional
      M_base.tif         # optional
      M_top.tif          # optional
      ppm_z.tif          # optional (or one GeoTIFF per timestep if needed)
      pdeath.tif         # optional
    geo/
      aoi.geojson
      lake_outline.geojson
      waterbodies.geojson
      osmnx.pmtiles      # vector tiles packaged as PMTiles
```

### `manifest.json` (minimum fields)
- `scenario_id`
- `created_utc`
- `aoi`: bbox in WGS84 + `maxBounds` for MapLibre
- `crs`: e.g., EPSG:32735
- `grid`:
  - `width`, `height`
  - `dx_m`, `dy_m`
  - affine transform coefficients OR bbox + resolution
- `time`:
  - `dt_report_min` (10)
  - `n_frames`
  - `timestamps_min` (array or computed)
- `parameters` (scenario knobs):
  - `release_fraction`, `released_volume_m3_stp`
  - `eruption_radius_m`
  - `wind_speed_mps`, `wind_dir_from_deg`
  - other key parameters you want surfaced
- `layers` (viewer configuration):
  - list of layers with `{name, type, pathTemplate, opacityDefault, legend}`
  - example: `frames/ppm_z/t{frame:04d}.png`
- `assets`:
  - paths to `metrics_10min.csv`, preview thumbnail, etc.

**Deliverable:** A Python function `write_manifest(run_dir, …)` that creates this consistently after each simulation.

---

## Stage 2 — Extend the simulation exporter (offline) to produce frames + raw data

### 2.1 Decide which layers you export
For POC, export at least:
- `ppm_z` (continuous heatmap, main overlay)
- `lethal_mask` (binary ≥100k ppm)
- `danger_mask` (binary ≥50k ppm)
- `Pdeath` (0..1)
- Optional: `h_base` and/or `M_base`, `M_top` for debugging and future 3D

### 2.2 Export strategy
- **Frames:** PNGs (RGBA) with transparent background outside plume.
- **Raw:** GeoTIFFs for selected fields (for analysis and 3D later).
- Metrics already saved.

### 2.3 Rendering PNG overlays (important)
You need consistent:
- colormap + min/max scaling strategy (fixed or percentile-based)
- alpha handling (transparent where value below threshold)
- georeferencing metadata in `manifest.json` so the frontend can place the image correctly:
  - either **bounds** (WGS84) for the PNG overlay
  - or provide the projected bbox + conversion (frontend will use WGS84)

**Recommended approach for MapLibre:**  
Store overlay bounds in **WGS84** per scenario:
- Convert scenario grid bbox (UTM) → WGS84 polygon/bounds.
- In the frontend, use MapLibre `ImageSource` with 4 corner coordinates.

### 2.4 Add a thumbnail
Create a `thumb.png` for each scenario (e.g., max ppm at peak frame).

**Deliverables:**
- `export_frames.py` (or notebook cell) that takes arrays per timestep and writes:
  - `frames/<layer>/t####.png`
- `export_raw_geotiff.py` that writes selected raw outputs (optional at first)
- `write_manifest.py`

---

## Stage 3 — OSMnx → vector tiles pipeline (offline preprocessing)

Goal: generate your own OSM-based context layer **once per AOI**, reusable across scenarios.

### 3.1 Extract features with OSMnx
From AOI polygon:
- roads (highway)
- water (natural=water, water=lake/river)
- buildings (optional)
- admin boundaries (optional, or keep your GADM)
- place labels (optional)

### 3.2 Normalize + simplify
- Reproject to WGS84
- Simplify geometry (Douglas-Peucker) to reduce tile size
- Keep only useful properties (e.g., highway type, name)

### 3.3 Build vector tiles
Pick one:
- **PMTiles** (recommended for GitHub Pages): single-file tile archive
- or MBTiles then convert to PMTiles

Tools:
- `tippecanoe` (classic) → MBTiles
- `pmtiles` tools → PMTiles

**Deliverable:**
- `geo/osmnx.pmtiles` + style JSON (line colors by highway type, water fill, etc.)

---

## Stage 4 — FastAPI “scenario server” (debug-first)

This stage is local/dev to iterate quickly.

### 4.1 API endpoints
- `GET /scenarios`  
  Returns list of scenario manifests (or a subset + summary).
- `GET /scenarios/{id}/manifest`  
- `GET /scenarios/{id}/metrics`  
- `GET /scenarios/{id}/frames/{layer}/{frame}`  
  Returns the PNG.
- `GET /tiles/osmnx.pmtiles` (or serve static file)

### 4.2 Scenario parameter dropdown support
Build a small index file generated offline:
- `scenarios/index.json` with:
  - scenario id
  - parameter values used in that scenario
  - path to manifest
This allows the UI to build dropdown options like:
- release fraction: [0.01, 0.05, 0.1]
- wind dir: [0, 90, 180]
- wind speed: [0, 2, 4]
…and map a selection to a scenario id.

**Deliverables:**
- `app/main.py` FastAPI server with static file mount + endpoints
- `scenarios/index.json` generator

---

## Stage 5 — Frontend viewer (MapLibre + overlays + timeline)

### 5.1 UI layout
- Scenario selector panel:
  - dropdowns populated from `/scenarios` or `/scenarios/index.json`
- Time controls:
  - slider (frame index)
  - play/pause, speed (x1/x2/x4)
- Layer toggles:
  - ppm_z, lethal, Pdeath, base_h, etc.
- Opacity sliders per layer

### 5.2 Map setup
- MapLibre with:
  - max bounds = AOI bbox
  - min/max zoom
- Base layer:
  - your PMTiles vector tile source
- Overlay:
  - MapLibre `ImageSource` for current frame PNG
  - on slider change: update the URL (or swap source image)

### 5.3 Performance
- Preload a small window of frames around current timestep (e.g., t-2…t+2)
- Cache images in browser memory
- Keep only AOI-sized images

**Deliverables:**
- `frontend/` (Vite + React) or simple static HTML/JS
- A “player” component that swaps frames by timestep

---

## Stage 6 — GitHub Pages deployment (static hosting)

Because frames + manifests are static, GH Pages can host everything.

### 6.1 Structure for GH Pages
- Put `frontend/dist` plus `scenarios/` bundles in `docs/` (or use a separate gh-pages branch).
- Use relative paths so it works on Pages.

### 6.2 Serving PMTiles
PMTiles can be fetched as a static file; MapLibre loads it with a PMTiles protocol handler (small JS library).

### 6.3 Replace FastAPI endpoints
Switch frontend from `http://localhost:8000/scenarios/...` to relative:
- `./scenarios/index.json`
- `./scenarios/scenario_001/manifest.json`
- `./scenarios/scenario_001/frames/ppm_z/t0000.png`

**Deliverable:** A build script that copies scenarios into the GH Pages output folder.

---

## Stage 7 — Collect data for 3D (render later)

You want 3D eventually; collect the data now.

### 7.1 What to store for 3D
Per timestep (or reduced timesteps):
- `h` (base layer thickness)
- `M_base` and `M_top`
- terrain DEM (already in processed)
- optional derived: slope magnitude, flow speed

Store as:
- GeoTIFF stacks (simple)
- or Zarr with `(time,y,x)` for quick slicing

### 7.2 3D rendering candidates (future)
- `deck.gl`:
  - `TerrainLayer` + a `GridCellLayer` for thickness
- `three.js`:
  - heightfield mesh for DEM
  - translucent surface for h / ppm

**Deliverable now:** Export Zarr or per-timestep GeoTIFFs for h/M fields.

---

## Stage 8 — Scaling scenario library

When you add more scenarios:
- generate them offline (grid search over parameters)
- add them to `scenarios/index.json`
- viewer dropdown updates automatically

Later enhancements:
- “closest match” selection when a dropdown combination doesn’t exist
- side-by-side scenario comparison
- click-on-map → time series at that pixel (requires raw arrays, not only PNG)

---

## Execution checklist (what to do next)

1. **Exporter**
   - Implement `write_manifest()`
   - Implement `export_png_frames()` for ppm/lethal/Pdeath
   - Implement optional `export_raw_geotiff()` for h/M
   - Generate `scenarios/index.json`

2. **OSM tiles**
   - OSMnx extract → simplify → tippecanoe → pmtiles
   - Save style JSON

3. **FastAPI**
   - Serve `scenarios/` as static + simple endpoints
   - Test locally

4. **Frontend**
   - MapLibre base map with PMTiles
   - Overlay player: slider + play + toggles

5. **GitHub Pages**
   - Build frontend to `docs/`
   - Copy `scenarios/` bundles into `docs/scenarios/`
   - Publish

---

## Notes / constraints

- PNG overlays keep the UI simple and fast. Storage grows with:
  - number of scenarios × timesteps × layers
- GeoTIFFs for *every* timestep may get large quickly; keep them:
  - only for selected layers
  - only for selected timesteps
  - or switch to Zarr later
- AOI-bounding in MapLibre + trimmed overlays avoids wasted compute/rendering.
