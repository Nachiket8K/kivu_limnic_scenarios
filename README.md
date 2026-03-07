# Kivu Limnic Scenarios  
**Scenario discovery + browser playback for dense CO₂ cloud hazards from a potential limnic eruption at Lake Kivu**

Live viewer (GitHub Pages): `https://nachiket8k.github.io/kivu_limnic_scenarios/`

## Why this project exists

**Lake Kivu** (Rwanda / DRC) contains very large amounts of dissolved gases (CO₂ and methane). A rare but catastrophic hazard associated with deep, gas-charged lakes is a **limnic eruption** (also called a lake overturn): a sudden exsolution and release of dissolved gas that forms a ground-hugging, oxygen-displacing CO₂ cloud.

This project was built as a **scenario discovery tool** to explore:

- **What fraction of a theoretical upper bound** of stored CO₂ could become **fatal** if released  
- **How the hazard footprint evolves in space and time** (each timestep as a “video frame”)  
- **How many people might be exposed** above defined concentration thresholds (danger / lethal) over time  

The key design choice: **treat the lake’s theoretical upper-limit CO₂ inventory as a reference “capacity”** and run scenarios as **percentages of that upper limit** (e.g., 0.1%, 0.5%, 1%, …). This makes the scenario space interpretable: it answers “how bad does it get if X% of the maximum plausible CO₂ inventory is released?”

## Limnic eruptions in plain language

A **limnic eruption** happens when deep water is charged with dissolved gas under high pressure (Henry’s law), and a trigger (landslide, earthquake, mixing event, etc.) lifts gas-rich water upward. As pressure drops, gas comes out of solution, bubbles form, and buoyancy drives a **runaway** upward flow—like opening a shaken soda bottle, but on a lake scale.

Why it’s deadly: **CO₂ is denser than air**, so it can flow downhill into valleys, displace oxygen, and cause asphyxiation.

Historical events in Cameroon demonstrate the mechanism and consequences:
- **Lake Monoun (1984):** dozens of fatalities reported (commonly ~37).  
- **Lake Nyos (1986):** ~1,700 fatalities and thousands of livestock deaths; a CO₂ cloud travelled into surrounding low-lying areas.

## What this model simulates (and what it does not)

### Included (current scope)
This code focuses on **CO₂ release and atmospheric/near-ground dispersion**, producing:
- time-stepped rasters (exported as PNG frames for fast playback)
- time series metrics (CSV) capturing hazard dynamics and exposed population

The physics are implemented as a **two-layer approximation**:

1) **Dense base layer (TWODEE)**  
   A shallow-layer, terrain-following dense gas flow that “settles” and spreads along topography (valley-seeking behavior).

2) **Wind-advected top layer (DISGAS)**  
   A more dilute component that is transported primarily by wind and turbulent diffusion, with coupling/decay terms.

This is intentionally built for **rapid what-if exploration**, not as a full CFD/LES solver.

### Not included (yet)
- full lake overturn / gas exsolution dynamics inside the lake  
- bathymetry-driven lake physics (explicitly deferred in this project phase)  
- ERA5-driven time-varying winds (optional future upgrade)  
- chemistry beyond CO₂ (e.g., methane ignition scenarios)  

## Academic foundations for the math

This implementation is a pragmatic “scenario discovery” approximation, but it is explicitly inspired by established dense-gas dispersion modeling approaches:

### Dense-gas, terrain-following flow: TWODEE
- **Folch et al. (2017)** present high-resolution dense gas dispersion modelling using **TWODEE-2.1**, with application to the **1986 Lake Nyos** disaster. This paper is the conceptual anchor for treating dense gas as a terrain-following shallow layer.

### Passive dispersion (wind + diffusion): DISGAS
- **DISGAS-2.0** is an Eulerian model for passive dispersion of gas and particles, describing wind field options and turbulent diffusion (K-theory), used widely in volcanic gas dispersion workflows.

### Lake Kivu risk context (background literature)
For the broader Lake Kivu hazard context (not necessarily the numerical scheme):
- **Bärenbold et al. (2020)**: discussion of gas concentrations and limnic eruption risk at Lake Kivu.  
- **Saboorian‑Jooybari & Hassanzadeh (2026)**: “On the risk of a dissolved gas-triggered limnic eruption in Lake Kivu.”  
- **Nature (2021)** feature article discussing Lake Kivu’s gas system and hazards.

## Repository structure (high level)

- `notebooks/` — step-by-step workflow (AOI prep, alignment, proof-of-concept runs)  
- `src/kivu/` — Python package structure (I/O, preprocessing, model, viz)  
- `src/kivu/model/sim_run.py` — offline scenario runner exporting web bundles  
- `docs/` — GitHub Pages site (static viewer + published `docs/scenarios/...`)  
- `web/` — local web assets for development  
- `osm/` — OSM context extraction outputs (roads, water, buildings)  
- `data/` — raster inputs + processed aligned grids (kept mostly out of git due to size)  

## Quickstart

### 1) Create environment
Conda (recommended):
```bash
conda env create -f environment.yml
conda activate kivu_co2
```

### 2) Prepare data (AOI + alignment)
Use the notebooks in `notebooks/` to:
- define AOI  
- align DEM + population rasters to a common grid  
- create masks (lake / water / land)  
- compute slopes/gradients needed by the solver  

### 3) Run a scenario and export a bundle
The offline runner exports:
- `manifest.json`  
- `metrics_10min.csv`  
- PNG frames per layer (for browser playback)

### 4) View results in the browser
- Local development: serve `docs/` with a static server:
```bash
python -m http.server 8000 --directory docs
```
- GitHub Pages: published from `docs/`

## Scenario discovery philosophy

This repo is designed around a clean separation:
- **Compute offline once** (expensive physics)  
- **Render cheaply in browser** (fast playback, repeatable comparison)  

That makes it feasible to:
- precompute a scenario library (different release fractions, wind fields)  
- select scenarios from dropdowns (parameter-driven)  
- compare impact metrics quickly without recomputing PDEs  


## References (key)
- Folch, A., Barcons, J., Kozono, T., & Costa, A. (2017). *High-resolution modelling of atmospheric dispersion of dense gas using TWODEE‑2.1: application to the 1986 Lake Nyos limnic eruption.* **Natural Hazards and Earth System Sciences**, 17, 861–879.  
- INGV / Datasim. *DISGAS-2.0: A model for passive DISpersion of GAS (User Manual).*  
- Bärenbold, F. (2020). *No increasing risk of a limnic eruption at Lake Kivu…* (PMC).  
- Saboorian‑Jooybari, H., & Hassanzadeh, H. (2026). *On the risk of a dissolved gas-triggered limnic eruption in Lake Kivu.* **Environmental Science: Processes & Impacts**.  
- NOAA NGDC hazard entry (Lake Nyos event).  
- Britannica overview of Lake Nyos disaster.  
- USGS media note on “Exploding lakes” and degassing context.  
