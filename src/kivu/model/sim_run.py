#!/usr/bin/env python3
"""
sim_run.py — Lake Kivu limnic eruption scenario runner (export bundles for web viewer)

What it does
- Loads processed AOI grid + masks from data/processed/static_layers.json
- Runs the "Better Physics" model (TWODEE-like dense base layer + DISGAS-like wind top layer)
- Exports a scenario bundle:
  scenarios/scenario_<id>/
    manifest.json
    metrics_10min.csv
    frames/<layer>/t####.png
    raw/<optional GeoTIFFs>

Design notes
- Intended for *offline precomputation*: browser viewer plays back PNG frames (fast) and reads metrics CSV.
- Keep compute separate from rendering: no physics in the browser.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import xy as rc_to_xy
from scipy.special import erf as ERF

# Headless matplotlib for PNG export
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyproj import Transformer


# -----------------------------
# Utilities
# -----------------------------
def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def stable_hash(obj: Dict[str, Any], n: int = 8) -> str:
    """Deterministic short id from a dict (sorted keys)."""
    s = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha1(s.encode("utf-8")).hexdigest()[:n]


def wind_to_uv(speed_mps: float, dir_from_deg: float) -> Tuple[float, float]:
    """
    Meteorological wind direction (coming FROM):
    0 = from North, 90 = from East.
    Convert to (u,v) where +x east, +y north.
    """
    dir_to = (dir_from_deg + 180.0) % 360.0
    rad = np.deg2rad(dir_to)
    u = speed_mps * np.sin(rad)   # +x east
    v = speed_mps * np.cos(rad)   # +y north
    return float(u), float(v)


# -----------------------------
# Raster loading
# -----------------------------
def load_raster(path: Path) -> Tuple[rasterio.io.DatasetReader, np.ndarray]:
    ds = rasterio.open(path)
    arr = ds.read(1)
    return ds, arr


def close_datasets(dss: Iterable[rasterio.io.DatasetReader]) -> None:
    for ds in dss:
        try:
            ds.close()
        except Exception:
            pass


# -----------------------------
# Model helpers (dense base + wind top)
# -----------------------------
H_FLOOR = 1e-4  # meters; below this we treat as dry


def primitives(h: np.ndarray, mx: np.ndarray, my: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Dry-cell safe velocities."""
    h = h.astype("float32")
    mx = mx.astype("float32")
    my = my.astype("float32")
    u = np.zeros_like(h, dtype="float32")
    v = np.zeros_like(h, dtype="float32")
    wet = h > H_FLOOR
    u[wet] = mx[wet] / h[wet]
    v[wet] = my[wet] / h[wet]
    return u, v


def enforce_positivity_and_consistency(
    h: np.ndarray, mx: np.ndarray, my: np.ndarray, M_base: np.ndarray, rho_co2: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Enforce:
    - h >= 0
    - M_base >= 0
    - M_base <= rho_co2 * h (cannot exceed pure CO2 in that layer)
    - dry cells have mx=my=0
    """
    h = np.maximum(h, 0.0).astype("float32")
    M_base = np.maximum(M_base, 0.0).astype("float32")
    cap = (rho_co2 * h).astype("float32")
    M_base = np.minimum(M_base, cap).astype("float32")
    dry = h <= H_FLOOR
    mx = mx.astype("float32")
    my = my.astype("float32")
    mx[dry] = 0.0
    my[dry] = 0.0
    return h, mx, my, M_base


def co2_massfrac(M_base: np.ndarray, h: np.ndarray, rho_co2: float) -> np.ndarray:
    denom = (rho_co2 * h + 1e-12)
    f = np.clip(M_base / denom, 0.0, 1.0).astype("float32")
    return f


def rho_bar_from_frac(f: np.ndarray, rho_air: float, rho_co2: float) -> np.ndarray:
    return (rho_air + f * (rho_co2 - rho_air)).astype("float32")


def compute_gp(M_base: np.ndarray, h: np.ndarray, rho_air: float, rho_co2: float, g: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reduced gravity field g' = g*(rho_bar - rho_air)/rho_air, rho_bar from local CO2 fraction.
    Clamped to >=0.
    """
    f = np.zeros_like(h, dtype="float32")
    wet = h > H_FLOOR
    f[wet] = np.clip(M_base[wet] / (rho_co2 * h[wet] + 1e-12), 0.0, 1.0)
    rho_bar = (rho_air + f * (rho_co2 - rho_air)).astype("float32")
    gp = (g * (rho_bar - rho_air) / (rho_air + 1e-12)).astype("float32")
    gp = np.maximum(gp, 0.0).astype("float32")
    return gp, rho_bar, f


def flux_x_state(h: np.ndarray, mx: np.ndarray, my: np.ndarray, M: np.ndarray, gp: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u, v = primitives(h, mx, my)
    p = 0.5 * gp * h * h
    Fh  = mx
    Fmx = mx * u + p
    Fmy = my * u
    FM  = M * u
    return Fh, Fmx, Fmy, FM, u


def flux_y_state(h: np.ndarray, mx: np.ndarray, my: np.ndarray, M: np.ndarray, gp: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u, v = primitives(h, mx, my)
    p = 0.5 * gp * h * h
    Gh  = my
    Gmx = mx * v
    Gmy = my * v + p
    GM  = M * v
    return Gh, Gmx, Gmy, GM, v


def rusanov_interface(F_L: np.ndarray, F_R: np.ndarray, Q_L: np.ndarray, Q_R: np.ndarray, amax: np.ndarray) -> np.ndarray:
    return 0.5 * (F_L + F_R) - 0.5 * amax * (Q_R - Q_L)


def diffuse(F: np.ndarray, K: float, dt: float, dx: float, dy: float) -> np.ndarray:
    if K <= 0:
        return F
    F = F.astype("float32")
    L = np.pad(F, ((0,0),(1,0)), mode="edge")[:, :-1]
    R = np.pad(F, ((0,0),(0,1)), mode="edge")[:, 1:]
    U = np.pad(F, ((1,0),(0,0)), mode="edge")[:-1, :]
    D = np.pad(F, ((0,1),(0,0)), mode="edge")[1:, :]
    lap = (L - 2*F + R) / (dx*dx) + (U - 2*F + D) / (dy*dy)
    return (F + (K * dt) * lap).astype("float32")


def apply_friction(mx: np.ndarray, my: np.ndarray, h: np.ndarray, Cf: float, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    u, v = primitives(h, mx, my)
    speed = np.sqrt(u*u + v*v).astype("float32")
    mx2 = mx - (Cf * u * speed * h * dt).astype("float32")
    my2 = my - (Cf * v * speed * h * dt).astype("float32")
    return mx2.astype("float32"), my2.astype("float32")


def shiftR(A: np.ndarray) -> np.ndarray:
    return np.pad(A, ((0,0),(0,1)), mode="edge")[:, 1:]


def shiftL(A: np.ndarray) -> np.ndarray:
    return np.pad(A, ((0,0),(1,0)), mode="edge")[:, :-1]


def shiftD(A: np.ndarray) -> np.ndarray:
    return np.pad(A, ((0,1),(0,0)), mode="edge")[1:, :]


def shiftU(A: np.ndarray) -> np.ndarray:
    return np.pad(A, ((1,0),(0,0)), mode="edge")[:-1, :]


def twodee_step(
    h: np.ndarray,
    mx: np.ndarray,
    my: np.ndarray,
    M_base: np.ndarray,
    dt: float,
    dx: float,
    dy: float,
    rho_air: float,
    rho_co2: float,
    g: float,
    dz_dx: np.ndarray,
    dz_dy: np.ndarray,
    Cf: float,
    K_h: float,
    K_m: float,
    K_M: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One explicit TWODEE-like step (Rusanov + topo forcing + friction + diffusion)."""
    h, mx, my, M_base = enforce_positivity_and_consistency(h, mx, my, M_base, rho_co2)
    gp, _, _ = compute_gp(M_base, h, rho_air, rho_co2, g)

    # X interfaces
    hR, mxR, myR, MR, gpR = shiftR(h), shiftR(mx), shiftR(my), shiftR(M_base), shiftR(gp)
    FhL, FmxL, FmyL, FML, uL = flux_x_state(h,  mx,  my,  M_base, gp)
    FhR, FmxR, FmyR, FMR, uR = flux_x_state(hR, mxR, myR, MR,     gpR)

    cL = np.sqrt(np.maximum(gp*h,   0.0)).astype("float32")
    cR = np.sqrt(np.maximum(gpR*hR, 0.0)).astype("float32")
    amax_x = np.maximum(np.abs(uL)+cL, np.abs(uR)+cR).astype("float32")
    amax_x = np.maximum(amax_x, 1e-6)

    Fx_h  = rusanov_interface(FhL,  FhR,  h,  hR,  amax_x)
    Fx_mx = rusanov_interface(FmxL, FmxR, mx, mxR, amax_x)
    Fx_my = rusanov_interface(FmyL, FmyR, my, myR, amax_x)
    Fx_M  = rusanov_interface(FML,  FMR,  M_base, MR, amax_x)

    h      = h      - (dt/dx) * (Fx_h  - shiftL(Fx_h))
    mx     = mx     - (dt/dx) * (Fx_mx - shiftL(Fx_mx))
    my     = my     - (dt/dx) * (Fx_my - shiftL(Fx_my))
    M_base = M_base - (dt/dx) * (Fx_M  - shiftL(Fx_M))

    # Y interfaces
    gp, _, _ = compute_gp(M_base, h, rho_air, rho_co2, g)
    hD, mxD, myD, MD, gpD = shiftD(h), shiftD(mx), shiftD(my), shiftD(M_base), shiftD(gp)
    GhL, GmxL, GmyL, GML, vL = flux_y_state(h,  mx,  my,  M_base, gp)
    GhD, GmxD, GmyD, GMD, vD = flux_y_state(hD, mxD, myD, MD,     gpD)

    cL = np.sqrt(np.maximum(gp*h,   0.0)).astype("float32")
    cD = np.sqrt(np.maximum(gpD*hD, 0.0)).astype("float32")
    amax_y = np.maximum(np.abs(vL)+cL, np.abs(vD)+cD).astype("float32")
    amax_y = np.maximum(amax_y, 1e-6)

    Gy_h  = rusanov_interface(GhL,  GhD,  h,  hD,  amax_y)
    Gy_mx = rusanov_interface(GmxL, GmxD, mx, mxD, amax_y)
    Gy_my = rusanov_interface(GmyL, GmyD, my, myD, amax_y)
    Gy_M  = rusanov_interface(GML,  GMD,  M_base, MD, amax_y)

    h      = h      - (dt/dy) * (Gy_h  - shiftU(Gy_h))
    mx     = mx     - (dt/dy) * (Gy_mx - shiftU(Gy_mx))
    my     = my     - (dt/dy) * (Gy_my - shiftU(Gy_my))
    M_base = M_base - (dt/dy) * (Gy_M  - shiftU(Gy_M))

    # Topography forcing and friction
    gp_now, _, _ = compute_gp(M_base, h, rho_air, rho_co2, g)
    mx = mx - (gp_now * h * dz_dx * dt).astype("float32")
    my = my - (gp_now * h * dz_dy * dt).astype("float32")
    mx, my = apply_friction(mx, my, h, Cf, dt)

    # Stabilizing diffusion
    h      = diffuse(h,      K_h, dt, dx, dy)
    mx     = diffuse(mx,     K_m, dt, dx, dy)
    my     = diffuse(my,     K_m, dt, dx, dy)
    M_base = diffuse(M_base, K_M, dt, dx, dy)

    h, mx, my, M_base = enforce_positivity_and_consistency(h, mx, my, M_base, rho_co2)
    return h, mx, my, M_base


def advect_upwind_scalar(F: np.ndarray, u: float, v: float, dt: float, dx: float, dy: float) -> np.ndarray:
    """Simple upwind advection for the top-layer mass (DISGAS-like)."""
    F = F.astype("float32")
    L = np.pad(F, ((0,0),(1,0)), mode="edge")[:, :-1]
    R = np.pad(F, ((0,0),(0,1)), mode="edge")[:, 1:]
    U = np.pad(F, ((1,0),(0,0)), mode="edge")[:-1, :]
    D = np.pad(F, ((0,1),(0,0)), mode="edge")[1:, :]
    dFdx = np.where(u >= 0, (F - L)/dx, (R - F)/dx)
    dFdy = np.where(v >= 0, (F - U)/dy, (D - F)/dy)
    return (F - dt*(u*dFdx + v*dFdy)).astype("float32")


# -----------------------------
# Paper-grounded hazard mapping (Eq 4–5) + lethality (Eq 7–9)
# -----------------------------
def rho_at_height_z(h_eff: np.ndarray, rho_bar: np.ndarray, rho_air: float, S1: float, z: float) -> np.ndarray:
    h = np.maximum(h_eff, 1e-6).astype("float32")
    term = (2.0 / (S1 + 1e-12)) * (rho_bar - rho_air)
    expo = np.exp(-(2.0 / (S1 + 1e-12)) * (z / h)).astype("float32")
    rho_z = (rho_air + term * expo).astype("float32")
    rho_z = np.where(h_eff <= 0, rho_air, rho_z).astype("float32")
    return rho_z


def ppm_from_rho(rho_z: np.ndarray, rho_air: float, rho_co2: float, c_background_ppm: float) -> np.ndarray:
    frac = (rho_z - rho_air) / ((rho_co2 - rho_air) + 1e-12)
    frac = np.clip(frac, 0.0, 1.0).astype("float32")
    ppm = (c_background_ppm + (1e6 - c_background_ppm) * frac).astype("float32")
    return ppm


def prob_death_from_ppm(ppm: np.ndarray, exposure_min: np.ndarray, pars: Dict[str, Any]) -> np.ndarray:
    c_percent = (ppm / 1e4).astype("float32")
    d_min = exposure_min.astype("float32")
    mu = pars["fatal_a0"] + pars["fatal_b0"] / (1.0 + np.power(d_min, pars["fatal_c0"]))
    sig = pars["fatal_a1"] + pars["fatal_b1"] / (1.0 + np.power(d_min, pars["fatal_c1"]))
    z = (c_percent - mu) / (np.sqrt(2.0) * (sig + 1e-12))
    P = (0.5 * (1.0 + ERF(z))).astype("float32")
    return np.clip(P, 0.0, 1.0).astype("float32")


# -----------------------------
# Export: PNG frames + optional GeoTIFFs
# -----------------------------
def array_to_rgba_png(
    arr: np.ndarray,
    out_path: Path,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: str = "viridis",
    alpha_threshold: Optional[float] = None,
    dpi: int = 140,
) -> None:
    """
    Save a raster as a PNG image.
    - If alpha_threshold is set, values below it become transparent.
    - vmin/vmax control normalization; if None, computed from finite values.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    a = arr.astype("float32")
    finite = np.isfinite(a)
    if vmin is None:
        vmin = float(np.nanmin(a[finite])) if finite.any() else 0.0
    if vmax is None:
        vmax = float(np.nanmax(a[finite])) if finite.any() else 1.0
    if vmax <= vmin:
        vmax = vmin + 1e-6

    norm = (a - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0.0, 1.0)

    cm = plt.get_cmap(cmap)
    rgba = cm(norm)  # (H,W,4) float 0..1

    if alpha_threshold is not None:
        mask = np.where(np.isfinite(a), a, -np.inf) < alpha_threshold
        rgba[..., 3] = np.where(mask, 0.0, rgba[..., 3])

    plt.figure(figsize=(6, 6), dpi=dpi)
    plt.axis("off")
    plt.imshow(rgba, interpolation="nearest")
    plt.tight_layout(pad=0)
    plt.savefig(out_path, bbox_inches="tight", pad_inches=0)
    plt.close()


def write_geotiff_like(template_ds: rasterio.io.DatasetReader, out_path: Path, array: np.ndarray, nodata: Optional[float] = None) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    meta = template_ds.meta.copy()
    meta.update(count=1, dtype="float32", nodata=nodata, compress="deflate")
    with rasterio.open(out_path, "w", **meta) as dst:
        dst.write(array.astype("float32"), 1)


# -----------------------------
# Scenario export: manifest + index
# -----------------------------
def utm_bounds_to_wgs84(bounds: Tuple[float, float, float, float], crs: str) -> Dict[str, Any]:
    """
    bounds: (left, bottom, right, top) in projected CRS.
    returns corners and bbox in WGS84.
    """
    transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    left, bottom, right, top = bounds
    # corners in order: top-left, top-right, bottom-right, bottom-left (MapLibre ImageSource expects 4 corners)
    tl = transformer.transform(left, top)
    tr = transformer.transform(right, top)
    br = transformer.transform(right, bottom)
    bl = transformer.transform(left, bottom)
    xs = [tl[0], tr[0], br[0], bl[0]]
    ys = [tl[1], tr[1], br[1], bl[1]]
    return {
        "corners_lonlat": [list(tl), list(tr), list(br), list(bl)],
        "bbox_lonlat": [min(xs), min(ys), max(xs), max(ys)],
    }


def write_manifest(
    scenario_dir: Path,
    scenario_id: str,
    static: Dict[str, Any],
    scenario: Dict[str, Any],
    grid_ds: rasterio.io.DatasetReader,
    n_frames: int,
    dt_report_min: int,
    layers: List[Dict[str, Any]],
) -> None:
    wgs84 = utm_bounds_to_wgs84(grid_ds.bounds, str(grid_ds.crs))
    manifest = {
        "scenario_id": scenario_id,
        "created_utc": _dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "aoi": {
            "crs": str(grid_ds.crs),
            "bounds_projected": [grid_ds.bounds.left, grid_ds.bounds.bottom, grid_ds.bounds.right, grid_ds.bounds.top],
            "bbox_lonlat": wgs84["bbox_lonlat"],
            "corners_lonlat": wgs84["corners_lonlat"],
            "maxBounds_lonlat": wgs84["bbox_lonlat"],
        },
        "grid": {
            "width": int(grid_ds.width),
            "height": int(grid_ds.height),
            "dx_m": float(static["dx_m"]),
            "dy_m": float(static["dy_m"]),
            "transform": list(grid_ds.transform),
        },
        "time": {
            "dt_report_min": int(dt_report_min),
            "n_frames": int(n_frames),
            "timestamps_min": [int(i * dt_report_min) for i in range(n_frames)],
        },
        "parameters": scenario,
        "assets": {
            "metrics": "metrics_10min.csv",
        },
        "layers": layers,
    }
    (scenario_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def update_index(index_path: Path, scenario_id: str, scenario: Dict[str, Any]) -> None:
    """
    Maintain scenarios/index.json with scenario list and available dropdown options.
    """
    index_path.parent.mkdir(parents=True, exist_ok=True)
    if index_path.exists():
        idx = json.loads(index_path.read_text(encoding="utf-8"))
    else:
        idx = {"scenarios": []}

    # Remove existing entry if present
    idx["scenarios"] = [s for s in idx.get("scenarios", []) if s.get("scenario_id") != scenario_id]
    idx["scenarios"].append({
        "scenario_id": scenario_id,
        "path": f"scenario_{scenario_id}/manifest.json",
        "parameters": scenario,
    })

    # Derive dropdown options from all scenarios present
    opts: Dict[str, List[Any]] = {}
    for s in idx["scenarios"]:
        for k, v in s.get("parameters", {}).items():
            # keep only JSON-serializable scalar-ish values for dropdowns
            if isinstance(v, (int, float, str, bool)):
                opts.setdefault(k, set()).add(v)
    idx["options"] = {k: sorted(list(v)) for k, v in opts.items()}

    index_path.write_text(json.dumps(idx, indent=2), encoding="utf-8")


# -----------------------------
# Disk init
# -----------------------------
def disk_mask(H: int, W: int, center_r: int, center_c: int, radius_m: float, dx: float, dy: float, mask_limit: Optional[np.ndarray] = None) -> np.ndarray:
    rr = np.arange(H)[:, None]
    cc = np.arange(W)[None, :]
    dr_m = (rr - center_r) * dy
    dc_m = (cc - center_c) * dx
    disk = ((dr_m**2 + dc_m**2) <= radius_m**2).astype(np.uint8)
    if mask_limit is not None:
        disk = (disk & mask_limit.astype(np.uint8)).astype(np.uint8)
    return disk


# -----------------------------
# Main simulation runner
# -----------------------------
def run_single_scenario(
    static: Dict[str, Any],
    scenario_in: Dict[str, Any],
    out_root: Path,
    save_geotiff: bool,
    save_frames: bool,
) -> Path:
    """
    Runs one scenario and exports a scenario bundle into out_root/scenario_<id>/
    Returns the scenario directory.
    """
    # Load rasters
    paths = static["paths"]
    dem_ds, dem = load_raster(Path(paths["dem"]))
    pop_ds, pop = load_raster(Path(paths["pop"]))
    mask_lake_ds, mask_lake = load_raster(Path(paths["mask_lake"]))
    grid_ds, _ = load_raster(Path(paths["grid_template"]))

    source_idx = json.loads(Path(paths["source_index"]).read_text(encoding="utf-8"))
    src_r, src_c = int(source_idx["row"]), int(source_idx["col"])

    H, W = dem.shape
    dx = float(static["dx_m"]); dy = float(static["dy_m"])
    cell_area = dx * dy

    # Terrain gradients
    dz_dy, dz_dx = np.gradient(dem.astype("float32"), dy, dx)

    # Scenario defaults + derived fields
    scenario = dict(scenario_in)  # copy
    scenario.setdefault("duration_h", 12.0)
    scenario.setdefault("dt_internal_s", 0.3)
    scenario.setdefault("dt_report_min", 10)
    scenario.setdefault("save_images_every_min", 60)
    scenario.setdefault("exposure_substep_min", 3)

    scenario.setdefault("co2_inventory_upper_m3_stp", 3.0e11)
    scenario.setdefault("release_fraction", 0.01)
    scenario.setdefault("eruption_radius_m", 3000.0)
    scenario.setdefault("constrain_disk_to_lake", True)

    scenario.setdefault("wind_speed_mps", 2.0)
    scenario.setdefault("wind_dir_deg_from", 90.0)

    scenario.setdefault("g", 9.81)
    scenario.setdefault("rho_air_ref", 1.20)
    scenario.setdefault("rho_co2_ref", 1.84)

    scenario.setdefault("Cf_friction", 0.002)
    scenario.setdefault("K_h_m2ps", 0.5)
    scenario.setdefault("K_M_m2ps", 0.5)
    scenario.setdefault("K_m_m2ps", 0.2)

    scenario.setdefault("K_top_m2ps", 20.0)
    scenario.setdefault("strip_coeff_1ps_per_mps", 2.0e-4)
    scenario.setdefault("top_background_loss_1ps", 0.0)

    scenario.setdefault("base_background_loss_1ps", 2.0e-6)
    scenario.setdefault("entrain_rate_1ps", 2.0e-5)
    scenario.setdefault("entrain_wind_gain", 1.0)

    scenario.setdefault("S1", 0.5)
    scenario.setdefault("c_background_ppm", 420.0)
    scenario.setdefault("z_human_m", 1.5)

    scenario.setdefault("fatal_a0", 5.056); scenario.setdefault("fatal_b0", 17.885); scenario.setdefault("fatal_c0", 0.357)
    scenario.setdefault("fatal_a1", 0.662); scenario.setdefault("fatal_b1", 2.421);  scenario.setdefault("fatal_c1", 0.354)

    scenario.setdefault("exposure_threshold_ppm", 30000.0)

    scenario["released_volume_STP_m3"] = float(scenario["co2_inventory_upper_m3_stp"]) * float(scenario["release_fraction"])

    # Determine scenario id from "UI-visible" parameters
    id_basis = {
        "release_fraction": scenario["release_fraction"],
        "eruption_radius_m": scenario["eruption_radius_m"],
        "wind_speed_mps": scenario["wind_speed_mps"],
        "wind_dir_deg_from": scenario["wind_dir_deg_from"],
        "duration_h": scenario["duration_h"],
        "dt_report_min": scenario["dt_report_min"],
    }
    scenario_id = stable_hash(id_basis)

    scenario_dir = out_root / f"scenario_{scenario_id}"
    ensure_dir(scenario_dir)
    ensure_dir(scenario_dir / "frames")
    ensure_dir(scenario_dir / "raw")

    # Write scenario config used
    (scenario_dir / "scenario.json").write_text(json.dumps(scenario, indent=2), encoding="utf-8")

    # Timing
    duration_s = int(round(float(scenario["duration_h"]) * 3600))
    dt = float(scenario["dt_internal_s"])

    dt_report_s = int(round(float(scenario["dt_report_min"]) * 60))
    save_img_s = int(round(float(scenario["save_images_every_min"]) * 60))
    sub_s = int(round(float(scenario["exposure_substep_min"]) * 60))

    n_steps = int(np.ceil(duration_s / dt))
    report_every = max(1, int(round(dt_report_s / dt)))
    save_every = max(1, int(round(save_img_s / dt)))
    sub_every = max(1, int(round(sub_s / dt)))

    # Wind
    u_wind, v_wind = wind_to_uv(float(scenario["wind_speed_mps"]), float(scenario["wind_dir_deg_from"]))
    wind_mag = float(np.hypot(u_wind, v_wind))

    # Init disk release
    limit_mask = mask_lake if bool(scenario["constrain_disk_to_lake"]) else None
    disk = disk_mask(H, W, src_r, src_c, float(scenario["eruption_radius_m"]), dx, dy, mask_limit=limit_mask)
    disk_px = int(disk.sum())
    if disk_px <= 0:
        close_datasets([dem_ds, pop_ds, mask_lake_ds, grid_ds])
        raise ValueError("Eruption disk is empty. Increase radius or disable constrain_disk_to_lake.")

    released_V = float(scenario["released_volume_STP_m3"])
    rho_co2 = float(scenario["rho_co2_ref"])
    rho_air = float(scenario["rho_air_ref"])
    g = float(scenario["g"])

    M_total_kg = released_V * rho_co2
    m_area = M_total_kg / (disk_px * cell_area + 1e-12)
    h0 = m_area / (rho_co2 + 1e-12)

    h = np.zeros((H, W), dtype="float32")
    mx = np.zeros((H, W), dtype="float32")
    my = np.zeros((H, W), dtype="float32")
    M_base = np.zeros((H, W), dtype="float32")
    M_top = np.zeros((H, W), dtype="float32")
    exposure_min = np.zeros((H, W), dtype="float32")
    Pdeath = np.zeros((H, W), dtype="float32")

    h[disk == 1] = h0
    M_base[disk == 1] = m_area
    h, mx, my, M_base = enforce_positivity_and_consistency(h, mx, my, M_base, rho_co2)

    # Metrics
    THRESH_DANGER = 50_000.0
    THRESH_LETHAL = 100_000.0
    pop_safe = np.where(np.isfinite(pop.astype("float32")), pop.astype("float32"), 0.0)

    metrics_rows: List[Dict[str, Any]] = []

    # Frame export config
    layer_defs = [
        {"name": "ppm_z", "type": "png", "pathTemplate": "frames/ppm_z/t{frame:04d}.png", "opacityDefault": 0.7},
        {"name": "lethal_mask", "type": "png", "pathTemplate": "frames/lethal_mask/t{frame:04d}.png", "opacityDefault": 0.5},
        {"name": "pdeath", "type": "png", "pathTemplate": "frames/pdeath/t{frame:04d}.png", "opacityDefault": 0.6},
        {"name": "base_h", "type": "png", "pathTemplate": "frames/base_h/t{frame:04d}.png", "opacityDefault": 0.5},
    ]

    # Prepare frame directories
    if save_frames:
        for ld in layer_defs:
            ensure_dir(scenario_dir / ld["pathTemplate"].split("/")[0] / ld["name"])

    # cumulative accounting (optional)
    cum_entr_kg = 0.0
    cum_strip_kg = 0.0
    cum_base_loss_kg = 0.0

    # Run loop
    ppm: Optional[np.ndarray] = None
    frame_idx = 0

    for step in range(n_steps + 1):
        # Base layer step
        h, mx, my, M_base = twodee_step(
            h, mx, my, M_base,
            dt, dx, dy,
            rho_air, rho_co2, g,
            dz_dx.astype("float32"), dz_dy.astype("float32"),
            float(scenario["Cf_friction"]),
            float(scenario["K_h_m2ps"]),
            float(scenario["K_m_m2ps"]),
            float(scenario["K_M_m2ps"]),
        )

        # Entrain base -> top
        entrain_eff = float(scenario["entrain_rate_1ps"]) * (1.0 + float(scenario["entrain_wind_gain"]) * wind_mag)
        dM = (entrain_eff * dt) * M_base
        dM = np.minimum(dM, M_base)
        M_base -= dM
        M_top += dM
        cum_entr_kg += float(np.sum(dM) * cell_area)
        h, mx, my, M_base = enforce_positivity_and_consistency(h, mx, my, M_base, rho_co2)

        # Base background loss
        base_loss = float(scenario["base_background_loss_1ps"])
        if base_loss > 0:
            dM_loss = (1.0 - np.exp(-base_loss * dt)) * M_base
            M_base -= dM_loss
            cum_base_loss_kg += float(np.sum(dM_loss) * cell_area)
            h, mx, my, M_base = enforce_positivity_and_consistency(h, mx, my, M_base, rho_co2)

        # Top layer: wind advection + diffusion + stripping
        M_top = advect_upwind_scalar(M_top, u_wind, v_wind, dt, dx, dy)
        M_top = diffuse(M_top, float(scenario["K_top_m2ps"]), dt, dx, dy)
        M_top = np.maximum(M_top, 0).astype("float32")

        strip_rate = float(scenario["strip_coeff_1ps_per_mps"]) * wind_mag + float(scenario["top_background_loss_1ps"])
        if strip_rate > 0:
            before = float(np.sum(M_top) * cell_area)
            M_top = (M_top * np.exp(-strip_rate * dt)).astype("float32")
            after = float(np.sum(M_top) * cell_area)
            cum_strip_kg += max(0.0, before - after)

        # Hazard substep
        if step % sub_every == 0:
            f_base = co2_massfrac(M_base, h, rho_co2)
            rho_bar_base = rho_bar_from_frac(f_base, rho_air, rho_co2)
            h_top_equiv = (M_top / (rho_co2 + 1e-12)).astype("float32")
            h_eff = (h + h_top_equiv).astype("float32")

            rho_z = rho_at_height_z(h_eff, rho_bar_base, rho_air, float(scenario["S1"]), float(scenario["z_human_m"]))
            ppm = ppm_from_rho(rho_z, rho_air, rho_co2, float(scenario["c_background_ppm"]))

            exposure_min += (ppm >= float(scenario["exposure_threshold_ppm"])).astype("float32") * float(scenario["exposure_substep_min"])
            Pdeath = prob_death_from_ppm(ppm, exposure_min, scenario)

        # Reporting + exports
        if step % report_every == 0:
            # Ensure ppm exists (if report_every < sub_every)
            if ppm is None:
                f_base = co2_massfrac(M_base, h, rho_co2)
                rho_bar_base = rho_bar_from_frac(f_base, rho_air, rho_co2)
                h_top_equiv = (M_top / (rho_co2 + 1e-12)).astype("float32")
                h_eff = (h + h_top_equiv).astype("float32")
                rho_z = rho_at_height_z(h_eff, rho_bar_base, rho_air, float(scenario["S1"]), float(scenario["z_human_m"]))
                ppm = ppm_from_rho(rho_z, rho_air, rho_co2, float(scenario["c_background_ppm"]))
                Pdeath = prob_death_from_ppm(ppm, exposure_min, scenario)

            danger = (ppm >= THRESH_DANGER)
            lethal = (ppm >= THRESH_LETHAL)

            area_danger_km2 = float(danger.sum() * cell_area / 1e6)
            area_lethal_km2 = float(lethal.sum() * cell_area / 1e6)
            pop_danger = float(np.sum(pop_safe[danger]))
            pop_lethal = float(np.sum(pop_safe[lethal]))
            expected_fatalities = float(np.sum(pop_safe * Pdeath))

            u, v = primitives(h, mx, my)
            vmax_flow = float(np.nanmax(np.sqrt(u*u + v*v)))

            row = {
                "time_min": float(frame_idx * float(scenario["dt_report_min"])),
                "time_hr": float(frame_idx * float(scenario["dt_report_min"]) / 60.0),
                "vol_base_geom_m3": float(np.sum(h) * cell_area),
                "mass_base_kg": float(np.sum(M_base) * cell_area),
                "mass_top_kg": float(np.sum(M_top) * cell_area),
                "cum_entrained_kg": float(cum_entr_kg),
                "cum_stripped_kg": float(cum_strip_kg),
                "cum_base_loss_kg": float(cum_base_loss_kg),
                "max_ppm_z": float(np.nanmax(ppm)),
                "p95_ppm_z": float(np.nanpercentile(ppm[np.isfinite(ppm)], 95)) if np.isfinite(ppm).any() else np.nan,
                "area_danger_km2": area_danger_km2,
                "area_lethal_km2": area_lethal_km2,
                "pop_danger": pop_danger,
                "pop_lethal": pop_lethal,
                "expected_fatalities": expected_fatalities,
                "max_flow_speed_mps": vmax_flow,
                "mean_exposure_min": float(np.mean(exposure_min)),
                "max_exposure_min": float(np.max(exposure_min)),
                "max_Pdeath": float(np.max(Pdeath)),
            }
            metrics_rows.append(row)

            # Export frames
            if save_frames:
                # ppm overlay (transparent under threshold)
                array_to_rgba_png(
                    ppm,
                    scenario_dir / "frames" / "ppm_z" / f"t{frame_idx:04d}.png",
                    vmin=0.0, vmax=200_000.0, cmap="magma",
                    alpha_threshold=1_000.0,  # hide near-background
                )
                lethal_png = lethal.astype("float32")
                array_to_rgba_png(
                    lethal_png,
                    scenario_dir / "frames" / "lethal_mask" / f"t{frame_idx:04d}.png",
                    vmin=0.0, vmax=1.0, cmap="Reds",
                    alpha_threshold=0.5,
                )
                array_to_rgba_png(
                    Pdeath,
                    scenario_dir / "frames" / "pdeath" / f"t{frame_idx:04d}.png",
                    vmin=0.0, vmax=1.0, cmap="inferno",
                    alpha_threshold=0.01,
                )
                array_to_rgba_png(
                    h,
                    scenario_dir / "frames" / "base_h" / f"t{frame_idx:04d}.png",
                    vmin=0.0, vmax=max(float(np.nanmax(h)), 1e-6), cmap="viridis",
                    alpha_threshold=H_FLOOR,
                )

            # Optional raw GeoTIFFs at report times (can grow large)
            if save_geotiff:
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"h_t{frame_idx:04d}.tif", h, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"M_base_t{frame_idx:04d}.tif", M_base, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"M_top_t{frame_idx:04d}.tif", M_top, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"ppm_t{frame_idx:04d}.tif", ppm, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"Pdeath_t{frame_idx:04d}.tif", Pdeath, nodata=0.0)

            frame_idx += 1

    # Save metrics
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(scenario_dir / "metrics_10min.csv", index=False)

    # Write manifest + update index
    write_manifest(
        scenario_dir=scenario_dir,
        scenario_id=scenario_id,
        static=static,
        scenario=id_basis | {"released_volume_STP_m3": scenario["released_volume_STP_m3"]},
        grid_ds=grid_ds,
        n_frames=frame_idx,
        dt_report_min=int(scenario["dt_report_min"]),
        layers=layer_defs,
    )
    update_index(out_root / "index.json", scenario_id, id_basis | {"released_volume_STP_m3": scenario["released_volume_STP_m3"]})

    # Close datasets
    close_datasets([dem_ds, pop_ds, mask_lake_ds, grid_ds])
    return scenario_dir


# -----------------------------
# CLI
# -----------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Run Lake Kivu scenarios and export web-viewer bundles.")
    ap.add_argument("--scenario", type=str, default=None,
                    help="Path to a scenario JSON file (single scenario dict). If omitted, uses defaults.")
    ap.add_argument("--out", type=str, default="scenarios",
                    help="Output root directory for scenario bundles (default: scenarios).")
    ap.add_argument("--save-geotiff", action="store_true",
                    help="Also write GeoTIFFs at each report time (storage-heavy).")
    ap.add_argument("--no-frames", action="store_true",
                    help="Skip PNG frame export (metrics + manifest only).")
    return ap.parse_args()


def main() -> None:
    project_root = Path(".").resolve()
    proc_dir = project_root / "data" / "processed"
    static_path = proc_dir / "static_layers.json"
    if not static_path.exists():
        raise FileNotFoundError(f"Missing {static_path}. Run preprocessing notebooks first.")
    static = load_json(static_path)

    scenario_in: Dict[str, Any] = {}
    args = parse_args()
    if args.scenario:
        scenario_in = load_json(Path(args.scenario))

    out_root = Path(args.out).resolve()
    ensure_dir(out_root)

    scenario_dir = run_single_scenario(
        static=static,
        scenario_in=scenario_in,
        out_root=out_root,
        save_geotiff=bool(args.save_geotiff),
        save_frames=not bool(args.no_frames),
    )

    # Minimal success print (no debug spam)
    print(json.dumps({
        "status": "ok",
        "scenario_dir": str(scenario_dir),
        "manifest": str(scenario_dir / "manifest.json"),
        "metrics": str(scenario_dir / "metrics_10min.csv"),
    }, indent=2))


if __name__ == "__main__":
    main()
