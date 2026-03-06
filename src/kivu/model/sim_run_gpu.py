#!/usr/bin/env python3
"""
sim_run_numba.py — Lake Kivu limnic eruption scenario runner (Numba-accelerated)

- CPU acceleration (recommended for NVIDIA GeForce MX330): big speedups without GPU porting.
- Exports the same scenario bundle structure:
  scenarios/scenario_<id>/
    manifest.json
    metrics_10min.csv
    frames/<layer>/t####.png
    raw/<optional GeoTIFFs>

Run examples (Windows cmd):
  python sim_run_numba.py
  python sim_run_numba.py --scenario scenariospec.json
  python sim_run_numba.py --save-geotiff
  python sim_run_numba.py --no-frames
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import math
from pathlib import Path
from tqdm import tqdm
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import rasterio
from scipy.special import erf as ERF

# Headless matplotlib for PNG export
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyproj import Transformer

try:
    from numba import njit, prange
except Exception as e:
    raise RuntimeError(
        "Numba is required for this runner. Install with: pip install numba\n"
        f"Import error: {e!r}"
    )


# -----------------------------
# Utilities
# -----------------------------
def load_json(path: Path) -> Dict[str, Any]:
    print(f"Loading JSON from {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def stable_hash(obj: Dict[str, Any], n: int = 8) -> str:
    s = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha1(s.encode("utf-8")).hexdigest()[:n]


def wind_to_uv(speed_mps: float, dir_from_deg: float) -> Tuple[float, float]:
    dir_to = (dir_from_deg + 180.0) % 360.0
    rad = np.deg2rad(dir_to)
    u = speed_mps * np.sin(rad)   # +x east
    v = speed_mps * np.cos(rad)   # +y north
    print("Breaking Wind")
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
# Numba kernels
# -----------------------------
H_FLOOR = 1e-4  # meters

@njit(fastmath=True)
def _enforce_consistency_inplace(h, mx, my, M_base, rho_co2):
    H, W = h.shape
    for i in range(H):
        for j in range(W):
            if h[i, j] < 0.0:
                h[i, j] = 0.0
            if M_base[i, j] < 0.0:
                M_base[i, j] = 0.0
            cap = rho_co2 * h[i, j]
            if M_base[i, j] > cap:
                M_base[i, j] = cap
            if h[i, j] <= H_FLOOR:
                mx[i, j] = 0.0
                my[i, j] = 0.0


@njit(fastmath=True)
def _compute_gp_and_uv(h, mx, my, M_base, rho_air, rho_co2, g, gp, u, v):
    H, W = h.shape
    inv_rho_air = 1.0 / (rho_air + 1e-12)
    for i in range(H):
        for j in range(W):
            hij = h[i, j]
            if hij > H_FLOOR:
                uij = mx[i, j] / hij
                vij = my[i, j] / hij
                f = M_base[i, j] / (rho_co2 * hij + 1e-12)
                if f < 0.0:
                    f = 0.0
                if f > 1.0:
                    f = 1.0
                rho_bar = rho_air + f * (rho_co2 - rho_air)
                gpj = g * (rho_bar - rho_air) * inv_rho_air
                if gpj < 0.0:
                    gpj = 0.0
                gp[i, j] = gpj
                u[i, j] = uij
                v[i, j] = vij
            else:
                gp[i, j] = 0.0
                u[i, j] = 0.0
                v[i, j] = 0.0


@njit(fastmath=True, parallel=True)
def _diffuse_laplacian(F, K, dt, dx, dy, out):
    if K <= 0.0:
        out[:, :] = F
        return
    H, W = F.shape
    idx2 = 1.0 / (dx * dx)
    idy2 = 1.0 / (dy * dy)
    for i in prange(H):
        i_up = i - 1 if i > 0 else 0
        i_dn = i + 1 if i < H - 1 else H - 1
        for j in range(W):
            j_lt = j - 1 if j > 0 else 0
            j_rt = j + 1 if j < W - 1 else W - 1
            lap = (F[i, j_lt] - 2.0 * F[i, j] + F[i, j_rt]) * idx2 + (F[i_up, j] - 2.0 * F[i, j] + F[i_dn, j]) * idy2
            out[i, j] = F[i, j] + (K * dt) * lap


@njit(fastmath=True, parallel=True)
def _advect_upwind_scalar(F, u_wind, v_wind, dt, dx, dy, out):
    H, W = F.shape
    invdx = 1.0 / dx
    invdy = 1.0 / dy
    for i in prange(H):
        i_up = i - 1 if i > 0 else 0
        i_dn = i + 1 if i < H - 1 else H - 1
        for j in range(W):
            j_lt = j - 1 if j > 0 else 0
            j_rt = j + 1 if j < W - 1 else W - 1
            if u_wind >= 0.0:
                dFdx = (F[i, j] - F[i, j_lt]) * invdx
            else:
                dFdx = (F[i, j_rt] - F[i, j]) * invdx
            if v_wind >= 0.0:
                dFdy = (F[i, j] - F[i_up, j]) * invdy
            else:
                dFdy = (F[i_dn, j] - F[i, j]) * invdy
            out[i, j] = F[i, j] - dt * (u_wind * dFdx + v_wind * dFdy)


@njit(fastmath=True, parallel=True)
def _apply_friction(mx, my, h, Cf, dt, out_mx, out_my):
    H, W = h.shape
    for i in prange(H):
        for j in range(W):
            hij = h[i, j]
            if hij > H_FLOOR:
                u = mx[i, j] / hij
                v = my[i, j] / hij
                sp = math.sqrt(u*u + v*v)
                out_mx[i, j] = mx[i, j] - (Cf * u * sp * hij * dt)
                out_my[i, j] = my[i, j] - (Cf * v * sp * hij * dt)
            else:
                out_mx[i, j] = 0.0
                out_my[i, j] = 0.0


@njit(fastmath=True, parallel=True)
def _twodee_step_numba(h, mx, my, M_base, dz_dx, dz_dy,
                       dt, dx, dy, rho_air, rho_co2, g,
                       Cf, K_h, K_m, K_M,
                       gp, u, v,
                       h_new, mx_new, my_new, M_new,
                       tmp0, tmp1, tmp2, tmp3):
    _enforce_consistency_inplace(h, mx, my, M_base, rho_co2)
    _compute_gp_and_uv(h, mx, my, M_base, rho_air, rho_co2, g, gp, u, v)

    H, W = h.shape
    invdx = 1.0 / dx
    invdy = 1.0 / dy

    h_new[:, :] = h
    mx_new[:, :] = mx
    my_new[:, :] = my
    M_new[:, :] = M_base

    for i in prange(H):
        i_up = i - 1 if i > 0 else 0
        i_dn = i + 1 if i < H - 1 else H - 1
        for j in range(W):
            j_lt = j - 1 if j > 0 else 0
            j_rt = j + 1 if j < W - 1 else W - 1

            # X-left interface between (j_lt) and (j)
            hL = h[i, j_lt]; hR = h[i, j]
            mxL = mx[i, j_lt]; mxR = mx[i, j]
            myL = my[i, j_lt]; myR = my[i, j]
            ML = M_base[i, j_lt]; MR = M_base[i, j]
            gpL = gp[i, j_lt]; gpR = gp[i, j]
            uL = u[i, j_lt]; uR = u[i, j]
            pL = 0.5 * gpL * hL * hL
            pR = 0.5 * gpR * hR * hR
            FhL = mxL; FhR = mxR
            FmxL = mxL * uL + pL
            FmxR = mxR * uR + pR
            FmyL = myL * uL
            FmyR = myR * uR
            FML = ML * uL
            FMR = MR * uR
            cL = math.sqrt(max(gpL * hL, 0.0))
            cR = math.sqrt(max(gpR * hR, 0.0))
            amax = abs(uL) + cL
            a2 = abs(uR) + cR
            if a2 > amax:
                amax = a2
            if amax < 1e-6:
                amax = 1e-6
            FxL_h  = 0.5*(FhL + FhR)   - 0.5*amax*(hR - hL)
            FxL_mx = 0.5*(FmxL+ FmxR)  - 0.5*amax*(mxR - mxL)
            FxL_my = 0.5*(FmyL+ FmyR)  - 0.5*amax*(myR - myL)
            FxL_M  = 0.5*(FML + FMR)   - 0.5*amax*(MR - ML)

            # X-right interface between (j) and (j_rt)
            hL2 = h[i, j]; hR2 = h[i, j_rt]
            mxL2 = mx[i, j]; mxR2 = mx[i, j_rt]
            myL2 = my[i, j]; myR2 = my[i, j_rt]
            ML2 = M_base[i, j]; MR2 = M_base[i, j_rt]
            gpL2 = gp[i, j]; gpR2 = gp[i, j_rt]
            uL2 = u[i, j]; uR2 = u[i, j_rt]
            pL2 = 0.5 * gpL2 * hL2 * hL2
            pR2 = 0.5 * gpR2 * hR2 * hR2
            FhL2 = mxL2; FhR2 = mxR2
            FmxL2 = mxL2 * uL2 + pL2
            FmxR2 = mxR2 * uR2 + pR2
            FmyL2 = myL2 * uL2
            FmyR2 = myR2 * uR2
            FML2 = ML2 * uL2
            FMR2 = MR2 * uR2
            cL2 = math.sqrt(max(gpL2 * hL2, 0.0))
            cR2 = math.sqrt(max(gpR2 * hR2, 0.0))
            amax2 = abs(uL2) + cL2
            a22 = abs(uR2) + cR2
            if a22 > amax2:
                amax2 = a22
            if amax2 < 1e-6:
                amax2 = 1e-6
            FxR_h  = 0.5*(FhL2 + FhR2)   - 0.5*amax2*(hR2 - hL2)
            FxR_mx = 0.5*(FmxL2+ FmxR2)  - 0.5*amax2*(mxR2 - mxL2)
            FxR_my = 0.5*(FmyL2+ FmyR2)  - 0.5*amax2*(myR2 - myL2)
            FxR_M  = 0.5*(FML2 + FMR2)   - 0.5*amax2*(MR2 - ML2)

            # Y-up interface between (i_up) and (i)
            hU = h[i_up, j]; hD0 = h[i, j]
            mxU = mx[i_up, j]; mxD0 = mx[i, j]
            myU = my[i_up, j]; myD0 = my[i, j]
            MU = M_base[i_up, j]; MD0 = M_base[i, j]
            gpU = gp[i_up, j]; gpD0 = gp[i, j]
            vU = v[i_up, j]; vD0 = v[i, j]
            pU = 0.5*gpU*hU*hU
            pD0 = 0.5*gpD0*hD0*hD0
            GhU = myU; GhD0 = myD0
            GmxU = mxU*vU
            GmxD0 = mxD0*vD0
            GmyU = myU*vU + pU
            GmyD0 = myD0*vD0 + pD0
            GMU = MU*vU
            GMD0 = MD0*vD0
            cU = math.sqrt(max(gpU*hU, 0.0))
            cD0 = math.sqrt(max(gpD0*hD0, 0.0))
            amaxY = abs(vU) + cU
            aY2 = abs(vD0) + cD0
            if aY2 > amaxY:
                amaxY = aY2
            if amaxY < 1e-6:
                amaxY = 1e-6
            GyU_h  = 0.5*(GhU + GhD0)     - 0.5*amaxY*(hD0 - hU)
            GyU_mx = 0.5*(GmxU+GmxD0)     - 0.5*amaxY*(mxD0 - mxU)
            GyU_my = 0.5*(GmyU+GmyD0)     - 0.5*amaxY*(myD0 - myU)
            GyU_M  = 0.5*(GMU + GMD0)     - 0.5*amaxY*(MD0 - MU)

            # Y-down interface between (i) and (i_dn)
            hU2 = h[i, j]; hD2 = h[i_dn, j]
            mxU2 = mx[i, j]; mxD2 = mx[i_dn, j]
            myU2 = my[i, j]; myD2 = my[i_dn, j]
            MU2 = M_base[i, j]; MD2 = M_base[i_dn, j]
            gpU2 = gp[i, j]; gpD2 = gp[i_dn, j]
            vU2 = v[i, j]; vD2 = v[i_dn, j]
            pU2 = 0.5*gpU2*hU2*hU2
            pD2 = 0.5*gpD2*hD2*hD2
            GhU2 = myU2; GhD2 = myD2
            GmxU2 = mxU2*vU2
            GmxD2 = mxD2*vD2
            GmyU2 = myU2*vU2 + pU2
            GmyD2 = myD2*vD2 + pD2
            GMU2 = MU2*vU2
            GMD2 = MD2*vD2
            cU2 = math.sqrt(max(gpU2*hU2, 0.0))
            cD2 = math.sqrt(max(gpD2*hD2, 0.0))
            amaxY2 = abs(vU2) + cU2
            aY22 = abs(vD2) + cD2
            if aY22 > amaxY2:
                amaxY2 = aY22
            if amaxY2 < 1e-6:
                amaxY2 = 1e-6
            GyD_h  = 0.5*(GhU2 + GhD2)    - 0.5*amaxY2*(hD2 - hU2)
            GyD_mx = 0.5*(GmxU2+GmxD2)    - 0.5*amaxY2*(mxD2 - mxU2)
            GyD_my = 0.5*(GmyU2+GmyD2)    - 0.5*amaxY2*(myD2 - myU2)
            GyD_M  = 0.5*(GMU2 + GMD2)    - 0.5*amaxY2*(MD2 - MU2)

            h_new[i, j]  = h_new[i, j]  - dt * ((FxR_h - FxL_h) * invdx + (GyD_h - GyU_h) * invdy)
            mx_new[i, j] = mx_new[i, j] - dt * ((FxR_mx- FxL_mx)* invdx + (GyD_mx- GyU_mx)* invdy)
            my_new[i, j] = my_new[i, j] - dt * ((FxR_my- FxL_my)* invdx + (GyD_my- GyU_my)* invdy)
            M_new[i, j]  = M_new[i, j]  - dt * ((FxR_M - FxL_M) * invdx + (GyD_M - GyU_M) * invdy)

    _compute_gp_and_uv(h_new, mx_new, my_new, M_new, rho_air, rho_co2, g, gp, u, v)
    for i in prange(H):
        for j in range(W):
            mx_new[i, j] = mx_new[i, j] - gp[i, j] * h_new[i, j] * dz_dx[i, j] * dt
            my_new[i, j] = my_new[i, j] - gp[i, j] * h_new[i, j] * dz_dy[i, j] * dt

    _apply_friction(mx_new, my_new, h_new, Cf, dt, tmp1, tmp2)
    mx_new[:, :] = tmp1
    my_new[:, :] = tmp2

    _diffuse_laplacian(h_new, K_h, dt, dx, dy, tmp0); h_new[:, :] = tmp0
    _diffuse_laplacian(mx_new, K_m, dt, dx, dy, tmp1); mx_new[:, :] = tmp1
    _diffuse_laplacian(my_new, K_m, dt, dx, dy, tmp2); my_new[:, :] = tmp2
    _diffuse_laplacian(M_new, K_M, dt, dx, dy, tmp3); M_new[:, :] = tmp3

    _enforce_consistency_inplace(h_new, mx_new, my_new, M_new, rho_co2)


# -----------------------------
# Hazard (paper-grounded) + lethality
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


def co2_massfrac(M_base: np.ndarray, h: np.ndarray, rho_co2: float) -> np.ndarray:
    denom = (rho_co2 * h + 1e-12)
    f = np.clip(M_base / denom, 0.0, 1.0).astype("float32")
    return f


def rho_bar_from_frac(f: np.ndarray, rho_air: float, rho_co2: float) -> np.ndarray:
    return (rho_air + f * (rho_co2 - rho_air)).astype("float32")


# -----------------------------
# Export helpers
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
    rgba = cm(norm)

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


def utm_bounds_to_wgs84(bounds: Tuple[float, float, float, float], crs: str) -> Dict[str, Any]:
    transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    left, bottom, right, top = bounds
    tl = transformer.transform(left, top)
    tr = transformer.transform(right, top)
    br = transformer.transform(right, bottom)
    bl = transformer.transform(left, bottom)
    xs = [tl[0], tr[0], br[0], bl[0]]
    ys = [tl[1], tr[1], br[1], bl[1]]
    return {"corners_lonlat": [list(tl), list(tr), list(br), list(bl)], "bbox_lonlat": [min(xs), min(ys), max(xs), max(ys)]}


def write_manifest(scenario_dir: Path, scenario_id: str, static: Dict[str, Any], scenario: Dict[str, Any],
                   grid_ds: rasterio.io.DatasetReader, n_frames: int, dt_report_min: int, layers: List[Dict[str, Any]]) -> None:
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
        "time": {"dt_report_min": int(dt_report_min), "n_frames": int(n_frames), "timestamps_min": [int(i * dt_report_min) for i in range(n_frames)]},
        "parameters": scenario,
        "assets": {"metrics": "metrics_10min.csv"},
        "layers": layers,
    }
    (scenario_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def update_index(index_path: Path, scenario_id: str, scenario: Dict[str, Any]) -> None:
    index_path.parent.mkdir(parents=True, exist_ok=True)
    if index_path.exists():
        idx = json.loads(index_path.read_text(encoding="utf-8"))
    else:
        idx = {"scenarios": []}

    idx["scenarios"] = [s for s in idx.get("scenarios", []) if s.get("scenario_id") != scenario_id]
    idx["scenarios"].append({"scenario_id": scenario_id, "path": f"scenario_{scenario_id}/manifest.json", "parameters": scenario})

    opts: Dict[str, Any] = {}
    for s in idx["scenarios"]:
        for k, v in s.get("parameters", {}).items():
            if isinstance(v, (int, float, str, bool)):
                opts.setdefault(k, set()).add(v)
    idx["options"] = {k: sorted(list(v)) for k, v in opts.items()}
    index_path.write_text(json.dumps(idx, indent=2), encoding="utf-8")


def disk_mask(H: int, W: int, center_r: int, center_c: int, radius_m: float, dx: float, dy: float, mask_limit: Optional[np.ndarray] = None) -> np.ndarray:
    rr = np.arange(H)[:, None]
    cc = np.arange(W)[None, :]
    dr_m = (rr - center_r) * dy
    dc_m = (cc - center_c) * dx
    disk = ((dr_m**2 + dc_m**2) <= radius_m**2).astype(np.uint8)
    if mask_limit is not None:
        disk = (disk & mask_limit.astype(np.uint8)).astype(np.uint8)
    print("Disk mask created.")
    return disk


# -----------------------------
# Scenario runner
# -----------------------------
def run_single_scenario(static: Dict[str, Any], scenario_in: Dict[str, Any], out_root: Path, save_geotiff: bool, save_frames: bool) -> Path:
    paths = static["paths"]
    dem_ds, dem = load_raster(Path(paths["dem"]))
    pop_ds, pop = load_raster(Path(paths["pop"]))
    mask_lake_ds, mask_lake = load_raster(Path(paths["mask_lake"]))
    grid_ds, _ = load_raster(Path(paths["grid_template"]))
    print("Rasters loaded.")

    source_idx = json.loads(Path(paths["source_index"]).read_text(encoding="utf-8"))
    src_r, src_c = int(source_idx["row"]), int(source_idx["col"])

    H, W = dem.shape
    dx = float(static["dx_m"]); dy = float(static["dy_m"])
    cell_area = dx * dy

    dz_dy, dz_dx = np.gradient(dem.astype("float32"), dy, dx)
    dz_dx = dz_dx.astype("float32"); dz_dy = dz_dy.astype("float32")

    scenario = dict(scenario_in)
    scenario.setdefault("duration_h", 24.0)
    scenario.setdefault("dt_internal_s", 0.3)
    scenario.setdefault("dt_report_min", 5)
    scenario.setdefault("save_images_every_min", 5)
    scenario.setdefault("exposure_substep_min", 3)

    scenario.setdefault("co2_inventory_upper_m3_stp", 3.0e11)
    scenario.setdefault("release_fraction", 0.1)
    scenario.setdefault("eruption_radius_m", 3000.0)
    scenario.setdefault("constrain_disk_to_lake", True)

    scenario.setdefault("wind_speed_mps", 0.0)
    scenario.setdefault("wind_dir_deg_from", 0.0)

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

    print("Running scenario ")

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

    (scenario_dir / "scenario.json").write_text(json.dumps(scenario, indent=2), encoding="utf-8")

    duration_s = int(round(float(scenario["duration_h"]) * 3600))
    dt = float(scenario["dt_internal_s"])
    dt_report_s = int(round(float(scenario["dt_report_min"]) * 60))
    sub_s = int(round(float(scenario["exposure_substep_min"]) * 60))

    n_steps = int(np.ceil(duration_s / dt))
    report_every = max(1, int(round(dt_report_s / dt)))
    sub_every = max(1, int(round(sub_s / dt)))

    u_wind, v_wind = wind_to_uv(float(scenario["wind_speed_mps"]), float(scenario["wind_dir_deg_from"]))
    wind_mag = float(np.hypot(u_wind, v_wind))

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
    _enforce_consistency_inplace(h, mx, my, M_base, rho_co2)

    # Scratch arrays for Numba
    gp = np.zeros((H, W), dtype="float32")
    u = np.zeros((H, W), dtype="float32")
    v = np.zeros((H, W), dtype="float32")
    h_new = np.zeros((H, W), dtype="float32")
    mx_new = np.zeros((H, W), dtype="float32")
    my_new = np.zeros((H, W), dtype="float32")
    M_new = np.zeros((H, W), dtype="float32")
    tmp0 = np.zeros((H, W), dtype="float32")
    tmp1 = np.zeros((H, W), dtype="float32")
    tmp2 = np.zeros((H, W), dtype="float32")
    tmp3 = np.zeros((H, W), dtype="float32")

    THRESH_DANGER = 50_000.0
    THRESH_LETHAL = 100_000.0
    pop_safe = np.where(np.isfinite(pop.astype("float32")), pop.astype("float32"), 0.0)

    layer_defs = [
        {"name": "ppm_z", "type": "png", "pathTemplate": "frames/ppm_z/t{frame:04d}.png", "opacityDefault": 0.7},
        {"name": "lethal_mask", "type": "png", "pathTemplate": "frames/lethal_mask/t{frame:04d}.png", "opacityDefault": 0.5},
        {"name": "pdeath", "type": "png", "pathTemplate": "frames/pdeath/t{frame:04d}.png", "opacityDefault": 0.6},
        {"name": "base_h", "type": "png", "pathTemplate": "frames/base_h/t{frame:04d}.png", "opacityDefault": 0.5},
    ]
    if save_frames:
        for ld in layer_defs:
            ensure_dir(scenario_dir / "frames" / ld["name"])

    cum_entr_kg = 0.0
    cum_strip_kg = 0.0
    cum_base_loss_kg = 0.0
    ppm: Optional[np.ndarray] = None
    frame_idx = 0
    metrics_rows: List[Dict[str, Any]] = []

    # Warm-up JIT (compilation)
    _twodee_step_numba(h, mx, my, M_base, dz_dx, dz_dy,
                       dt, dx, dy, rho_air, rho_co2, g,
                       float(scenario["Cf_friction"]), float(scenario["K_h_m2ps"]), float(scenario["K_m_m2ps"]), float(scenario["K_M_m2ps"]),
                       gp, u, v, h_new, mx_new, my_new, M_new, tmp0, tmp1, tmp2, tmp3)
    # swap
    h, mx, my, M_base = h_new.copy(), mx_new.copy(), my_new.copy(), M_new.copy()

    for step in tqdm (range(1, n_steps + 1)):
        _twodee_step_numba(h, mx, my, M_base, dz_dx, dz_dy,
                           dt, dx, dy, rho_air, rho_co2, g,
                           float(scenario["Cf_friction"]), float(scenario["K_h_m2ps"]), float(scenario["K_m_m2ps"]), float(scenario["K_M_m2ps"]),
                           gp, u, v, h_new, mx_new, my_new, M_new, tmp0, tmp1, tmp2, tmp3)
        h, mx, my, M_base = h_new, mx_new, my_new, M_new

        entrain_eff = float(scenario["entrain_rate_1ps"]) * (1.0 + float(scenario["entrain_wind_gain"]) * wind_mag)
        dM = (entrain_eff * dt) * M_base
        dM = np.minimum(dM, M_base)
        M_base = (M_base - dM).astype("float32")
        M_top = (M_top + dM).astype("float32")
        cum_entr_kg += float(np.sum(dM) * cell_area)
        _enforce_consistency_inplace(h, mx, my, M_base, rho_co2)

        base_loss = float(scenario["base_background_loss_1ps"])
        if base_loss > 0:
            dM_loss = (1.0 - np.exp(-base_loss * dt)) * M_base
            M_base = (M_base - dM_loss).astype("float32")
            cum_base_loss_kg += float(np.sum(dM_loss) * cell_area)
            _enforce_consistency_inplace(h, mx, my, M_base, rho_co2)

        _advect_upwind_scalar(M_top, u_wind, v_wind, dt, dx, dy, tmp0)
        M_top = tmp0.copy()
        _diffuse_laplacian(M_top, float(scenario["K_top_m2ps"]), dt, dx, dy, tmp1)
        M_top = tmp1.copy()
        M_top = np.maximum(M_top, 0).astype("float32")

        strip_rate = float(scenario["strip_coeff_1ps_per_mps"]) * wind_mag + float(scenario["top_background_loss_1ps"])
        if strip_rate > 0:
            before = float(np.sum(M_top) * cell_area)
            M_top = (M_top * np.exp(-strip_rate * dt)).astype("float32")
            after = float(np.sum(M_top) * cell_area)
            cum_strip_kg += max(0.0, before - after)

        if step % sub_every == 0:
            f_base = co2_massfrac(M_base, h, rho_co2)
            rho_bar_base = rho_bar_from_frac(f_base, rho_air, rho_co2)
            h_top_equiv = (M_top / (rho_co2 + 1e-12)).astype("float32")
            h_eff = (h + h_top_equiv).astype("float32")

            rho_z = rho_at_height_z(h_eff, rho_bar_base, rho_air, float(scenario["S1"]), float(scenario["z_human_m"]))
            ppm = ppm_from_rho(rho_z, rho_air, rho_co2, float(scenario["c_background_ppm"]))

            exposure_min += (ppm >= float(scenario["exposure_threshold_ppm"])).astype("float32") * float(scenario["exposure_substep_min"])
            Pdeath = prob_death_from_ppm(ppm, exposure_min, scenario)

        if step % report_every == 0:
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

            # velocity diagnostic (CPU)
            u0 = np.zeros_like(h); v0 = np.zeros_like(h)
            wet = h > H_FLOOR
            u0[wet] = mx[wet] / h[wet]; v0[wet] = my[wet] / h[wet]
            vmax_flow = float(np.nanmax(np.sqrt(u0*u0 + v0*v0)))

            metrics_rows.append({
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
            })

            if save_frames:
                array_to_rgba_png(ppm, scenario_dir / "frames" / "ppm_z" / f"t{frame_idx:04d}.png",
                                  vmin=0.0, vmax=200_000.0, cmap="magma", alpha_threshold=1_000.0)
                array_to_rgba_png(lethal.astype("float32"), scenario_dir / "frames" / "lethal_mask" / f"t{frame_idx:04d}.png",
                                  vmin=0.0, vmax=1.0, cmap="Reds", alpha_threshold=0.5)
                array_to_rgba_png(Pdeath, scenario_dir / "frames" / "pdeath" / f"t{frame_idx:04d}.png",
                                  vmin=0.0, vmax=1.0, cmap="inferno", alpha_threshold=0.01)
                array_to_rgba_png(h, scenario_dir / "frames" / "base_h" / f"t{frame_idx:04d}.png",
                                  vmin=0.0, vmax=max(float(np.nanmax(h)), 1e-6), cmap="viridis", alpha_threshold=H_FLOOR)
                #print("Saved png.")

            if save_geotiff:
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"h_t{frame_idx:04d}.tif", h, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"M_base_t{frame_idx:04d}.tif", M_base, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"M_top_t{frame_idx:04d}.tif", M_top, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"ppm_t{frame_idx:04d}.tif", ppm, nodata=0.0)
                write_geotiff_like(grid_ds, scenario_dir / "raw" / f"Pdeath_t{frame_idx:04d}.tif", Pdeath, nodata=0.0)
                #print("Saved Geotiff")

            frame_idx += 1

    pd.DataFrame(metrics_rows).to_csv(scenario_dir / "metrics_10min.csv", index=False)

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

    close_datasets([dem_ds, pop_ds, mask_lake_ds, grid_ds])
    return scenario_dir


# -----------------------------
# CLI
# -----------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Run Lake Kivu scenarios and export bundles (Numba accelerated).")
    ap.add_argument("--scenario", type=str, default=None, help="Path to scenario JSON file (single dict).")
    ap.add_argument("--out", type=str, default="scenarios", help="Output root directory (default: scenarios).")
    ap.add_argument("--save-geotiff", action="store_true", help="Write GeoTIFFs at report times (storage-heavy).")
    ap.add_argument("--no-frames", action="store_true", help="Skip PNG frame export.")
    return ap.parse_args()


def main() -> None:
    project_root = Path(".").resolve()
    proc_dir = project_root / "data" / "processed"
    static_path = proc_dir / "static_layers.json"
    if not static_path.exists():
        print(static_path)
        raise FileNotFoundError(f"Missing {static_path}. Run preprocessing notebooks first.")
    static = load_json(static_path)
    print(f"Loaded static layers from {static_path}")

    args = parse_args()
    scenario_in: Dict[str, Any] = {}
    if args.scenario:
        scenario_in = load_json(Path(args.scenario))

    out_root = Path(args.out).resolve()
    ensure_dir(out_root)
    print(f"Output root directory: {out_root}")

    scenario_dir = run_single_scenario(static, scenario_in, out_root, bool(args.save_geotiff), not bool(args.no_frames))

    print(json.dumps({
        "status": "ok",
        "scenario_dir": str(scenario_dir),
        "manifest": str(scenario_dir / "manifest.json"),
        "metrics": str(scenario_dir / "metrics_10min.csv"),
    }, indent=2))


if __name__ == "__main__":
    main()
