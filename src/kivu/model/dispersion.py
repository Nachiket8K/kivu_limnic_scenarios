"""Dispersion numerical core scaffold."""
from __future__ import annotations
import numpy as np

def compute_slope_velocity(dem: np.ndarray, dx: float, dy: float, slope_strength: float, vmax: float) -> tuple[np.ndarray, np.ndarray]:
    dz_dy, dz_dx = np.gradient(dem, dy, dx)
    u = -slope_strength * dz_dx
    v = -slope_strength * dz_dy
    mag = np.sqrt(u*u + v*v) + 1e-12
    scale = np.minimum(1.0, vmax / mag)
    return u*scale, v*scale

def step_advection_diffusion(C: np.ndarray, u: np.ndarray, v: np.ndarray, K: float, dt: float, dx: float, dy: float) -> np.ndarray:
    # TODO: implement stable upwind advection + diffusion
    return C.copy()

def apply_loss(C: np.ndarray, loss_rate_1ps: float, dt: float) -> np.ndarray:
    return C * np.exp(-loss_rate_1ps * dt)
