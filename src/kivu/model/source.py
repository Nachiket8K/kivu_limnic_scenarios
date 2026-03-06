"""Source term model scaffold."""
from __future__ import annotations
import numpy as np

def release_rate_kgps(t_s: float, total_mass_kg: float, duration_s: float, profile: str) -> float:
    if profile == "pulse":
        return total_mass_kg if t_s <= 0 else 0.0
    if profile == "top_hat":
        return total_mass_kg / duration_s if 0 <= t_s <= duration_s else 0.0
    if profile == "exp_decay":
        tau = max(1.0, duration_s / 3.0)
        return (total_mass_kg / tau) * np.exp(-t_s / tau) if t_s >= 0 else 0.0
    raise ValueError(f"Unknown profile: {profile}")

def make_source_mask_placeholder(shape) -> np.ndarray:
    return np.zeros(shape, dtype=np.float32)
