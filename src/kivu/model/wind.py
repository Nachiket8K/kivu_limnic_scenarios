"""Wind field model scaffold."""
from __future__ import annotations
import numpy as np

def wind_uv(speed_mps: float, dir_deg: float) -> tuple[float, float]:
    rad = np.deg2rad(dir_deg)
    u = speed_mps * np.cos(rad)
    v = speed_mps * np.sin(rad)
    return float(u), float(v)
