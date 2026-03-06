"""Metric helpers scaffold."""
from __future__ import annotations
import numpy as np

def area_km2(mask: np.ndarray, cell_area_m2: float) -> float:
    return float(mask.sum() * cell_area_m2 / 1e6)
