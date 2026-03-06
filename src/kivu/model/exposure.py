"""Exposure calculations scaffold."""
from __future__ import annotations
import numpy as np

def mask_exceed(C_ppm: np.ndarray, threshold_ppm: float) -> np.ndarray:
    return (C_ppm >= threshold_ppm)

def population_exposed(mask: np.ndarray, pop: np.ndarray) -> float:
    return float(np.nansum(pop[mask]))
