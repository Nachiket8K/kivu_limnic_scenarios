"""Experiment design scaffold."""
from __future__ import annotations
import numpy as np

def generate_scenarios_placeholder(n: int, seed: int = 42):
    rng = np.random.default_rng(seed)
    return [{"scenario_id": i, "seed": int(rng.integers(0, 2**31-1))} for i in range(n)]
