"""Typed configuration and result containers."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Literal, Dict

ReleaseProfile = Literal["pulse", "top_hat", "exp_decay", "multi_pulse"]
SourceType = Literal["point", "polygon"]

@dataclass
class GridConfig:
    cell_size_m: float
    crs: str

@dataclass
class TimeConfig:
    duration_h: float
    dt_internal_s: float
    dt_report_h: float

@dataclass
class DispersionConfig:
    eddy_diffusivity_m2ps: float
    loss_rate_1ps: float
    slope_strength: float
    slope_vmax_mps: float

@dataclass
class DataConfig:
    dem: str
    pop: str
    lake_outline: str
    settlements: Optional[str] = None

@dataclass
class ModelConfig:
    grid: GridConfig
    time: TimeConfig
    dispersion: DispersionConfig
    data: DataConfig

@dataclass
class SourceConfig:
    source_type: SourceType
    release_profile: ReleaseProfile
    released_mass_kg: float
    release_duration_s: float
    initial_spread_radius_m: float = 0.0
    source_lon: Optional[float] = None
    source_lat: Optional[float] = None

@dataclass
class MetConfig:
    wind_speed_mps: float
    wind_dir_deg: float

@dataclass
class Thresholds:
    danger_ppm: float
    lethal_ppm: float

@dataclass
class RunMetrics:
    values: Dict[str, float]

@dataclass
class RunResult:
    metrics: RunMetrics
    map_paths: Optional[Dict[str, str]] = None
