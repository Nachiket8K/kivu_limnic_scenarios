"""Vector I/O helpers (GeoJSON/Shapefile) using geopandas."""
from __future__ import annotations
from pathlib import Path
import geopandas as gpd

def read_vector(path: str | Path) -> gpd.GeoDataFrame:
    return gpd.read_file(path)
