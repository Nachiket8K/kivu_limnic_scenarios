"""Raster I/O helpers using rasterio."""
from __future__ import annotations
from pathlib import Path
import rasterio
import numpy as np

def read_raster(path: str | Path):
    path = Path(path)
    with rasterio.open(path) as ds:
        arr = ds.read(1)
        return arr, ds.transform, ds.crs, ds.nodata

def write_raster(path: str | Path, arr: np.ndarray, transform, crs, nodata=None, dtype=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if dtype is None:
        dtype = arr.dtype
    height, width = arr.shape
    with rasterio.open(
        path, "w", driver="GTiff",
        height=height, width=width, count=1,
        dtype=dtype, crs=crs, transform=transform,
        nodata=nodata, compress="deflate"
    ) as ds:
        ds.write(arr, 1)
