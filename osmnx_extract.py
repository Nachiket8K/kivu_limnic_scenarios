#!/usr/bin/env python3
"""
osmnx_extract.py — Extract OSM context layers for the AOI using OSMnx.

Outputs (in ./osm/):
- aoi.geojson
- roads.geojson
- water.geojson
- buildings.geojson (optional; can be large)
"""

from pathlib import Path
import json

import geopandas as gpd
import shapely.geometry as geom
import osmnx as ox
from pyproj import Transformer


PROJECT_ROOT = Path(".").resolve()
PROC = PROJECT_ROOT / "data" / "processed"
STATIC = PROC / "static_layers.json"
OUT = PROJECT_ROOT / "osm"
OUT.mkdir(parents=True, exist_ok=True)

static = json.loads(STATIC.read_text(encoding="utf-8"))
paths = static["paths"]

# Load grid_template bounds (projected CRS)
import rasterio
grid_path = Path(paths["grid_template"])
with rasterio.open(grid_path) as ds:
    crs_proj = ds.crs.to_string()
    b = ds.bounds  # left,bottom,right,top in projected CRS

# Convert AOI bounds to WGS84
tf = Transformer.from_crs(crs_proj, "EPSG:4326", always_xy=True)
minlon, minlat = tf.transform(b.left, b.bottom)
maxlon, maxlat = tf.transform(b.right, b.top)

print("[AOI] WGS84 bounds:", (minlon, minlat, maxlon, maxlat))

# AOI polygon in WGS84
aoi_poly = geom.box(minlon, minlat, maxlon, maxlat)
aoi_gdf = gpd.GeoDataFrame({"name": ["aoi"]}, geometry=[aoi_poly], crs="EPSG:4326")
aoi_gdf.to_file(OUT / "aoi.geojson", driver="GeoJSON")
print("[WRITE]", OUT / "aoi.geojson")

# OSMnx config
ox.settings.use_cache = True
ox.settings.log_console = True

# ---- Roads (as graph -> edges) ----
print("[OSM] Download roads graph...")
G = ox.graph_from_polygon(aoi_poly, network_type="drive")  # change to 'all' if you want footpaths etc.
gdf_nodes, gdf_edges = ox.graph_to_gdfs(G)

# Keep useful attributes; drop super noisy columns
keep_cols = [c for c in ["highway", "name", "oneway", "lanes", "maxspeed"] if c in gdf_edges.columns]
roads = gdf_edges[keep_cols + ["geometry"]].reset_index(drop=True)
roads = roads.to_crs("EPSG:4326")
roads.to_file(OUT / "roads.geojson", driver="GeoJSON")
print("[WRITE]", OUT / "roads.geojson", "| features:", len(roads))

# ---- Water bodies ----
print("[OSM] Download water features...")
water_tags = {
    "natural": ["water"],
    "water": True,
    "waterway": True
}
water = ox.features_from_polygon(aoi_poly, tags=water_tags)
water = water[water.geometry.notnull()].copy()
water = water.to_crs("EPSG:4326")

# Keep only polygons/lines (exclude points)
water = water[water.geometry.type.isin(["Polygon", "MultiPolygon", "LineString", "MultiLineString"])].copy()
water.to_file(OUT / "water.geojson", driver="GeoJSON")
print("[WRITE]", OUT / "water.geojson", "| features:", len(water))

# ---- Buildings (optional, can be huge) ----
print("[OSM] Download buildings (optional)...")
bldg_tags = {"building": True}
bldg = ox.features_from_polygon(aoi_poly, tags=bldg_tags)
bldg = bldg[bldg.geometry.notnull()].copy()
bldg = bldg.to_crs("EPSG:4326")
bldg = bldg[bldg.geometry.type.isin(["Polygon", "MultiPolygon"])].copy()

# You may want to simplify buildings to reduce size:
bldg["geometry"] = bldg["geometry"].simplify(0.0001, preserve_topology=True)

bldg.to_file(OUT / "buildings.geojson", driver="GeoJSON")
print("[WRITE]", OUT / "buildings.geojson", "| features:", len(bldg))

print("[DONE] OSM extraction complete. Next step: vector tiles (PMTiles).")