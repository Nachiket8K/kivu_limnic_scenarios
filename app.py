# app.py
from __future__ import annotations

import json
from pathlib import Path
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

PROJECT_ROOT = Path(__file__).resolve().parent
SCENARIOS_DIR = PROJECT_ROOT / "scenarios"
INDEX_PATH = SCENARIOS_DIR / "index.json"

app = FastAPI(title="Kivu Scenario Server")

@app.get("/api/index")
def get_index():
    if not INDEX_PATH.exists():
        raise HTTPException(status_code=404, detail=f"Missing {INDEX_PATH}")
    return json.loads(INDEX_PATH.read_text(encoding="utf-8"))

@app.get("/api/scenario/{scenario_id}/manifest")
def get_manifest(scenario_id: str):
    manifest = SCENARIOS_DIR / f"scenario_{scenario_id}" / "manifest.json"
    if not manifest.exists():
        raise HTTPException(status_code=404, detail="Manifest not found")
    return json.loads(manifest.read_text(encoding="utf-8"))

# Serve scenario files (frames, metrics, raw, manifest)
app.mount("/scenarios", StaticFiles(directory=str(SCENARIOS_DIR)), name="scenarios")

# Serve the frontend
FRONTEND_DIR = PROJECT_ROOT / "web"
app.mount("/", StaticFiles(directory=str(FRONTEND_DIR), html=True), name="web")