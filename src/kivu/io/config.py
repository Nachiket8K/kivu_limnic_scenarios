"""Configuration utilities (YAML parsing)."""
from __future__ import annotations
from pathlib import Path
import yaml

def load_yaml(path: str | Path) -> dict:
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)
