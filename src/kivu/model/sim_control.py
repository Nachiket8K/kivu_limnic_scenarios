#!/usr/bin/env python3
r"""
Batch controller for Lake Kivu limnic eruption simulations.

This preserves `sim_run_gpu.py` and reuses its existing scenario runner to
execute a parameter sweep across release fraction and wind settings.

Example (Windows cmd):
  python sim_control.py --scenario scenariospec.json --out ..\..\..\docs\scenarios ^
    --release-fraction 0.1 0.25 0.5 ^
    --wind-speed-mps 0 2 5 ^
    --wind-dir-deg-from 0 90 180 270
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

import pandas as pd

from sim_run_gpu import ensure_dir, load_json, run_single_scenario


def _float_values(values: Sequence[float] | None, default: float) -> List[float]:
    if not values:
        return [float(default)]
    return [float(v) for v in values]


def build_scenarios(
    base_scenario: Dict[str, Any],
    release_fractions: Iterable[float],
    wind_speeds: Iterable[float],
    wind_dirs: Iterable[float],
) -> List[Dict[str, Any]]:
    scenarios: List[Dict[str, Any]] = []
    for release_fraction, wind_speed_mps, wind_dir_deg_from in itertools.product(
        release_fractions,
        wind_speeds,
        wind_dirs,
    ):
        scenario = dict(base_scenario)
        scenario["release_fraction"] = float(release_fraction)
        scenario["wind_speed_mps"] = float(wind_speed_mps)
        scenario["wind_dir_deg_from"] = float(wind_dir_deg_from)
        scenarios.append(scenario)
    return scenarios


def scenario_label(scenario: Dict[str, Any]) -> str:
    return (
        f"rf={float(scenario['release_fraction']):g}"
        f" | ws={float(scenario['wind_speed_mps']):g} m/s"
        f" | wd={float(scenario['wind_dir_deg_from']):g} deg"
    )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Batch-run Lake Kivu GPU scenarios across parameter combinations."
    )
    ap.add_argument(
        "--scenario",
        type=str,
        default=None,
        help="Path to a base scenario JSON file with shared settings.",
    )
    ap.add_argument(
        "--out",
        type=str,
        default="scenarios_batch",
        help="Output root directory for the batch run.",
    )
    ap.add_argument(
        "--release-fraction",
        type=float,
        nargs="+",
        default=None,
        help="One or more release_fraction values.",
    )
    ap.add_argument(
        "--wind-speed-mps",
        type=float,
        nargs="+",
        default=None,
        help="One or more wind_speed_mps values.",
    )
    ap.add_argument(
        "--wind-dir-deg-from",
        type=float,
        nargs="+",
        default=None,
        help="One or more wind_dir_deg_from values.",
    )
    ap.add_argument(
        "--save-geotiff",
        action="store_true",
        help="Write GeoTIFFs for each scenario at report times.",
    )
    ap.add_argument(
        "--no-frames",
        action="store_true",
        help="Skip PNG frame export for each scenario.",
    )
    ap.add_argument(
        "--backend",
        choices=["cpu", "gpu"],
        default="cpu",
        help="Execution backend for supported simulation phases.",
    )
    return ap.parse_args()


def main() -> None:
    project_root = Path(".").resolve()
    proc_dir = project_root / "data" / "processed"
    static_path = proc_dir / "static_layers.json"
    if not static_path.exists():
        raise FileNotFoundError(f"Missing {static_path}. Run preprocessing notebooks first.")

    static = load_json(static_path)
    args = parse_args()

    base_scenario: Dict[str, Any] = {}
    if args.scenario:
        base_scenario = load_json(Path(args.scenario))

    release_fractions = _float_values(args.release_fraction, base_scenario.get("release_fraction", 0.25))
    wind_speeds = _float_values(args.wind_speed_mps, base_scenario.get("wind_speed_mps", 0.0))
    wind_dirs = _float_values(args.wind_dir_deg_from, base_scenario.get("wind_dir_deg_from", 0.0))

    scenarios = build_scenarios(base_scenario, release_fractions, wind_speeds, wind_dirs)
    out_root = Path(args.out).resolve()
    ensure_dir(out_root)

    batch_summary: List[Dict[str, Any]] = []
    total = len(scenarios)

    print(f"Loaded static layers from {static_path}")
    print(f"Output root directory: {out_root}")
    print(f"Running {total} scenario(s) in batch.")

    for idx, scenario in enumerate(scenarios, start=1):
        label = scenario_label(scenario)
        print(f"[{idx}/{total}] {label}")
        scenario_dir = run_single_scenario(
            static=static,
            scenario_in=scenario,
            out_root=out_root,
            save_geotiff=bool(args.save_geotiff),
            save_frames=not bool(args.no_frames),
            backend=str(args.backend),
        )
        batch_summary.append(
            {
                "batch_index": idx,
                "label": label,
                "scenario_dir": str(scenario_dir),
                "release_fraction": float(scenario["release_fraction"]),
                "wind_speed_mps": float(scenario["wind_speed_mps"]),
                "wind_dir_deg_from": float(scenario["wind_dir_deg_from"]),
            }
        )

    summary_csv = out_root / "batch_runs.csv"
    summary_json = out_root / "batch_runs.json"
    pd.DataFrame(batch_summary).to_csv(summary_csv, index=False)
    summary_json.write_text(json.dumps(batch_summary, indent=2), encoding="utf-8")

    print(
        json.dumps(
            {
                "status": "ok",
                "runs": total,
                "output_root": str(out_root),
                "batch_runs_csv": str(summary_csv),
                "batch_runs_json": str(summary_json),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()