"""Batch runner scaffold (writes placeholder metrics.csv)."""
from __future__ import annotations
import argparse
import pandas as pd
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=str, default="runs/outputs/metrics.csv")
    args = p.parse_args()

    df = pd.DataFrame([{"scenario_id": i} for i in range(args.n)])
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Wrote placeholder metrics to {out}")

if __name__ == "__main__":
    main()
