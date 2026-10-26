"""Synthetic multi-plate screening data with a KNOWN best analysis configuration.

Structure mirrors a real high-content run: several plates, control wells on
every plate, ~150 cells per well. Built-in ground truth:

  * six candidate readouts (zone_1..zone_6); control separation peaks at zone_3
  * a plate effect shared by both controls (day-to-day drift), which paired
    SSMD should remove and unpaired SSMD should be hurt by
  * cells with extreme expression (red_norm outside ~0.5-2.5) are mis-measured,
    so a sensible expression gate should improve quality

    python scripts/make_synthetic_data.py --out data/synthetic_cells.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

SEPARATION = {"zone_1": 0.05, "zone_2": 0.20, "zone_3": 0.60,
              "zone_4": 0.30, "zone_5": 0.15, "zone_6": 0.05}
BEST_READOUT = "zone_3"


def make(n_plates: int = 8, wells_per_control: int = 3, cells_per_well: int = 300,
         plate_sd: float = 0.3, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for p in range(1, n_plates + 1):
        plate_shift = rng.normal(0, plate_sd)                 # shared by both controls
        for cond, is_pos, col in (("WT_DMSO", True, 2), ("F508del_DMSO", False, 3)):
            ctrl_shift = rng.normal(0, 0.05)                  # small plate x control interaction
            for w in range(wells_per_control):
                well = f"{'BCDEFG'[w]}{col:02d}"
                n = cells_per_well
                red_norm = rng.lognormal(0, 0.6, n)
                in_range = (red_norm > 0.5) & (red_norm < 2.5)
                cell = {"plate": f"P{p}", "well": well, "condition": cond,
                        "red_norm": red_norm.round(4)}
                for zone, sep in SEPARATION.items():
                    signal = plate_shift + ctrl_shift + (sep if is_pos else 0.0)
                    good = signal + rng.normal(0, 0.5, n)
                    bad = plate_shift + rng.normal(0, 5.0, n)  # mis-measured: no signal, much more noise
                    cell[zone] = np.where(in_range, good, bad).round(4)
                rows.append(pd.DataFrame(cell))
    return pd.concat(rows, ignore_index=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("data/synthetic_cells.csv"))
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    make(seed=a.seed).to_csv(a.out, index=False)
    print(f"Wrote {a.out}")
