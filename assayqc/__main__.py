"""Command-line interface: python -m assayqc CELLS.csv --pos P --neg N [options]"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import pandas as pd

from .sweep import best, sweep, to_grid


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(prog="assayqc",
                                description="Choose analysis parameters by assay quality (SSMD / Z').")
    p.add_argument("cells", type=Path, help="per-cell table (CSV or Excel)")
    p.add_argument("--pos", required=True, help="positive-control condition")
    p.add_argument("--neg", required=True, help="negative-control condition")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--readouts", nargs="+", help="candidate readout columns")
    g.add_argument("--readout-prefix", help="use every column starting with this prefix")
    p.add_argument("--gate-col", default=None, help="column to gate cells on, e.g. red_norm")
    p.add_argument("--gate-lower", nargs="+", type=float, default=[0.05, 0.25, 0.5, 0.75, 1.0])
    p.add_argument("--gate-upper", nargs="+", type=float, default=[2.5, 3.0, 3.5, 4.0])
    p.add_argument("--level", choices=["plate", "well", "cell"], default="plate")
    p.add_argument("--min-fraction-kept", type=float, default=0.5)
    p.add_argument("--min-zprime", type=float, default=None, help="optional robust Z' floor")
    p.add_argument("--condition-col", default="condition")
    p.add_argument("--plate-col", default="plate")
    p.add_argument("--well-col", default="well")
    p.add_argument("-o", "--outdir", type=Path, default=Path("results"))
    a = p.parse_args(argv)

    df = pd.read_excel(a.cells) if a.cells.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(a.cells)
    readouts = a.readouts or [c for c in df.columns if c.startswith(a.readout_prefix)]
    gates = (list(itertools.product(a.gate_lower, a.gate_upper)) if a.gate_col
             else [(float("-inf"), float("inf"))])

    res = sweep(df, readouts, a.pos, a.neg, gate_col=a.gate_col, gates=gates, level=a.level,
                condition_col=a.condition_col, plate_col=a.plate_col, well_col=a.well_col)
    a.outdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.outdir / "sweep_results.csv", index=False)

    metric = "ssmd_paired_robust" if a.level == "plate" else "ssmd_robust"
    top = best(res, metric, min_fraction_kept=a.min_fraction_kept, min_zprime=a.min_zprime)
    (a.outdir / "best_configuration.json").write_text(top.to_json(indent=2))
    print(f"{len(res)} configurations scored ({len(gates)} gates x {len(readouts)} readouts)")
    print(f"Best: readout={top['readout']}  gate=[{top['gate_low']:g}, {top['gate_high']:g}]  "
          f"{metric}={top[metric]:.2f} ({top['ssmd_quality']})  "
          f"Z'={top['zprime_robust']:.2f} ({top['zprime_quality']})  "
          f"cells kept={top['fraction_cells_kept']:.0%}")

    from .plot import heatmap
    best_cell = (f"{top['gate_low']:g}–{top['gate_high']:g}", top["readout"])
    for m in (metric, "zprime_robust"):
        ax = heatmap(to_grid(res, m), f"{m}: {a.pos} vs {a.neg} ({a.level} level)",
                     best=best_cell, xlabel="readout", ylabel=f"{a.gate_col or 'gate'} window")
        ax.figure.tight_layout()
        ax.figure.savefig(a.outdir / f"heatmap_{m}.png", dpi=200)
    print(f"Results in {a.outdir.resolve()}")


if __name__ == "__main__":
    main()
