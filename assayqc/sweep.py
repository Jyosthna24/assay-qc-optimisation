"""Grid search over analysis parameters, scored by assay quality.

Image-analysis pipelines have free parameters (which readout, which cells to
keep) that change the answer. Instead of choosing them by eye, score every
combination by how well it separates a positive from a negative control and
pick the best one, with the full grid kept as evidence.

    readouts : candidate columns, e.g. membrane-zone definitions or AAC time points
    gates    : (low, high) windows on a gating column, e.g. normalised expression
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from .metrics import classify_ssmd, classify_z_prime, signal_window, ssmd, ssmd_paired, z_prime


def control_values(
    df: pd.DataFrame, value: str, control: str,
    condition_col: str = "condition", plate_col: str = "plate", well_col: str = "well",
    level: str = "plate",
) -> pd.Series:
    """Values for one control at the chosen unit of analysis.

    plate - mean per plate (index = plate), so controls can be paired
    well  - mean per well
    cell  - every cell
    """
    d = df.loc[df[condition_col] == control, [plate_col, well_col, value]].dropna(subset=[value])
    if level == "plate":
        return d.groupby(plate_col)[value].mean()
    if level == "well":
        return d.groupby([plate_col, well_col])[value].mean()
    if level == "cell":
        return d[value].reset_index(drop=True)
    raise ValueError("level must be plate, well or cell")


def score(
    df: pd.DataFrame, value: str, pos: str, neg: str, level: str = "plate", **cols
) -> dict:
    p = control_values(df, value, pos, level=level, **cols)
    n = control_values(df, value, neg, level=level, **cols)
    out = {"n_pos": len(p), "n_neg": len(n)}
    for robust, tag in ((True, "robust"), (False, "classical")):
        out[f"zprime_{tag}"] = z_prime(p, n, robust)
        out[f"ssmd_{tag}"] = ssmd(p, n, robust)
        if level == "plate":
            common = p.index.intersection(n.index)
            out[f"ssmd_paired_{tag}"] = ssmd_paired(p[common], n[common], robust)
    out["signal_window"] = signal_window(p, n)
    return out


def sweep(
    df: pd.DataFrame,
    readouts: Sequence[str],
    pos: str,
    neg: str,
    gate_col: str | None = None,
    gates: Sequence[tuple[float, float]] = ((-np.inf, np.inf),),
    level: str = "plate",
    condition_col: str = "condition",
    plate_col: str = "plate",
    well_col: str = "well",
) -> pd.DataFrame:
    """Long table: one row per (gate, readout) with every quality metric."""
    cols = dict(condition_col=condition_col, plate_col=plate_col, well_col=well_col)
    n_total = len(df)
    rows = []
    for lo, hi in gates:
        sub = df if gate_col is None else df[df[gate_col].between(lo, hi)]
        for r in readouts:
            row = {"gate_low": lo, "gate_high": hi, "readout": r,
                   "fraction_cells_kept": len(sub) / n_total if n_total else np.nan}
            row |= score(sub, r, pos, neg, level=level, **cols)
            rows.append(row)
    return pd.DataFrame(rows)


def best(results: pd.DataFrame, metric: str | None = None, min_fraction_kept: float = 0.5,
         min_zprime: float | None = None) -> pd.Series:
    """Best configuration by `metric`, subject to keeping enough cells (and optionally a
    minimum robust Z').

    The retention constraint stops the search from 'winning' by discarding most of the data.
    """
    if metric is None:
        metric = "ssmd_paired_robust" if "ssmd_paired_robust" in results else "ssmd_robust"
    ok = results["fraction_cells_kept"] >= min_fraction_kept
    if min_zprime is not None:
        ok &= results["zprime_robust"] >= min_zprime
    cand = results[ok & results[metric].notna()]
    if cand.empty:
        raise ValueError("no configuration satisfies the constraints")
    row = cand.loc[cand[metric].idxmax()].copy()
    row["ssmd_quality"] = classify_ssmd(row[metric])
    row["zprime_quality"] = classify_z_prime(row["zprime_robust"])
    return row


def to_grid(results: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Pivot to a heatmap-ready grid: rows = gate, columns = readout (input order kept)."""
    r = results.assign(gate=results["gate_low"].map("{:g}".format) + "–"
                       + results["gate_high"].map("{:g}".format))
    grid = r.pivot(index="gate", columns="readout", values=metric)
    return grid.loc[pd.unique(r["gate"]), pd.unique(r["readout"])]
