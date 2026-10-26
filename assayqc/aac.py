"""Area above the curve (AAC) of normalised quench traces, vectorised over cells.

A cell's YFP signal y(t) is normalised so y(t0) = 1; anion influx quenches it.
More channel activity -> faster quench -> larger area above the curve:

    AAC(t) = (t - t0) - integral_{t0}^{t} y dt        (trapezoid rule)

Missing readings are skipped: the trapezoid bridges to the next valid point,
which equals integrating over the valid subset only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def aac_matrix(traces: pd.DataFrame | np.ndarray, times, t0: float = 0.0) -> np.ndarray:
    """(n_cells, n_times) traces -> (n_cells, n_times) AAC up to each time point."""
    y_all = np.asarray(traces, float)
    times = np.asarray(times, float)
    n = y_all.shape[0]
    last_x, last_y = np.full(n, t0), np.ones(n)
    auc = np.zeros(n)
    out = np.empty_like(y_all)
    for j, t in enumerate(times):
        y = y_all[:, j]
        valid = np.isfinite(y)
        auc = np.where(valid, auc + (t - last_x) * (y + last_y) / 2.0, auc)
        last_x = np.where(valid, t, last_x)
        last_y = np.where(valid, y, last_y)
        out[:, j] = (t - t0) - auc
    return out
