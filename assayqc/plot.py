"""Annotated heatmap of a quality metric across the parameter grid."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def heatmap(grid: pd.DataFrame, title: str, best: tuple[str, str] | None = None,
            xlabel: str = "readout", ylabel: str = "gate", ax=None):
    arr = grid.to_numpy(float)
    ax = ax or plt.subplots(figsize=(max(6, 0.55 * arr.shape[1]), max(4, 0.45 * arr.shape[0])))[1]
    im = ax.imshow(arr, cmap="RdYlGn", aspect="auto")
    ax.set_xticks(range(arr.shape[1]), grid.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(arr.shape[0]), grid.index, fontsize=7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if np.isfinite(arr[i, j]):
                ax.text(j, i, f"{arr[i, j]:.2f}", ha="center", va="center", fontsize=6)
    if best is not None:
        i, j = list(grid.index).index(best[0]), list(grid.columns).index(best[1])
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=2, ec="black"))
    ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return ax
