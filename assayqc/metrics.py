"""Assay-quality metrics: Z'-factor and SSMD, classical and robust, unpaired and paired.

P = positive control, N = negative control.

    MAD            = 1.4826 * median(|x - median(x)|)       (consistent with SD for normal data)
    Z'             = 1 - 3 (s_P + s_N) / |m_P - m_N|        Zhang et al. 1999
    SSMD           = (m_P - m_N) / sqrt(s_P^2 + s_N^2)      Zhang 2007
    paired SSMD    = mean(D) / SD(D),  D_i = P_i - N_i      Zhang 2011, per-plate differences

Classical versions use mean / SD; robust versions use median / MAD, so a few
outlying plates or wells cannot dominate the score.

The paired estimator is the right one when both controls sit on every plate:
plate-to-plate drift moves P and N together, and the unpaired formula wrongly
counts that shared drift as noise.
"""
from __future__ import annotations

import numpy as np

MAD_SCALE = 1.4826


def _clean(x) -> np.ndarray:
    x = np.asarray(x, float)
    return x[np.isfinite(x)]


def mad(x) -> float:
    x = _clean(x)
    return MAD_SCALE * float(np.median(np.abs(x - np.median(x)))) if x.size else np.nan


def _loc_scale(x, robust: bool) -> tuple[float, float]:
    x = _clean(x)
    if x.size < 2:
        return np.nan, np.nan
    return (float(np.median(x)), mad(x)) if robust else (float(np.mean(x)), float(np.std(x, ddof=1)))


def z_prime(pos, neg, robust: bool = True) -> float:
    mp, sp = _loc_scale(pos, robust)
    mn, sn = _loc_scale(neg, robust)
    spread = abs(mp - mn)
    return float(1 - 3 * (sp + sn) / spread) if spread > 0 else np.nan


def ssmd(pos, neg, robust: bool = True) -> float:
    """Signed: negative when the 'positive' control falls below the negative."""
    mp, sp = _loc_scale(pos, robust)
    mn, sn = _loc_scale(neg, robust)
    denom = np.sqrt(sp ** 2 + sn ** 2)
    return float((mp - mn) / denom) if denom > 0 else np.nan


def ssmd_paired(pos, neg, robust: bool = True) -> float:
    """SSMD from per-plate differences; pos[i] and neg[i] must come from the same plate."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if pos.shape != neg.shape:
        raise ValueError("paired SSMD needs one positive and one negative value per plate")
    d = _clean(pos - neg)
    m, s = _loc_scale(d, robust)
    return float(m / s) if s and s > 0 else np.nan


def signal_window(pos, neg, robust: bool = True) -> float:
    """Fold separation of the control means (signal-to-background)."""
    mp, _ = _loc_scale(pos, robust)
    mn, _ = _loc_scale(neg, robust)
    return float(mp / mn) if mn else np.nan


def classify_z_prime(z: float) -> str:
    if not np.isfinite(z):
        return "undefined"
    if z >= 0.5:
        return "excellent"
    if z > 0:
        return "acceptable"
    return "unusable"


def classify_ssmd(beta: float) -> str:
    """Quality bands for moderate-effect positive controls (Zhang 2007)."""
    if not np.isfinite(beta):
        return "undefined"
    b = abs(beta)
    if b >= 3:
        return "excellent"
    if b >= 2:
        return "good"
    if b >= 1:
        return "inferior"
    return "poor"
