import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from assayqc import (aac_matrix, best, classify_ssmd, classify_z_prime, mad, ssmd,  # noqa: E402
                     ssmd_paired, sweep, to_grid, z_prime)
from make_synthetic_data import BEST_READOUT, SEPARATION, make  # noqa: E402

POS, NEG = [9, 10, 11], [0, 1, 2]


def test_classical_metrics_by_hand():
    assert z_prime(POS, NEG, robust=False) == pytest.approx(1 - 6 / 9)
    assert ssmd(POS, NEG, robust=False) == pytest.approx(9 / np.sqrt(2))
    assert ssmd(NEG, POS, robust=False) == pytest.approx(-9 / np.sqrt(2))  # sign kept


def test_robust_metrics_by_hand():
    s = 1.4826
    assert mad(POS) == pytest.approx(s)
    assert z_prime(POS, NEG) == pytest.approx(1 - 3 * 2 * s / 9)
    assert ssmd(POS, NEG) == pytest.approx(9 / (s * np.sqrt(2)))


def test_mad_matches_sd_for_normal_data():
    x = np.random.default_rng(0).normal(0, 2, 200_000)
    assert mad(x) == pytest.approx(2, rel=0.01)


def test_robust_z_prime_resists_an_outlier_plate():
    pos = [10, 10.2, 9.8, 10.1, 9.9, 30]   # one broken plate
    neg = [1, 1.2, 0.8, 1.1, 0.9, 1.0]
    assert z_prime(pos, neg, robust=False) < 0 < z_prime(pos, neg, robust=True)


def test_paired_ssmd_removes_shared_plate_drift():
    drift = np.array([-1.0, 0.5, 2.0, -0.5, 1.0])
    pos = 1.0 + drift + np.array([0.05, -0.05, 0.02, -0.02, 0.0])
    neg = 0.0 + drift
    assert ssmd_paired(pos, neg, robust=False) > 10
    assert ssmd(pos, neg, robust=False) < 1


def test_quality_bands():
    assert classify_z_prime(0.6) == "excellent" and classify_z_prime(-0.1) == "unusable"
    assert classify_ssmd(3.5) == "excellent" and classify_ssmd(1.5) == "inferior"


def test_aac_trapezoid_and_missing_points():
    times = [1, 2, 3]
    assert aac_matrix(np.ones((1, 3)), times).ravel() == pytest.approx([0, 0, 0])
    assert aac_matrix(np.zeros((1, 3)), times).ravel() == pytest.approx([0.5, 1.5, 2.5])
    # a missing reading is bridged: same as integrating the valid points only
    assert aac_matrix(np.array([[np.nan, 0.0]]), [1, 2]).ravel() == pytest.approx([1.0, 1.0])


@pytest.fixture(scope="module")
def results():
    df = make(seed=1)
    gates = [(lo, hi) for lo in (0.05, 0.25, 0.5) for hi in (2.5, 4.0, np.inf)]
    return sweep(df, list(SEPARATION), "WT_DMSO", "F508del_DMSO",
                 gate_col="red_norm", gates=gates, level="plate")


def test_sweep_finds_the_planted_best_readout(results):
    assert best(results)["readout"] == BEST_READOUT


def test_gating_improves_quality(results):
    z3 = results[results.readout == BEST_READOUT]
    ungated = z3[(z3.gate_low == 0.05) & (z3.gate_high == np.inf)]["ssmd_paired_classical"].item()
    gated = z3[(z3.gate_low == 0.5) & (z3.gate_high == 2.5)]["ssmd_paired_classical"].item()
    assert gated > ungated


def test_paired_beats_unpaired_with_plate_drift(results):
    r = results[results.readout == BEST_READOUT]
    assert (r["ssmd_paired_robust"] > r["ssmd_robust"]).all()


def test_best_respects_cell_retention_constraint(results):
    with pytest.raises(ValueError):
        best(results, min_fraction_kept=1.01)


def test_grid_shape(results):
    g = to_grid(results, "zprime_robust")
    assert g.shape == (9, 6) and list(g.columns) == list(SEPARATION)
