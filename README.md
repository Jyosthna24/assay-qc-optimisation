# assay-qc-optimisation

![tests](https://github.com/Jyosthna24/assay-qc-optimisation/actions/workflows/tests.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

**Chooses image-analysis parameters by statistics instead of by eye. Every candidate configuration is scored by how well it separates a positive from a negative control (Z′-factor and SSMD), and the best one is kept, with the full grid as evidence.**

During my PhD (UCL, Vergani Lab) I used this approach to choose the **membrane zone** and **expression (red) thresholds** for a single-cell CFTR imaging assay, and the best **time point** for the functional readout. It builds on the lab's assay-optimisation framework and the published CFTR-SCAN pipeline. This repo is a general, open re-implementation that works with any per-cell table.

## The problem
A high-content imaging pipeline has free parameters, and each changes the result:
- *Which readout?* For example, how many pixels from the cell edge count as "membrane", or which time point to integrate a quench trace to.
- *Which cells?* For example, very low or very high expressers are often mis-measured.

If you tune these by looking at plots, the answer is subjective and can't be reproduced. If you tune them on the test conditions themselves, you risk overfitting. The principled alternative is to tune them **only on the controls**, by maximising assay quality.

## Method

For each combination of *gate* × *readout*, the controls are summarised at the chosen unit of analysis (plate means by default) and scored:

| Metric | Formula | Notes |
|---|---|---|
| Z′-factor | $1 - \dfrac{3(\sigma_P + \sigma_N)}{\lvert\mu_P - \mu_N\rvert}$ | ≥ 0.5 excellent, > 0 acceptable (Zhang et al., 1999) |
| SSMD | $\dfrac{\mu_P - \mu_N}{\sqrt{\sigma_P^2 + \sigma_N^2}}$ | effect size; ≥ 3 excellent for moderate controls (Zhang, 2007) |
| **Paired SSMD** | $\dfrac{\bar D}{s_D},\; D_i = P_i - N_i$ per plate | **the right estimator when both controls are on every plate** (Zhang, 2011) |
| Robust variants | median and MAD ($1.4826\cdot\text{median}\lvert x-\tilde{x}\rvert$) instead of mean and SD | one bad plate can't dominate the score |

**Why paired?** Plate-to-plate drift moves both controls together. The unpaired formula counts that shared drift as noise and **underrates good assays**. The paired formula works on per-plate differences, so the drift cancels out.

The best configuration is the one with the highest robust paired SSMD, **subject to keeping ≥ 50 % of cells**. Without that constraint, the search could "win" by throwing most of the data away. An optional Z′ floor is also available.

Also included: a vectorised, NaN-tolerant **area-above-the-curve (AAC)** calculation for quench time courses, so every time point can be swept as a candidate readout.

## Validation
The synthetic data generator plants a **known** best readout, a plate drift shared by both controls, and a population of mis-measured extreme expressers. The tests check that:
- the sweep recovers the planted best readout
- gating out extreme expressers improves SSMD
- paired SSMD beats unpaired SSMD whenever there is plate drift
- the robust Z′ survives a broken plate that drives the classical Z′ below zero
- every metric matches hand-calculated values, and MAD matches the SD on normal data

## Usage

```bash
pip install -r requirements.txt
python scripts/make_synthetic_data.py            # data/synthetic_cells.csv

python -m assayqc data/synthetic_cells.csv --pos WT_DMSO --neg F508del_DMSO \
    --readout-prefix zone_ --gate-col red_norm \
    --gate-lower 0.05 0.25 0.5 0.75 1.0 --gate-upper 2.5 3.0 3.5 4.0
```

```text
120 configurations scored (20 gates x 6 readouts)
Best: readout=zone_3  gate=[0.5, 2.5]  ssmd_paired_robust=... (excellent)  Z'=...  cells kept=82%
```

Outputs in `results/`:
- `sweep_results.csv`: every configuration, with all metrics (robust and classical, paired and unpaired), n and the fraction of cells kept
- `best_configuration.json`: the chosen configuration and its quality bands
- `heatmap_ssmd_paired_robust.png` and `heatmap_zprime_robust.png`: annotated grids with the best configuration outlined

Input: one row per cell, with `plate`, `well` and `condition` columns, candidate readout columns and (optionally) a gating column. The unit of analysis can be `--level plate | well | cell`.

## Skills demonstrated
Experimental design and assay validation · effect-size statistics · robust statistics · paired vs unpaired estimators · constrained hyperparameter search · avoiding overfitting by tuning on controls only · simulation-based testing · vectorised NumPy

## Acknowledgements
The membrane-proximity and AAC readouts, and the use of SSMD and Z′ to optimise them, come from the Vergani Lab's assay-development work at UCL ([CFTR-SCAN](https://github.com/vergani-lab-ucl/cftr-scan); Prins et al., 2020). This repository contains no lab data.

---
*Part of my PhD work on high-content CFTR screening. See my [GitHub profile](https://github.com/Jyosthna24).*
