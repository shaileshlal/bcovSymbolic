# Quintic near-boundary PT pulses — empirical residual learner v1.2

**Results:** terminal pulses **18/18 exact** using 13 calibration integers; diagonal residuals **17/19 exact** on a fixed 24-site grid (5 unevaluable); conditional post-transient r=6 residuals **32/32 exact**. All comparisons are exact integer/rational equality, not floating-point tolerance.

## Run (standard library only)

```bash
python run.py
python validation/score.py
python tests.py
```

Optional interactive notebook:

```bash
pip install marimo matplotlib
marimo edit residual_marimo.py
```

## Contents

- `predict_terminal.py`: exact symbolic quadratic recurrence in terminal pulses; 13 calibration integers, 18 predictions.
- `predict_terminal_values.py`: next normalized residual `Z` from predicted pulse and three preceding observations.
- `predict_diagonal.py`: two-dimensional `r,ℓ` pulse-recurrence test (17/19 exact).
- `predict_r6.py`: transient-to-stable interpolation (32/32 on `Z`, conditional on the scale).
- `empirical_coordinates.py`: period-four sign and quadratic envelope from earlier exploratory runs.
- `data/`: **only model inputs** (no withheld answers, PT solver, or modular-form routines).
- `validation/`: independent reference and scoring, plus observed mature PT scale anchors; the optional DT/PT backend is strictly isolated from predictors.
- `results/`: exact predictions and audit.
- `RESEARCH_NOTE.md`: definitions, equations, score breakdown, failures, and limitations.

## Provenance / important limitations

The exact PT observations are **conditional** on the rank-zero DT modularity and wall-crossing reference from the previous version. The **empirical predictor** never uses those equations. Normalized pulse labels are formed exclusively using exact PT integers and *observed mature PT values at the same depth*, but these anchors sometimes occur at **larger curve degrees than the prediction targets**. Consequently the strong exact-pulse and `Z` scores are **not** causal predictions of full PT coefficients from earlier degrees alone. The small family of recurrence models was developed through exploratory analysis, not selected before any inspection of the reference structure. This is a numerical/arithmetic research baseline, not an all-genus theorem or complete BCOV ambiguity solver.
