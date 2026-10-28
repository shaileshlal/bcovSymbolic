# Quintic empirical PT finite-difference defect baseline v0.6

**Purpose:** discover and predict the sparse failures of the quadratic finite-difference rule in near-Castelnuovo quintic PT integers. A pure-data exercise: **no modular forms, wall-crossing solver, DT generating series or BCOV equations are loaded by the predictors.**

## Quick start

```bash
python train_and_forecast.py    # predictions: Q21–25 from Q<=20, Q26–30 from Q<=25
python event_models.py         # small classifiers; needs numpy and scikit-learn
python validate.py             # reads withheld Q26–30 reference and scores numeric predictions
python validate_event_models.py
python plots.py                # needs matplotlib
python future_forecast.py      # unverified one-step Q31–35 candidates
python tests.py                # includes clean-no-holdout reproduction test
```

For interactive exploration: `pip install marimo pandas matplotlib scikit-learn` then `marimo edit defect_marimo.py`.

The **code executes the predictors in a clean directory with no holdout file present**, and its output hashes match the packaged predictions. All PT and finite-difference arithmetic is exact Python integer arithmetic; floating point is only used for classifier features and error diagnostics.

## Data boundaries

- `data/calibration_Q1_Q20.csv`: earlier PT integers, primary training data.
- `data/update_Q21_Q25.csv`: observed update for the second forecast, and labels for defect classifiers.
- `data/sites_Q21_Q25.csv`, `sites_Q26_Q30.csv`: unlabeled forecast-site coordinates.
- `data/reference_only_Q26_Q30.csv`: **reference-only** labels; only `validate.py`, `validate_event_models.py` read this file. The training scripts do not reference it.
- `results/block_*.csv`: exact integer predictions and causal model provenance.
- `results/event_risk_predictions_Q26_Q30.csv`: logistic, forest, and front activity-risk forecasts.
- `results/defect_analysis_only.csv`: labelled third differences (for diagnostics, **not model training**).
- `results/future_Q31_Q35.csv`: new, **unverified** candidate PT integers after observing all degree-30 data.

## Headline comparison

| Method | Forecast window | Exact PT integers |
|---|---|---:|
| v0.5 fixed training through Q20 | Q21–30 | 28/48 |
| v0.5 old first-block forecast | Q21–25 | 14/24 |
| Rolling quadratic after observing Q21–25 | Q26–30 | 15/24 |
| Mirror-curvature defect correction after observing Q21–25 | Q26–30 | **17/24** |
| Sequential update procedure | Q21–30 | **31/48** |

**Do not confuse** the update-assisted 31/48 with the fixed-data 28/48; their information budgets differ. The clean comparison is 17/24 vs 15/24 on the same later window with the same observed data.

Defect support is nested (13 active states become 9 with no reactivation). The moving-front *hypothesis* gets 23/24 defect statuses; a trained positional logistic classifier gets 21/24, outperforming larger engineered feature sets. See [RESEARCH_NOTE.md](RESEARCH_NOTE.md) for the full audit, all limitations and prospective tests.

Reference label source: Alexandrov–Feyzbakhsh–Klemm–Pioline–Schimannek, *Quantum geometry, stability and modularity*, arXiv:2301.08066v2 (published PT table). Its provenance is disclosed, and no modular theory is used in model fitting.
