# Quintic near-extremal arithmetic: empirical recurrence experiment

**Date:** 28 October 2026. **Model inputs:** only published integer PT/GV values at degrees 1–20. **Forward holdouts:** the published degree 21–30 integers, revealed only by independent validation. **No modular forms, DT Fourier series, or PT/DT wall-crossing relations are called by the fitted predictor.** The reference table *used for evaluation* is itself modularity-conditional, a crucial limitation.

## Main outcomes

- Smallest common sign period found on degree-1–20 calibration: **4** in `k=floor(degree/5)`.
- 48/48 held-out PT **signs** predicted correctly; GV signs: 43/48 by analogous rule.
- Independent signed-quadratic PT fit: **22/48** exact held-out integers.
- Mirror-residue affine transfer + signed-quadratic fallback: **28/48** exact integers (6/6 of the transferred coefficients); median relative error over the 20 misses: **0.0181%**.
- Learned local PT→GV arithmetic recurrence: **26/30** exact GV integers in first three relative shells (not independent from PT).
- PT second differences reproduce **nine** leading rank-zero DT q-coefficients *using training PT only*; three additional coefficients initially estimated to within 0.04%, 0.09%, 0.19%, and become exact once withheld degree data are included.
- With degree-21–30 integers subsequently added as calibration, **38/48** degree-31–40 near-boundary PT forecast sites meet a *stabilized-curvature or mirror-identity* criterion. These are **unverified forecasts**, except degree 35 and 40 leading terms, which independently agree with the established extremal formula.

## Reproduce (Python ≥3.10)

Core numerical/arithmetic scripts use only the standard library. Install `numpy`, `scikit-learn`, and `matplotlib` to reproduce the ablations and figures, or marimo to open the notebook.

```bash
python observations_raw.py                     # re-typesets source tables into CSV/JSON
python prepare_splits.py                       # isolated calibration/holdout sites
python run_experiments.py --train data/calibration.csv --sites data/sites_unlabeled.csv --out results/pt_predictions.csv
python finite_difference_fingerprint.py        # train only; DT reference NOT read
python fit_empirical_transform.py              # robust arithmetic PT/GV corrections
python validate.py                             # FIRST reads reference-only Q21..30
python validate_dt_fingerprint.py              # FIRST reads independent DT list
python future_forecast.py                      # extrapolate stabilized rows Q31..40
python tests.py                                # clean isolated no-oracle smoke test
python ml_ablation.py                          # requires scikit-learn + numpy
python plots.py                                # requires matplotlib
marimo edit eda_marimo.py                      # optional GUI
```

## Important file boundaries

- `data/calibration.csv`: predictor-visible degrees 1–20, including a few partial earliest rows.
- `data/sites_unlabeled.csv`: degree/shell pairs at Q21–30, **no integer labels**.
- `data/reference_only_Q21_Q30.csv`: physically withheld from `run_experiments.py`.
- `results/pt_predictions.csv`: exact integer predictions, with selected model for each site.
- `results/scored_predictions.csv`, `results/gv_predictions.csv`, `results/metrics.json`: post-hoc holdout scoring.
- `results/calibration_dt_fingerprint.csv`: PT-derived curvature coefficients, no DT input.
- `results/dt_fingerprint_comparison.csv`: DT correspondence validated after predictions.
- `results/forward_Q31_Q40.csv`: exploratory recurrence continuation; blanks indicate no stabilized pattern.
- `RESEARCH_NOTE.md`: equations, empirical findings, error analysis, limitations.
- `figures/`: exact-hit counts and finite-difference stabilization.

## Provenance and epistemic status

Transcribed (using integers rather than floating-point approximations) from Tables 2 and 3 of S. Alexandrov, S. Feyzbakhsh, A. Klemm, B. Pioline, T. Schimannek, *Quantum geometry, stability and modularity*, arXiv:2301.08066v2, https://arxiv.org/html/2301.08066. The paper's Q>20 integers depend on its modularity analysis; we use them solely for *evaluation*, not for *features, priors, loss functions, or fitted parameters*. The previous project work influenced the choice of residue-phase and difference candidates, so this is not a historically preregistered model discovery. No theorem is claimed from a finite number of exact matches.
