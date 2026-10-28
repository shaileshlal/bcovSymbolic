# Quintic BCOV intermediate-only ML experiment

Reproduces the original v0.1 anomaly-feature and simple sequence regression
baselines on **only** the undetermined YY ambiguity coefficients. Includes a
physics-motivated mixed-source magnitude normalization and a separate
within-genus shape experiment.

## Setup and run

From this directory:

```bash
python -m pip install numpy scikit-learn pandas matplotlib marimo
python run_experiment.py
python tests.py
marimo edit intermediate_ml_marimo.py
```

Exact rational section/boundary tests from the previous intermediate notebook:

```bash
cd source
python tests_intermediate.py
# Optional slower independent kernel solve through genus 9:
python tests_intermediate.py --boundary --max-genus 9
```

## Contents

- `Research_Report.md`: all inputs, ablations, forward-holdout results, caveats.
- `run_experiment.py`: reproducible scikit-learn benchmarks, exact rational label audit.
- `intermediate_ml_marimo.py`: interactive displays, model selection, plots.
- `tests.py`: exact label consistency and leakage/holdout checks.
- `results/benchmark.json`, `results/predictions.csv`: generated machine-readable results.
- `data/intermediate_training.csv`, `data/intermediate_exact.json`: the 22 exact labels and original residual-basis data.
- `source/`: reference YY amplitude archive and the exact residual/boundary implementation from the previous release, plus its independent tests. No changes to the reference data.

The reference-derived zero-pivot anchor must **not** be supplied to the predictor.
The learned models condition on a target-genus BCOV anomaly computed with the
true, already fixed lower genera; their successes do not solve the next physical
genus without additional input.

All ML predictions are approximate; exact rational reconstruction, conifold
gap, orbifold regularity and GV integrality are not replaced by these fits.
