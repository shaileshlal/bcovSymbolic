#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python observations_raw.py
python prepare_splits.py
python run_experiments.py --train data/calibration.csv --sites data/sites_unlabeled.csv --out results/pt_predictions.csv
python finite_difference_fingerprint.py
python fit_empirical_transform.py
python validate.py
python validate_dt_fingerprint.py
python future_forecast.py
python tests.py
if python -c 'import numpy, sklearn' 2>/dev/null;then OPENBLAS_NUM_THREADS=1 python ml_ablation.py;fi
if python -c 'import matplotlib' 2>/dev/null;then python plots.py;fi
