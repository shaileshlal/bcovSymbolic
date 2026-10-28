#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python train_and_forecast.py
python event_models.py
python validate.py > results/validation_log.json
python validate_event_models.py > results/event_model_log.json
python future_forecast.py > results/future_log.txt
python plots.py
python tests.py
