# Quintic GV–YY two-stream reconstruction experiment

This is a self-contained exact-arithmetic physics audit and small-data ML ablation of the upper-degree YY ambiguity in the mirror quintic. Python 3.11+ and NumPy are sufficient (no GPU).

Run from this directory:

```
python run_experiment.py
python statistical_controls.py
python exact_response_audit.py
python verify_high_boundary.py  # exact high-X conifold/orbifold/MUM, g10–12; takes time
python high_degree_phase_study.py # historical stronger rolling-window baseline
```

Outputs:
- `results.json`: all exact MUM/GV checks, loss metrics, and test predictions.
- `statistical_controls.json`: historical-window sensitivity and a reproducible 199-permutation label-shuffle control.
- `exact_response_audit.json`: exact q Jacobian zeros and diagonal response values.
- `boundary_verification.json`: exact B-model high-X reconstruction at g10–12.
- `results_high_degree_phase.json`: historic inverse-genus baseline from the earlier high-X study.
- `Research_Report.md`: theory, experimental protocol, results, and limitations.

`quintic_ai_predictor/` contains the exact rational BCOV/period/GV library and corrected reference amplitudes g2..12. Training data are separated from evaluation by forward genus; predicting the *shape* of the high-X tail currently uses the true top coefficient at evaluation, explicitly described as an oracle. The exact boundary solver does **not** require that oracle or positive-degree current-genus GV inputs.
