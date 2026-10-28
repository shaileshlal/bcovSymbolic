# Quintic empirical GV/PT/DT prediction, exact solver v0.4

**Focus:** Can small DT/PT patterns predict the unknown quintic genus-54 GV integer 319125, and how far can they extrapolate?

**Short answer:** The charged-residue learner trained on PT degrees 6–10 predicts all degree-11–30 leading PT integers correctly. More strongly, exact modular reconstruction from seven polar rank-zero DT coefficients + wall crossing predicts 319125 and the next four PT values at degree 21 with no degree-21 labels in the fit. Exact PT/GV inversion then recovers five GV values including the first disconnected-pair correction. See `RESEARCH_REPORT_v04.md` for qualifications.

## Install and run

```bash
pip install sympy numpy scikit-learn
python modular_core.py
python empirical_benchmark.py
python period_scan.py
python model_ablation.py
python pt_from_dt_solver.py
python pt_to_gv.py
python bcov_exact_response.py
python tests.py
```

Optional interactive exploration: `pip install marimo pandas` then `marimo edit dt_pt_marimo.py`.

## The programs

- `modular_core.py`: exact rational 9×9 solver for rank-zero DT series from seven polar coefficients and two allowed-exponent constraints.
- `empirical_benchmark.py`: PT-leading residue/kinematic experiment; blind-within-protocol integer predictions, future forecast Q31–40.
- `period_scan.py`: rolling-origin training-only modulo scan.
- `model_ablation.py`: naive ML and interpolation baselines.
- `pt_from_dt_solver.py`: exact bounded wall enumeration and DT→PT for Q=21, five coefficients.
- `pt_to_gv.py`: binomial/plethystic conversion, including the \(175\times2875\) product at fifth order.
- `bcov_exact_response.py`: exact genus-54 BCOV homogeneous lattice perturbation, PF/MUM response.
- `tests.py`: reproducibility assertions, target label only in checks.
- `results/`: machine-readable benchmark scores, exact predictions, and forecast files.
- `data/residue_family_benchmark.csv`: reference data from the previous project package, training vs modular comparison-only provenance attached.

**Not claimed:** Complete BCOV genus-54 free energy, an all-genus proof of the initial empirically selected period, an independent physical confirmation of 319125 using BCOV alone, or a prospective discovery predating the published modular tables.

### Isolated no-oracle run

The following physically separates the degree-21–30 comparison labels from all files opened by the predictor:

```bash
python prepare_blind_splits.py
python blind_predict.py
python validate_blind.py
```

Only the last command reads the comparison-only labels. The input SHA-256s and predicted-file hash are in `results/isolated_run_manifest.json`.

`python extremal_check.py` independently confirms the Q=35 and Q=40 prospective predictions against the known d=5k geometry.
