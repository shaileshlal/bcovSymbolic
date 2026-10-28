# Quintic empirical pulse depth stress-test v1.1

**Task:** Test the previously developed near-extremal PT pulse recurrences at shell depths 11–15 and examine whether the empirical recursion might be finitely generated.

**Result:** On 65 PT targets at degrees 41–53, a calibrated exact-rational, affine-profile, paired-residue learner predicts **25/25 integers exactly**, abstaining on 40. On a separately scored forward-degree block 54–70 it predicts **49/49 exactly** and abstains on 36. A finite-rank audit shows that 46 modeled targets are generated from eight additional integer anchors under the working ansatz, while 19 lie beyond its profile range. There is no demonstrated finite-depth termination.

The learner receives only exact PT integer calibration data through degree 40, not modular forms or DT wall formulas. However, the calibration integers and the separately computed validation targets are **conditional on** a DT/modular and wall-crossing reference calculation. All mathematical assumptions and oracle-assisted diagnostics are disclosed in `RESEARCH_NOTE.md`.

## Reproduce

Install Python 3.10+ and SymPy (`python -m pip install sympy`). Predictor-only scripts use the Python standard library. From this directory:

```bash
python predictor/baselines.py
python predictor/linear_extension_pairs.py
python predictor/prospective.py
python predictor/anchor_rank.py
python validator/extended_reference.py
python validator/benchmark.py
python validator/profile_audit.py
python validator/future_profile_test.py
python validator/closure_tests.py
python validator/anchor_reconstruction.py
python tests.py
```

The modular reference generates `h_mu` through q-order 26 and requires SymPy. Its arithmetic check against the previous conditional reference is part of `tests.py`.

## Important folders

- `data/`: PT integers at degrees <=40 (calibration) and unlabeled target site lists.
- `predictor/`: phase, envelope, exact affine-pair solver, baselines, rank analysis; **never imports validator**.
- `validator/`: separate conditional DT/PT reference, exact checks, structural diagnostics. `anchor_reconstruction.py` uses eight **withheld oracle integers** and is explicitly not scored as a blind prediction.
- `results/`: predictions, calibration rank audit, finite-recurrence nullspace checks, metrics, plots.
- `depth_marimo.py`: interactive exploration; `pip install marimo matplotlib`, then `marimo edit depth_marimo.py`.

The flat `validator/reference_deep_11_15.csv` contains both calibration and held-out reference integers. Predictor scripts never read it. The extra degree-54–70 scores are generated separately by the benchmark validator and are not used to tune model parameters.
