# Quintic residual GV displacement — exploratory study v0.3

Self-contained machine-readable lattice positions through genus 1000, exact mathematics, exploratory feature scans, integer-label holdouts, and a marimo notebook.

## Quick start

```bash
python -m pip install numpy scikit-learn sympy matplotlib
python build_data.py
python eda.py
python learn_top.py
python permutation_control.py
python exact_solver.py
python neighbor_regulator.py
python frontier_solver.py
python tests.py
```

For marimo: `python -m pip install marimo` then `marimo edit displacement_marimo.py`. The notebook requires precomputed data/results (already bundled). The source was syntax-checked; marimo itself was not installed in the execution environment.

## Dataset and targets

- `data/all_lattice_sites_g2_g1000.csv.gz`: 199,400 residual lattice sites, with exact geometry features, three-way status and **blank** displacement for unknown values. The file name denotes scope genus ≤1000; there are no residual slots at genus 2 or 3.
- `data/certified_labels.csv`: 63,560 geometric labels — **63,544 zero values and 16 known nonzero extremal values**. Do **not** interpret this as 63,560 independent measurements.
- `data/extremal_displacements.csv`: the 16 exact nonzero labels, with train k=4..11 and test k=12..19.
- `data/exact_g4_g12_calibration.csv`: 22 exact lower-genus GV residual checks.
- `data/frontier_g49_g80.csv`, `data/frontier_status_summary.csv`: transparent look at the transition to unprotected charges.
- `results/model_benchmarks.json` and `results/heldout_predictions.csv`: genuine forward-index holdout scores on the *known extremal family*. These data are generated from the established extremal theorem; success is not a new unknown-GV prediction.
- `results/frontier_predictions.json`: exact g51 and structurally determined g54 section responses, with **abstention** on unknown charges.
- `results/neighbor_regulator.json`: exact mirror-map-derived nearest-neighbor transport verified against prior BCOV data.

A solver predicting an integer vector `n` must first apply the proven vanishing and extremal constraints, then use `exact_solver.py` for `delta_c = H^{-1} delta_n`. The full H matrix is bundled only for genera 4–12. For g=51/54 we derive the terminal-direction response analytically; we do **not** claim to have a full g51/g54 BCOV boundary particular.

## Main observations

1. The first nonzero residual charge is **n(51,20)=175**; all positions up to genus 53 are fixed by geometry. The first unresolved value is **n(54,21)**.
2. For extremal d=5k, signs are period-4 in k and magnitudes quadratic. A symbolic search and a 2D Fourier-phase model each get 8/8 withheld extremal integers, while naive black-box regressors get 0/8 exactly. This rediscovers a known extremal curve theorem, not physics beyond it.
3. Any last-free-degree displacement has `delta f_g = 5^(4g-4-5D) * n(g,D) * x^D`.
4. Its next-degree transport is `delta n(g,D+1) = (240(g-1)-770D) n(g,D)` for D≥2, with the known D=1 twice-cover correction. All nine bundled exact cases pass.
5. **n(54,21) is still unknown.** The nearest-height choice zero is an unverified fallback, and the model should abstain rather than claim to know it.

See `RESEARCH_REPORT.md` for methods, EDA, tables, derivation and caveats.

Sources: https://arxiv.org/abs/2210.13411 and historical BCOV exact audits bundled under `source/`.
