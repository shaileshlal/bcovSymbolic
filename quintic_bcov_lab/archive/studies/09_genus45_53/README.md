# Quintic GV displacement: genus 45–53 bootstrap and PT/DT periodicity

An exact-arithmetic proof of concept in the original Yamaguchi–Yau quintic gauge.

## Two distinct experiments

1. `code/bootstrap_45_53.py` computes the *169 residual GV charge slots* at genera 45–53, using only the proven Castelnuovo vanishing and the proven extremal genus formula. Exactly one of these charges is nonzero: `n[51,20]=175`. The value 175 is produced from the projective-bundle Euler characteristic formula, not read from the target genus. It is an **analytic/geometric** prediction rather than a novel machine-learning discovery.
2. `code/near_extremal_benchmark.py` trains on just three early, independently known values in each `degree mod 5` family and predicts two subsequent values per family. The integer-valued predictor factors into a quadratic magnitude envelope and a phase of period four. Its **8/8 exact held-out matches** are to GV values in Table 2 of arXiv:2301.08066v2, which **assumes modularity**. It predicts `n[54,21]=319125` without reading that held-out label.

A third component `code/bcov_sensitivity_51_53.py` actually propagates the *difference* caused by `n[51,20]=175` through the exact BCOV/YY anomaly equations and B-model boundaries to genera 52 and 53. Its full rational polynomials are in `data/delta_*.json`. This is **not a full recomputation of physical F_45,...,F_53**, because the current exact high-genus full amplitudes were not supplied in the available repository. It is a certified nontrivial sequential sensitivity calculation.

`code/pt_dt_arithmetic.py` tests the GV-to-PT and rank-one DT leading coefficients over integers. PT/rank-one DT predictions are algebraic transforms of GV, so agreement does **not** give an independent lattice constraint; the rank-zero DT modular form is a separate additional input.

## Run (Python 3.10+, no nonstandard dependencies for exact tests)

```bash
python code/bootstrap_45_53.py
python code/near_extremal_benchmark.py
python code/pt_dt_arithmetic.py
python code/bcov_sensitivity_51_53.py
python code/tests.py
```

The optional interactive notebook (`displacement_marimo.py`) requires `pip install marimo pandas matplotlib` and can be opened using `marimo edit displacement_marimo.py`. This notebook was syntax-checked but not launched in marimo.

## Outputs

- `data/sequential_g45_g53.csv`: each genus/degree residual charge, with its derivation.
- `data/bootstrap_results.json`: 45–53 charge count and exceptional degree-20 charge.
- `data/residue_family_benchmark.csv`: *training/holdout provenance*, predictions and baselines.
- `data/near_extremal_results.json`: exact fitted integer model, scores, and the genus-54 candidate.
- `data/pt_dt_checks.json`: low-order PT/rank-one DT integer responses and checks.
- `data/delta_P51_exact.json`, `data/delta_P52_exact.json`, `data/delta_P53_exact.json`: normalized sparse-polynomial sensitivity in YY variables `(u,v1,v2,v3,X)`.
- `data/delta_R52_exact.json`, `data/delta_R53_exact.json`: the exact anomaly-source changes.
- `data/physical_sensitivity.json`: integration and conifold/orbifold/MUM validation audit.
- `data/reference_P2.json`: the small *exact* g2 seed used to produce the genus-53 source.
- `RESEARCH_NOTE.md`: mathematical analysis and caveats.

All data are exact integers or numerators/denominators encoded as strings. No target-genus 54 label is needed to fit the near-boundary model.

## Source/assumptions

- Zhiyu Liu & Yongbin Ruan, *Castelnuovo bound and higher genus Gromov-Witten invariants of quintic 3-folds*, https://arxiv.org/abs/2210.13411. Proven bound, extremal curve result.
- Alexandrov, Feyzbakhsh, Klemm, Pioline & Schimannek, *Quantum geometry, stability and modularity*, https://arxiv.org/html/2301.08066 (v2, Tables 2 & 3, §5.2). The high-genus entries of Tables 2 & 3 explicitly **assume modularity**.

**Not established:** an all-genus theorem for the period-four *near*-extremal formulas, independent validation of the published genus-54 integer from a new BCOV run, or self-sustained complete free energies through genus 53.
