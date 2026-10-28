# Research audit and evaluation protocol

1. Install package: `pip install -e '.[science,test]'`.
2. `python -m quintic_bcov_lab backend-check` confirms the pinned YY gauge and
   exact genus-2 anomaly equation.
3. `python -m quintic_bcov_lab generate-intermediate` produces 22 exact labels.
4. `python -m pytest -q tests` checks exact section/rational utilities and
   stage source-code identity.
5. `python -m quintic_bcov_lab prepare` re-creates one reference amplitude
   copy as symlinks in historical script paths. Shared sources are not training
   holdouts for late empirical models.
6. `python -m quintic_bcov_lab all --action test` runs all 18 historical
   test programs *in chronological order*. Reports are in `run_logs/`.
7. To inspect a stage, run its `--action run` command and then `--action test`;
   the original README documents the split between prediction and validation.
8. `python -m quintic_bcov_lab audit` flags any original file modified by a
   reproduction run (old snapshots can be overwritten in place).

### Research claims and their limits

- **Analytic:** YY ring closure, anomaly particular solutions, exact boundary
  constraints, low-degree GV vanishings, and the rational section transform.
- **Reference-derived:** the 22 residual coefficient targets and their
  zero-pivot affine representative at g2–12.
- **Empirical:** periodic phases, defect-front activity, affine inverse-envelope
  pulse profiles, and quadratic terminal-pulse fits.
- **Conditional comparisons:** the reference PT/DT predictions obtained from
  modularity / wall crossing in source archives. The empirical models do not
  ingest those formulas, but those mathematical assumptions remain embedded
  in the reference values.
- **Oracle-assisted diagnostics:** supplying an exact mature scale or a withheld
  PT anchor to reconstruct other coefficients is *not* scored as a blind
  non-oracle prediction.

### Dataset partitioning

Late-stage studies keep their `predictor/` (or top-level `predict_*.py`) and
`validation/`/`validator/` directories separate. In original stage test
programs, the predictor is executed with the latter directories excluded.
The unified runner uses these tests without weakening them. Do not import a
reference module from a new prediction model.

### Exact versus numerical arithmetic

- Use `fractions.Fraction` for YY coefficients, lattice steps, PT integers,
  and differences. The shared code **rejects float-to-rational coercion**.
- Use `numpy` and `sklearn` only for diagnostic continuous regressors, never
  to decide exact equality for a large GV/PT integer.
- Store rational numerator/denominator as base-10 strings in JSON/CSV.
- Record performance both for exact matches and for logarithmic/relative
  errors, explicitly reporting any abstentions.

### Repository layout

```text
src/quintic_bcov_lab/       Unified exact BCOV adapter and math utilities
src/quintic_bcov_lab/vendor/engine/  Pinned exact YY/BCOV backend
shared/reference_amplitudes.json   Single corrected exact reference, g2..12
shared/PROVENANCE.json      Source/seed/checkpoint and SHA256 information
studies/01...18/            Archived experimental sources and data
notebooks/overview_marimo.py     Interactive stage index
run_logs/                 Reproduction logs (generated on demand)
tests/                    Unified pytest suite
```

Optional genus-1000 exact dataset is distributed separately as noted in README.
