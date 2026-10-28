# Quintic residual ambiguity dataset and marimo notebook

## Start

From this directory:

```bash
python -m pip install marimo matplotlib
marimo edit intermediate_marimo.py
# or: marimo run intermediate_marimo.py
```

For a non-interactive data rebuild, **only Python's standard library** is needed:

```bash
python intermediate.py --source data/reference_amplitudes.json --output data
python tests_intermediate.py
python tests_intermediate.py --boundary --max-genus 9
```

The checked-in CSV and JSON have already been generated, so you can use them without installing marimo.

## What the coefficients mean

We use `x=psi^(-5)=3125*z_PF` and `X=x/(x-1)`, in the **original YY gauge**. The corrected genus-four reference is used.

At genus `g>=2`:

```
K=3*g-3
m=ceil(K/5)
D=floor((2*g-2)/5)=g-1-m
physical: a_g(X)=sum_{j=0}^K a[g,j] X^j

physical = anchor + delta
anchor: a^*_{g,j}=0 at j=m+1,...,g-1
homogeneous residual: delta(X)=X*(1-X)^m*sum_{r=0}^{D-1} c[g,r]*X^r
```

Restore the reduced one-dimensional section factor (suppressing the common Hodge-period factor):

```
f_g(x) = (x-1)^(g-1)*a_g(x/(x-1))
delta f_g(x) = sum_{r=0}^{D-1} b[g,r] * x^(r+1)*(1-x)^(D-r-1)
b[g,r] = (-1)^(g-2-r)*c[g,r]
           = x*Q_g(x), deg(Q_g)<=D-1
```

The homogeneous kernel is fixed by conifold gap, orbifold regularity, and equality of the MUM constant-map term. The `anchor` here is a **reference-derived** canonical representative selected by setting the intermediate pivot entries to zero; it is **not** a separate solution of those boundary conditions from independently known source data at every genus. As a result `anchor` and `a_pivot_phys` must not be provided as test-time ML features for a blind prediction. For blind inference, instead compute the zero-pivot anchor from the BCOV anomaly and the physical boundary RHS independently.

In the **full physical section** the period/Hodge-line frame has to be restored. `x -> 1-x` is a coordinate transformation, not a quintic involution symmetry giving additional physical equations.

## Provenance

The current user-named `bcov0.2` package is the **dual-corpus structural-learner v0.2** (`bcov_dual_corpus_v02.zip`), which has `tasks.jsonl`, a `structure_registry.json`, and separate BCOV/Lean corpora. This library archive could not be materialized during this run. More importantly, that corpus is not itself a verified physical genus-by-genus amplitude provider.

This package therefore uses the accessible exact **BCOV YY amplitude archive** from the earlier quintic predictor/midpoint experiment, including `reference_amplitudes.json`, and the unchanged exact `bcov_yy.py`, `physical_boundaries.py`, `period_frames.py`, `exact_series.py`, and `hybrid/boundary.py` to verify the known boundary kernel. The marimo notebook has an optional v0.2 path checker (`./bcov_dual_corpus_v02`) for structural-corpus integration without overwriting exact physical targets.

The source references contain exact rational nonholomorphic polynomial terms and physical pure-X ambiguity vectors for genera 2–12. The boundary-gauge kernel has been independently checked with the exact conifold/orbifold/MUM linear solver for genera **2–9** here; the exact 11-genus basis conversion was checked for all genera **2–12**. This release does not claim it ran the expensive full boundary solve at genera 10–12 anew.

## Files

- `intermediate_marimo.py`: reactive, selectable-genus **marimo notebook** that re-runs the exact transformation and regenerates the dataset.
- `intermediate.py`: standard-library arithmetic implementation, JSON/CSV exporter, v0.2-manifest discovery.
- `data/intermediate_training.csv`: **22 rows**, one per intermediate coefficient across 11 genera; includes exact rational labels, YY pivot, section-basis label, and simple mixed-anomaly features. Numerical signed-log fields are visualization/training features only; all rational labels remain exact strings.
- `data/intermediate_exact.json`: per-genus full reference, canonical zero-pivot anchor, homogeneous residual, `X`-basis and `x(1-x)`-basis vectors, with all rationals as strings.
- `data/validation.json`: 11 dimension/round-trip records.
- `data/reference_amplitudes.json`: corrected YY reference inputs, including exact anomaly polynomials (no floats).
- `quintic_ai_predictor/`: **unmodified** exact BCOV and physical boundary modules copied from the previous reference package, for independent audits.
- `tests_intermediate.py`: exact rational tests and optional physical boundary-system checks.

## Sample counts

| Genus | Residual parameters |
|---:|---:|
| 2 | 0 |
| 3 | 0 |
| 4 | 1 |
| 5 | 1 |
| 6 | 2 |
| 7 | 2 |
| 8 | 2 |
| 9 | 3 |
| 10 | 3 |
| 11 | 4 |
| 12 | 4 |

Total: **22** supervised intermediate directions, 11 genera.

## Limitations

The package contains a small **exact experimental dataset**, not a new AI model or a prediction of unknown higher-genus invariants. In particular, using `anchor` derived from `physical` would introduce label leakage. The difference between reference-derived coordinate labels and independently computed boundary-fixed features is explicitly preserved throughout. `marimo` was unavailable in the generation container (package installation was blocked); the notebook was syntax-compiled, all 16 cells were exercised with a simple API stub, and its computational engine was executed and independently verified. The first actual marimo GUI run should be performed in an environment with `marimo` installed.
