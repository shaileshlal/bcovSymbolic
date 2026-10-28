# Quintic near-boundary PT: recursive defect-shell experiment, v0.7

Python-standard-library **exact-integer** prediction pipeline, continuing v0.6. The aim is to model **where quadratic finite differences fail**, not to use modularity, wall-crossing, or the DT generating function as predictive inputs.

## Run

```bash
python recursive_defects.py                 # train on degree <=25; predict 26-30, with no references opened
python validate.py                          # separate program opens the held-out degree-26-30 reference
python recursive_defects.py --calibration data/observed_Q6_Q30.csv --update data/empty.csv --sites data/sites_Q31_Q35.csv --out results/forecasts_Q31_Q35.csv
python tests.py                             # isolated no-oracle replay + exact identity checks
python plots.py                             # optional matplotlib/numpy
pip install marimo pandas matplotlib       # optional GUI
marimo edit recursive_marimo.py
```

## Key idea

For fixed residue `s`, shell depth `t`, and after removing the period-four sign, let `T(k)` be the PT integer sequence. Fit a degree-two polynomial `Q_2` to the first three observations. Define `eps(j)=Delta^3 T(j-3)` in the convention `eps(j)=T(j)-3T(j-1)+3T(j-2)-T(j-3)`. Then the identity

```
T(k) = Q_2(k) + sum_{j>=k0+3}^{k} eps(j) * binomial(k-j+2,2)
```

holds in exact integer/rational arithmetic. The empirically interesting part is that `eps` is sparse and often switches off; a finite-lived impulse induces a new **piecewise quadratic** branch. A global cubic fit treats the impulse as persisting and performs substantially worse.

## Strict evaluation

- Calibration observations degrees 6–20 plus extra low-degree reference rows: `data/calibration_Q1_Q20.csv`.
- Sequentially revealed degree 21–25 observations: `data/update_Q21_Q25.csv`.
- Unlabeled sites 26–30: `data/sites_Q26_Q30.csv`.
- **Only `validate.py`** reads `data/reference_only_Q26_Q30.csv`.
- All 24 streams have five positions in the degree 6–30 panel (3 calibration + 2 later). There are 11 with no observed third-difference pulses in either later block, 4 with one pulse, and 9 with pulses in both.

### Exact held-out matches, degrees 26–30

| Model | Matches / eligible sites |
|---|---:|
| Rolling quadratic | 15 / 24 |
| Global cubic | 11 / 24 |
| `A*k*Q_leading + B*Q_leading + C` | 10 / 19 |
| Paired-residue curvature / bracket rule | **17 / 24** |

The last method is carried over unchanged from v0.6. v0.7's main contribution is the precise pulse/locking interpretation and a separate negative control for polynomial-product extrapolation. **No new improvement over 17/24 is claimed.**

### Recursive shell observations

- Previously zero third difference stays zero in 11/11 observed next-block transitions.
- Defect support shrinks from 13 active sites to 9, with zero new activations.
- Pairing `s=1 ↔ 4` gives identical latest curvatures through `t=0` in the earlier block, through `t=1` in the later block; pairing `s=2 ↔ 3` advances from `t<=1` to `t<=2`.
- Persistent pulses show attenuation in 9/9 and nondecreasing 5-adic valuation in 9/9; signs reverse in 8/9. These statistics are retrospective and not guaranteed at the next degree.
- A zero-or-attenuation defect interval constructed from the preceding pulse covers all 24 held-out cases, but is not a calibrated probability interval and can contain very many integers.

## Interpretation and limitations

A local correction should usually be **a gated quadratic spline** `eps*binomial(k-j+2,2)*1[k>=j]`, rather than a global higher-degree polynomial. Predicting the future impulse `eps` remains the actual open problem. A promising small-data model has a sparse activity/front stage and then a residue-paired integer residual stage; keep arithmetic over Python integers/fractions.

The published PT integers used for labels have **modularity-conditional provenance**. None of those modular expressions, PT-DT wall-crossing constraints, or rank-zero DT coefficients is used by the prediction algorithm. Retrospective model/feature choices have been informed by this small dataset: this is an **empirical numerical baseline**, not a preregistered prospective test or a proof of all-degree recurrence.

`results/forecasts_Q31_Q35.csv` consists of one-step, **unverified** predictions using all previously observed degree<=30 integers; no degree 36–40 self-feeding is performed.

## Additional negative control

Run `python shell_transfer.py` to fit a `y^shift(A+B*y)` transfer between consecutive **defect shell polynomials**, conditional on the first two new defect coefficients. The deeper-shell checks fail for the pinned residue-1 family; see the last section of the research note. This is deliberately **not** counted as a blind integer benchmark.
