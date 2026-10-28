# Quintic empirical PT scale-pulse learner v1.0

**Results:** **59 / 59 exact PT integers** predicted (6 abstentions) on the 65 depth-6–10 sites at degrees 41–53; the frozen previous v0.9 model predicted 32 / 32 (33 abstentions). Accuracy is with respect to a **conditional** DT/wall-crossing numerical reference, not an unconditional enumerative proof.

**Research hygiene:** The entire predictor uses PT integer calibration data at degrees <=40 and no modular forms, wall-crossing rules, or target PT values. The *feature family* (affine inverse-envelope residuals and anti-pair transfer) was chosen during exploratory inspection of the broader data, so the 59/59 retrospective benchmark is **not a pristine preregistered blind test**. The no-reference replay proves computational separation of calibration and validation files, not independence of prior model-selection decisions.

## Run with only standard-library Python

From this directory:

```bash
python predictor/pulse_profile.py
# Prints: PT predictions 59 abstentions 6 out of 65

# OPTIONAL: scorer reads validator/reference_Q20_Q53.csv
python validator/validate.py
python validator/one_seed_audit.py
python tests.py

# Optional plots
pip install matplotlib
python figures.py

# Optional notebook GUI
pip install marimo matplotlib
marimo edit pulse_marimo.py
```

`tests.py` constructs a temporary directory containing ONLY `data/` and `predictor/`, then reproduces the forecasts byte-for-byte without any reference files.

## Exact arithmetic and representation

Write `degree = 5k+s` and shell depth `t`, with `s=0,1,2,3,4`. The period-four sign `sigma_s(k)` is fitted from shallow PT labels. Define

```text
W_{s,t}(k) = 5k^2 + (2s-5)k + c_s - 2t
(c_0,c_1,c_2,c_3,c_4) = (10,6,4,4,6)
A_{s,t}(k) = PT(5k+s,t) / (sigma_s(k) * W_{s,t}(k))
K(s,t) = t-1 if s=0; t+3-s otherwise
r = K(s,t) - k
```

The three independent mature-scale classes are `{0}`, `{1,4}`, and `{2,3}`. The learned **pulse profile** is

```text
A_{s,t}(k) = B_{group(s),t} + (-1)^t [P_{s,r} + Q_{s,r}/W_{s,t}(k)]
```

All `B`, `P`, `Q`, and PT predictions are exact rational numbers / integers. No floats enter the predictor.

The source of each profile and the number of independent training support points are recorded in `results/fitted_pulses.json`. This is an actual fitting program, not a hardcoded table of the 59 answers.

## Empirical finding hierarchy

1. **Constant mature scale:** for `k >= K`, `A = B`. This is the frozen v0.9 result.
2. **Terminal pulse transfer:** the last active scale pulse satisfies `D_{s,t+1}(K+1) = -D_{s,t}(K)`. This increases exact predictions from 32 to 41 in an autonomous exact-propagation learner.
3. **Affine-in-1/W profiles:** at fixed distance `r`, fit `P + Q/W` to earlier residuals. A conservative version increases exact predictions to 48.
4. **Cross-residue anti-pair:** the `r=1` calibration reveals opposite pulse amplitudes in `s=0` versus `2`, and `s=3` versus `4`; the `s=0,2` depth-two profiles independently verify the sign relation, and a joint two-point depth-two fit transfers it to `s=3,4`. This plus one additional depth-four profile gives **59 exact predictions**.

For example, the learned `s=1,r=2` parameters are `P=-1135625` and `Q=28750`. At degree 41, depth 8, `k=8`, `W=286`, and `B=1867777794427065287875`:

```text
PT(41, depth=8) = 286 * (1867777794427065287875 - 1135625 + 28750/286)
                = 534184449206140347572250
```

This agrees with the withheld reference integer exactly.

## Secondary residuals

Define `D_{s,t}(k) = A_{s,t}(k) - A_{s,t}(k-1)` and `J_{s,t}(k) = D_{s,t}(k) + D_{s,t+1}(k+1)`. With `r=K(s,t)-k`, the fitted `P` terms cancel identically, leaving

```text
J_{s,t}(k) = (-1)^t [
  Q_{s,r}  (1/W_{s,t}(k)   - 1/W_{s,t+1}(k+1))
 +Q_{s,r+1}(1/W_{s,t+1}(k) - 1/W_{s,t}(k-1))
]
```

The current fitted `Q` profiles reproduce **14/14** testable nonzero `J` values in the reference. This is a derived consistency check on the affine hypothesis, not 14 additional independent measurements.

## Still unresolved

Six sites, all at depth 10:

```text
(degree,depth) = (42,10), (43,10), (47,10), (48,10), (52,10), (53,10)
```

The common mature scale for the pair `{2,3}` at depth 10 has no calibration anchor. One *withheld* integer, `PT(43,10)`, fixes that scale and then reconstructs four of the other five integers exactly; the remaining `PT(42,10)` is a genuinely new deeper-profile input in this empirical ansatz. `validator/one_seed_audit.py` demonstrates this **oracle sensitivity**, not an independent AI prediction.

## Directory structure

- `data/`: earlier PT labels and 65 unlabeled test sites; training max degree 40.
- `predictor/pulse_profile.py`: exact fitting and no-oracle predictions.
- `predictor/deep_envelope.py`, `predictor/recursive_defects.py`: reused frozen v0.9 empirical utilities.
- `results/pulse_predictions_Q41_Q53.csv`: all predictions and explicit abstentions.
- `results/fitted_pulses.json`: complete learned rational parameters and inference provenance.
- `validator/`: conditional reference labels, scores, front tests, closed secondary-defect checks and one-seed sensitivity.
- `pulse_marimo.py`: optional interactive notebook.
- `figures.py`: optional plots; `figures/`: generated images.

This is a restricted PT arithmetic experiment, **not** a complete holomorphic ambiguity reconstruction, and it does not certify the underlying DT/wall-crossing reference assumptions.
