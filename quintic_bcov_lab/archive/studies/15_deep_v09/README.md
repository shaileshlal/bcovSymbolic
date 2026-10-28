# Quintic PT empirical deep-shell experiment v0.9

**Focus:** extend the frozen v0.8 recursive/finite-difference experiments from shallow shell depths to **t=6,7,8,9,10**, with curve degree **Q=31,...,53**. `Q` is curve degree, **not** genus. This is a reproducible, conditional-reference numerical experiment, not an independent enumerative theorem.

## Main results

- Rank-zero DT reference truncated at 14 q-orders and independent degree-one PT wall inputs reproduces all **109** previously released shallow reference values (Q=31,...,53, earlier shells).
- The deep-reference solver evaluates **106/115** candidate sites at Q=31,...,53. Nine sites at Q=31,...,40 require additional certified degree-2+ lower PT input and are **excluded**, never assigned zero. **All 65** sites Q=41,...,53 are available for evaluation.
- The **quadratic envelope** was inferred solely from shallow training integers (Q<=30):
  
  `W(s,t,k) = 5*k*k + (2*s-5)*k + c[s]-2*t`, `c=[10,6,4,4,6]`, `Q=5k+s`.
- After removing the period-four sign `sigma_s(k)` learned in v0.8, `A(s,t,k) = sigma_s(k)*PT(Q,t)/W(s,t,k)` has a sharp last possible scale-change position `K(s,t) = t-1` for `s=0` and `t+3-s` for `s=1,2,3,4`. This **pre-existing** diagonal-front rule predicts all **91/91** observed scale-change states at new depths (64/64 at Q>=41). The equivalent third-difference front matches **41/41** tests (39/39 at Q>=41).
- Ungated, a one-anchor envelope calibrated on Q<=40 yields **28/65** exact held-out PT integers; some misses have extraordinarily small relative errors but are *still incorrect integers*.
- A no-oracle **confidence-gated** predictor (only predict beyond the front using an already-mature same/paired-residue anchor) makes **32 predictions, all 32 exact**, and abstains on 33 not independently selected. No future data enter the prediction. This is **conditional on the learned front and paired-scale identification**; it does not fix the unsolved active pulse amplitudes.
- Of 24 persistent quotient-defect pairs, all decrease in magnitude, 19 reverse sign. In contrast, the prior 5-adic monotonicity is **not** universal at deeper shells.

## Package organization

- `predictor/`: pure integer/rational arithmetic, **no oracle**. `deep_envelope.py` generates forced forecasts, `gated_recurrence.py` abstains outside its learned stabilization regime, `recursive_defects.py` is the frozen v0.7/v0.8 predictor implementation.
- `data/`: shallow signed-phase calibration, earlier empirical PT data, **deep calibration only through Q=40**, and unlabeled target sites Q=41,...,53. No withheld PT labels are in the predictor input folder.
- `validator/`: conditional reference and exact metrics, derived from a rank-zero DT ansatz plus PT wall-crossing. These functions are never imported by predictor code. `shallow_reference_Q31_Q53.csv` is the independent earlier baseline for 109 exact cross-checks. `reference_Q20_Q53.csv` marks unavailable rows distinctly.
- `results/`: predictions, exact third/quotient-difference tables, metrics.
- `deep_shell_marimo.py`: interactive exploration of residue, depth, scale and discrete defects.
- `figures/`: defect-support and exact prediction plots.

## Run (Python 3.10+, SymPy for the validator)

From the extracted project directory:

```sh
# Purely empirical: NO validation/reference imports
python predictor/deep_envelope.py
python predictor/gated_recurrence.py

# Separate reference and scorer; SymPy needed here only
python validator/deep_reference.py
PYTHONPATH=validator python validator/audit.py
python tests.py

# Optional GUI
pip install marimo matplotlib
marimo edit deep_shell_marimo.py
```

**No-oracle replay:** `tests.py` also builds a temporary copy containing only `predictor/`, the calibration inputs, and the unlabeled target sites, and checks that the gated predictions are identical without the `validator/` directory.

## Provenance and caveats

The conditional exact reference uses the rank-zero DT modular series and the wall-crossing PT relation of Alexandrov, Feyzbakhsh, Klemm, Pioline and Schimannek: <https://arxiv.org/abs/2301.08066>. Neither the DT modular ansatz nor wall crossing is a model feature. The degree-one PT input follows from the independently known 2875 quintic lines, `PT_1(q)=2875*q/(1+q)^2`.

The front rule is **not a new out-of-the-blue machine-learned law**: it was extracted retrospectively in the previous v0.8 shallow benchmark. Its success at *new depths* is a real conditional cross-depth test, not a proof for all depths or degrees. The true active-zone pulse **values** are largely unresolved by the current empirical model.
