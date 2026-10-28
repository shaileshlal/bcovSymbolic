# Quintic PT recursive-defect survival, v0.8

**Date:** 28 October 2026. **All degree indices below refer to curve degree, not genus.**

This evaluates the **unchanged v0.7 numerical PT predictor** on degrees **31–53**. The predictor uses *only observed PT integers through degree 30*, the learned four-periodic phase, local quadratic finite differences, and paired-residue curvature. It does **not** call modular forms, DT coefficients, or wall-crossing.

Two forecasts are produced:

1. `autonomous_predictions_Q31_Q53.csv`: predicted degrees 31–53, appending **only the model's predictions** after every five-degree block (no true future labels at all).
2. `teacher_forced_predictions_Q31_Q53.csv`: predict each five-degree block before revealing its true values, then allow those previous values into the next block. This is explicitly a rolling-update benchmark, **not** blind 23-degree extrapolation.

The reference is generated **separately** by the quintic rank-zero DT modular ansatz followed by the applicable wall-crossing relation, as in Alexandrov–Feyzbakhsh–Klemm–Pioline–Schimannek (arXiv:2301.08066). It is therefore *conditional on those assumptions*. This theory enters **only the validator**. The reference module first reproduces 48 previously published PT values at degrees 21–30, then checks four extremal degree-35/40/45/50 leading values against Liu–Ruan's independent geometric formula (arXiv:2210.13411).

## Results

| Degree block | Rolling update (exact/total) | Autonomous (exact/total) | Nonzero third-difference defects |
|---|---:|---:|---:|
| 26–30, earlier window | — | — | 9 |
| 31–35 | 22/24 | 22/24 | 4 |
| 36–40 | 23/24 | 22/24 | 2 |
| 41–45 | 24/24 | 22/24 | 1 |
| 46–50 | 24/24 | 22/24 | 0 |
| 51–53 | 13/13 | 11/13 | 0 |
| **Total 31–53** | **106/109** | **99/109** | **7 actual nonzero pulses** |

On the rolling-update protocol, quadratic, cubic and paired-curvature models recover 102/109, 93/109, and 106/109 exact PT integers respectively. The `linear*quadratic` candidate recovers 79/86 of the eligible cases (ineligible leading shells not counted).

All 109 new reference sites satisfy the domain of the wall-crossing formula. Only one possible non-vacuum wall candidate occurs in the tested strip, at degree 31 shell 4, and its lower PT coefficient vanishes; hence all evaluated reference values are effectively from single rank-zero DT contributions. **This makes the test especially favorable to quasi-polynomial extrapolation; deeper PT shells are not covered.**

## Reproduce

Requires Python >=3.10 and `sympy`. `matplotlib` is optional for figures.

```bash
# Predictor (does not read or import the validation reference):
python predictor/rollout.py

# Independent reference and scoring:
python validator/generate_reference.py
python validator/score_and_front.py
python make_figures.py
python tests.py

# Optional GUI:
pip install marimo pandas
marimo edit survival_marimo.py
```

`predictor/recursive_defects.py` is a byte-for-byte copy of v0.7, SHA-256 `e180b2d314fe91111c0ea9dabdc5aeb9472c375708ea29903f0d99ac697223d0`.
The isolated autonomous predictor can be run with the entire `validator/` folder and `results/reference*` files absent.

## Files

- `predictor/`: frozen predictor and pre-31 inputs. No reference code.
- `validator/`: independent DT-derived reference, exact integer scoring, pulse audit.
- `results/`: predictions, reference, exact misses, defect front, geometry spot checks.
- `figures/`: rolling/autonomous accuracy and defect support.
- `RESEARCH_NOTE.md`: mathematical summary, negative controls and interpretation.

## Limits

These results establish **survival of a particular arithmetic prediction rule against a conditional reference**, not an unconditional theorem, not uniqueness of the topological-string ambiguity, and not predictive performance across every PT shell at degree 53. Structural features were chosen in earlier retrospective analysis; no retuning was allowed after v0.7. Only the shallow 24-per-five-degree near-extremal shell strip is tested.
