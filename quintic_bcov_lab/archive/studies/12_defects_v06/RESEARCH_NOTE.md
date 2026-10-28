# A sparse-defect baseline for near-Castelnuovo quintic PT invariants

**Version 0.6 · 28 October 2026 · exploratory numerical/AI study**

## Summary

We investigate the **breakdown** of the unexpectedly accurate finite-difference laws found in v0.5. We use only the *integers* in previously tabulated PT invariants, their degree/shell indexing, and engineered features of their difference arrays. **No modular form, wall-crossing identity, rank-zero DT series, or BCOV boundary equation is a predictive input.** The published reference table itself was derived partly with modularity: source provenance must not be mistaken for independent physical evidence.

The main finding is that the third-difference defect is **sparse and temporally nested** in the two available forward-degree blocks. A model of the **location of the defect front** performs better than attempting to predict each enormous PT coefficient directly. A simple mirror-curvature correction is also a highly effective *numerical* model of the nonzero defect magnitudes.

## 1. Definition and exact arithmetic

Use degree `d=5k+s`, residue `s=0,1,2,3,4`, first physical PT shell `δ_s=(0,1,2,2,1)` and shell depth `t≥0`. Denote the PT integer by

\[
P_{s,t}(k)=\mathrm{PT}(5k+s,m_C(5k+s)+\delta_s+t).
\]

The period-four signs `σ_s(k)`, selected on degrees `≤20`, produce `T_{s,t}(k)=σ_s(k)P_{s,t}(k)`. All arithmetic involving PT data is exact Python integer arithmetic, including the differences

\[
C_{s,t}(k)=\Delta_k^2T_{s,t}(k),\qquad
\varepsilon_{s,t}(k)=\Delta_k^3T_{s,t}(k).
\]

Thus the exact recurrence is

\[
T_{k+1}=3T_k-3T_{k-1}+T_{k-2}+\varepsilon_{s,t}(k-2).
\]

The old v0.5 quadratic method sets the last term to zero. Our new study predicts the defect's (a) support, (b) size and (c) confidence of suppression.

## 2. Data and information barriers

We use a deliberately small, controlled sample: 24 residue–shell streams; first *physical* calibration values at degrees 6–20 (3 points per stream); the first genuinely forward block Q21–25 (24 integers); and a later Q26–30 block (24 integers).

- **Initial forecast:** `calibration_Q1_Q20.csv` → `block_Q21_Q25.csv`. The original v0.5 sign/quadratic/mirror-affine rule is fitted without Q>20.
- **Sequential update:** reveal the 24 observed integers at Q21–25, then forecast Q26–30 with a curvature/mirror rule. This is a **different information budget** from v0.5's fixed Q≤20 forecast of all 48 integers.
- **Supervised defect classifier:** 24 labels describing whether the third difference at Q21–25 is nonzero. Train a small classifier on these and predict the 24 defect statuses at Q26–30.
- **Validation:** all 24 Q26–30 target integers are in a separate `reference_only_Q26_Q30.csv` read **only** by `validate.py` and `validate_event_models.py`. A clean-copy reproduction runs the predictor **with that reference file entirely absent** and compares the saved forecast hashes.

These are retrospective degree-ordered holdouts, not prospectively registered experiments. Some architecture choices (notably the front-speed hypothesis) were developed while examining the earlier v0.5 results and these data. Such hypotheses need an independent Q31–40 test.

## 3. The defect-support geometry

Let `A_{s,t}(j)=1` when the next third difference is nonzero in block `j`. In the Q21–25 and Q26–30 windows:

| Residue `s` | Nonzero shells Q21–25 | Nonzero shells Q26–30 |
|---:|---|---|
| 0 | 4,5 | 5 |
| 1 | 1,2,3,4 | 1,2,3,4 |
| 2 | 1,2,3 | 2,3 |
| 3 | 2,3 | 3 |
| 4 | 3,4 | 4 |

There are **13 active locations** in the earlier block, **9** in the later block, and **no activations at a previously inactive site**. In particular, all 11 locations with zero observed previous defect remain exactly at zero in the next block. This is an exact descriptive statement for the two windows, not an all-genus theorem.

For each residue, the support forms an interval `t≥h_s(j)` truncated at the maximum available shell. Except for residue 1, the lower frontier moves by one shell. A **unit-speed moving-front hypothesis**, taking the lower frontier from Q21–25 and adding one, classifies **23/24** Q26–30 statuses correctly; its sole miss is residue 1, shell 1. The +1 rule is exploratory/post-hoc, not a validated extrapolation law. A conservative persistence rule (once defect zero, keep zero; otherwise remain active) gets **20/24**.

Under a purely artificial exchangeable-set null with 13 earlier and 9 later active sites, complete nesting has probability `binom(13,9)/binom(24,9)=0.0005468`. This is *not* a legitimate formal significance test for correlated arithmetic sequences; it quantifies only how restrictive the observed nesting would be under that null.

## 4. Engineered arithmetic/positional features

For each time-causal feature row, we compute:

- `(s,t,k_last,t−k_last)` and the earlier period-four phase;
- sign/log magnitude of the PT integer and of its observed second difference;
- the preceding third difference, if observable, and its scale relative to curvature;
- `v_2` and `v_5` of the observed curvature and PT integer;
- the paired-residue curvature at the same shell (`1↔4`, `2↔3`), its gap, and whether it brackets the last two local curvature estimates;
- an exact stabilization indicator `C_new=C_old`.

The **observed** (not predicted) defects show interesting extra arithmetic:

- All 9 defects that remain nonzero have nonincreasing *absolute sizes*: the median attenuation factor `|ε_new/ε_old|` is about **0.02165** (some families decay much faster).
- The defect **sign flips in 8/9** persistent streams.
- The `5`-adic valuation of the defect **never decreases** in these nine cases: 9/9 have `v_5(ε_new)≥v_5(ε_old)`.

These are useful candidates for sequence-model inductive biases, but the evidence is nine paired observations, and the patterns were examined after observing them. Their ability to predict the *next* nonzero defect remains untested.

## 5. A curvature-neighbour predictor of the defect magnitude

Let `C_old` and `C_recent` be the two latest *observed* second differences for `(s,t)`. Let `C_peer` be the latest second difference at its paired residue and the same depth.

We propose a simple deterministic estimator:

\[
\widehat C_{next}=\begin{cases}
C_{peer},& C_{old}\ne C_{recent},\quad C_{peer}\in [\min(C_{old},C_{recent}),\max(C_{old},C_{recent})],\\
C_{recent},&\text{otherwise}.
\end{cases}
\]

Then `\widehat T_next=2 T_last−T_penult+\widehat C_next`. This is equivalent to predicting a sparse defect correction `\widehat ε=\widehat C_next−C_recent`. The important feature is that **all three curvature inputs are observable before the new PT coefficient**. The chosen mirror pair and bracketing rule are structural hypotheses, not trained on the hidden integer.

### Exact-number benchmark

| Comparison | Train/update information | Evaluated PT integers | Exact | Mean relative error | Maximum relative error |
|---|---|---:|---:|---:|---:|
| v0.5 fixed sign/quadratic/mirror-affine | ≤20, no update | Q21–30 (48) | **28/48** | 0.289% | 8.504% |
| v0.5 first forecast only | ≤20 | Q21–25 (24) | **14/24** | 0.204% | 4.639% |
| Rolling quadratic | ≤25 | Q26–30 (24) | **15/24** | 0.008758% | 0.1994% |
| **Mirror-curvature defect predictor** | ≤25 | Q26–30 (24) | **17/24** | **0.0000650%** | **0.000741%** |
| **Sequential two-block protocol** | ≤20, then ≤25 | Q21–30 (48) | **31/48** | 0.1019% | 4.639% |

**Important:** the `31/48` score benefits from being allowed to observe Q21–25 before forecasting Q26–30. It is not a head-to-head replacement for the 28/48 *fixed-training* score. The stronger apples-to-apples result is the **second-block 17/24 versus 15/24**, and the reduction in numerical error on misses.

Across Q26–30 the mirror correction fixes two additional integers exactly, while bringing the worst relative error down by a factor ≈269 compared with the rolling quadratic. In particular, near-equality of a target curvature with its mirrored partner is a good predictive *numerical* feature even when it fails to produce the exact integer.

### Risk/coverage rather than forced exactness

When an observed third difference is already zero, continuing with zero is correct at **11/11** such sites in Q26–30. These can be labelled `stable_observed`; this is a high-confidence **empirical status**, not a mathematical certificate. A separate 90%-target split-conformal interval based on 24 Q21–25 errors has relative half-width **0.1575%** and covers all **24/24** Q26–30 integers; temporal correlations mean the formal iid conformal guarantee does not apply.

## 6. Small supervised classifiers: position wins over overengineering

The target is binary: **will the next third difference be nonzero?** Train on the 24 labelled events Q21–25; evaluate on the subsequent 24 events Q26–30. Fixed hyperparameters, no tuning on Q26–30.

| Classifier | Next-block accuracy | AUC | Precision | Recall |
|---|---:|---:|---:|---:|
| Positional logistic (residue, depth, depth-minus-k) | **21/24** | **0.970** | 0.800 | 0.889 |
| Logistic + log magnitudes, neighbouring curvature, `v₂,v₅` | 19/24 | 0.933 | 1.000 | 0.444 |
| ExtraTrees, engineered features | 19/24 | 0.889 | 0.700 | 0.778 |
| Random forest, engineered features | 19/24 | 0.859 | 0.700 | 0.778 |
| Persistent-front heuristic | 20/24 | 0.867 | 0.692 | 1.000 |
| Unit-speed-front heuristic (retrospectively proposed) | **23/24** | 0.944 | 1.000 | 0.889 |

The lesson is an important negative result: **additional arithmetic features have not demonstrated an out-of-sample advantage** here, even though some display intriguing exact regularities *after* the fact. With 24 labelled training events, feature proliferation readily overfits. The simplest positional classifier is the strongest genuinely trained model in this sample.

## 7. Prospective coefficients (unverified)

After incorporating the full published degree-30 table, the same time-causal mirror-curvature rule provides **candidate**, not independently validated, predictions for the next five degrees:

| Degree | Leading PT | Next PT | Following PT |
|---:|---:|---:|---:|
| 31 | −724,500 | −94,495,500 | −309,752,868,000 |
| 32 | −108,446,500 | −38,893,294,000 | +82,972,988,478,750 |
| 33 | −115,757,500 | −41,545,109,500 | +88,695,263,546,250 |
| 34 | −879,750 | −114,988,500 | −377,747,400,000 |
| 35 | −550 | −87,200 | −6,318,000 |

The known extremal degree-35 leading number `−550` is a useful external consistency check; the subleading numbers should not be regarded as verified. We **do not automatically iterate Q36–40** without observing the intervening degree block: doing so would silently replace missing observations with model predictions and compound uncertainty.

## 8. Recommended empirical program

The next version should target three outputs per `(s,t,k)` in a *hurdle* architecture:

1. **Defect activity:** sparse front/change-point prediction using `(s,t,k)` and recent stability, with calibrated abstention.
2. **Defect arithmetic:** sign/valuation and a **small integer residual** about the paired-residue curvature, not the full PT value.
3. **Exact PT output:** apply the finite-difference recurrence with the predicted defect using Python integers, not floats.

A mixture-of-recurrences or small discrete state-space learner is appropriate: two or three states (`transient`, `paired-curvature`, `stable`). This is a stronger small-data inductive bias than a generic transformer trained on 24 examples. A future transformer becomes sensible once independent higher-degree observations and more CY targets are available.

**The next decisive test** is a newly isolated degree Q31–35 block or data from a second Calabi–Yau, with all rules fixed before comparison. In particular test: (i) no reactivation, (ii) front velocity, (iii) sign reversal and `5`-adic monotonicity, (iv) cross-residue curvature repair, and (v) exact PT integers at deeper shells.

## Source and caveat

The source labels are a transcription of PT Table 3 from S. Alexandrov et al., *Quantum geometry, stability and modularity*, arXiv:2301.08066v2. The modularity-conditional provenance of the held-out data is not a modular input into any predictor. The experiment uses only tabulated integers and arithmetic operations. This is one quintic geometry and a small number of degree blocks; possible arithmetic identities should be considered **conjectures** until either independently reproduced or established symbolically.
