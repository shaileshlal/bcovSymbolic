# Quintic PT near-extremal depth 11–15: recurrence stress test

**Experiment v1.1 · 28 October 2026**

## Executive summary

The empirical pulse model developed through shell depth 10 was subjected to a harder test at depths **11–15**. We reconstructed 65 conditional PT integers at degrees 41–53 and an independent 85-site reference at degrees 54–70. The **prediction model** uses exact-rational PT integer observations at degrees at most 40, the original period-four sign, the quadratic envelope, affine inverse-envelope pulses, and residue-pair antisymmetries. It does **not** call modular forms, DT coefficient generators, or wall crossing. The separate **validation** program does use the published rank-zero DT/wall-crossing framework and low-degree GV data, so all comparisons are conditional on that framework.

The original v1.0 affine-profile fit with no extra depth data has no identified stabilized scales beyond depth 10. We therefore supplied a fresh integer-only calibration table at depths 11–15, restricted to curve degrees <=40, without revealing the test degrees. The baseline and fitted models then gave:

| Model | Exact PT integers | Sites predicted | Abstentions | Median relative error on predicted sites |
|---|---:|---:|---:|---:|
| Linear extrapolation in k | 0 | 65 | 0 | 7.28% |
| Rolling quadratic in k | 0 | 36 | 29 | 1.24 × 10^-5 |
| Cubic interpolation | 0 | 6 | 59 | 1.82 × 10^-5 |
| Affine quadratic-envelope extrapolation | 0 | 65 | 0 | 1.73 × 10^-6 |
| Last normalized scale | 0 | 65 | 0 | 3.36 × 10^-8 |
| **Exact affine residual + paired residues** | **25** | **25** | **40** | **0** |

The last result is an **exploratory model**: residue-pair antisymmetries beyond the originally verified r<=2 were inspected in this research program. It is not 25 independently preregistered discoveries. With its parameters unchanged after this exploration, however, the model predicts **49/49** attempted exact integers on the subsequent degree-54–70 window (36 abstentions). This is a cleaner forward validation, albeit at sites closer to stabilization.

## 1. Representation

Write degree `d=5k+s`, `s ∈ {0,1,2,3,4}`, shell `t`, and

\[
W_{s,t}(k)=5k^2+(2s-5)k+c_s-2t, \qquad c=(10,6,4,4,6).
\]

Let `sigma_s(k)` denote the previously calibrated four-periodic sign. Then

\[
A_{s,t}(k)=\frac{\sigma_s(k)\operatorname{PT}(5k+s,t)}{W_{s,t}(k)}.
\]

Set `K_s(t)=t-1` for `s=0`, and `K_s(t)=t+3-s` for `s>0`; put `r=K_s(t)-k`. The working profile model is

\[
\boxed{A_{s,t}(k)=B_{\alpha(s),t}+(-1)^t\left(P_{s,r}+\frac{Q_{s,r}}{W_{s,t}(k)}\right),\quad 1\le r\le5.}
\]

There are three empirical scale groups: `α(0)=0`, `α(1)=α(4)=1`, `α(2)=α(3)=2`. The profile pairs show antisymmetry `(P,Q)_{0,r}=-(P,Q)_{2,r}` and `(P,Q)_{3,r}=-(P,Q)_{4,r}` through `r=5`. These are not assumed symmetries of the mirror quintic.

All parameters are exact `fractions.Fraction` objects. Before predicting a site, the solver reduces its linear profile equation against the calibration equations over Q. It emits an integer **only if the candidate is uniquely determined**. All other sites are marked `unidentified` or `outside_trained_profile`. There is no rounding and no optimization against held-out integers.

## 2. Reference and input provenance

The extended validator reconstructs the three rank-zero DT series through q-order 26 from the nine exact polar constraints in the existing validation code. It evaluates the PT/DT wall relation and all lower-degree wall states relevant to the tests, using independently known degree-one and degree-two genus-zero GV numbers `(2875,609250)` and the degree-three genus-one value `609250`. The degree-two PT function is independently assembled from the degree-two GV partition-function coefficient, including the double cover and disconnected degree-one contribution. Terms with a lower-degree PT input outside that list are marked unavailable; they are never set to zero.

For depths 11–15, every site at degrees 31–53 used by the calibration/holdout is calculable under the stated conditional reference. The new validator reproduces **116/116** previously certified shallow/deep PT reference values before the extension. Exact integer labels for degrees <=40 are copied into `data/observed_calibration_Q20_Q40_shell6_15.csv`. The learner never opens `validator/`.

Training data is *not* mathematically independent of the modular assumptions used to construct the conditional reference. The experiment is an evaluation of which features can **predict** that arithmetic without being told the corresponding analytic formulas.

## 3. Held-out results by shell

| Shell depth | Exact predictions / 13 targets | Eligible sites with r<=5 | Non-eligible sites (r>5) |
|---|---:|---:|---:|
| 11 | 11/13 | 13 | 0 |
| 12 | 9/13 | 12 | 1 |
| 13 | 3/13 | 10 | 3 |
| 14 | 2/13 | 7 | 6 |
| 15 | 0/13 | 4 | 9 |
| **Total** | **25/65** | **46** | **19** |

Example exact outputs (unseen by the training solver): `PT(42,11)=47080023927983572801070472750`, `PT(44,13)=-1254339017218216232117786519500`, `PT(45,14)=274572580194177671098994672175`.

On **degrees 54–70**, without updating any learned parameter or adding reference labels, the same solver identifies `17/17` integers at depth 11, `17/17` at depth 12, `11/11` at depth 13, and `4/4` at depth 14; it abstains on all 17 depth-15 sites. Thus **49/49 exact and 36 abstentions**. These later sites are more likely to lie at or near stabilization than the original deep-shell targets.

## 4. Structural profile tests and the first new obstruction

Using the *reference-derived* depth scales `B_{α,t}` only for **post-hoc structural diagnosis**, each of the 25 residue/front-distance families at `1<=r<=5` is precisely constant or affine in `1/W`, with at least one independent excess data point beyond the two parameters. The two residue antisymmetries also hold for all 25 profiles.

In contrast, at `r=6` the affine profile fails in residues `s=0,2,3,4`; only `s=1` still fits. The first explicit obstruction is at `s=0,r=6`: fitting the first two points (at envelope values 86 and 134) misses the third (at 192) by

\[
\boxed{\frac{18837667296875}{96}}\ne0
\]

in the signed normalized profile. This is a small *relative* defect (about `2.54 × 10^-9`), but an unambiguous obstruction to an exact integer formula. Larger `r` requires new correction terms or new shared scale sequences.

As an additional **conditional shape-only test**, we froze the 25 profiles inferred at depths <=15 and evaluated them at depths **16–20** at degrees 55–110. When the correct `B_{α,t}` values were **supplied by the validator**, the profile relation passed **123/123 exact** rational checks. This is evidence for a stable local *grammar*, **not 123 predictions of PT integers**: those new scales are oracles in this diagnostic.

## 5. Independent integer anchors and finite-window rank

Represent the unknowns by one scale variable `B_{α,t}` for each group and depth and two profile variables `(P_{s,r},Q_{s,r})` per independent residue/distance pair. Fix previously calibrated parameters from the v1.0 archive, and restrict the model to `r<=5` where exact calibration consistency has been verified.

The degree<=40 calibration equations have exact rational row rank **22**. Add the 46 eligible withheld sites as unlabeled linear functionals. Their span grows by only **eight** dimensions. One possible greedy set of integer anchors is

`(degree,depth) = (41,11),(43,11),(43,13),(44,14),(48,14),(45,15),(49,15),(53,15)`.

An explicitly **oracle-assisted**, separate sensitivity test supplies those eight reference integers and reconstructs all 46 eligible PT coefficients exactly. This is a meaningful finite-window compression: **46 sites require only eight new independent labels under the proposed profile law**. But it is not eight blind predictions, and it leaves 19 off-profile sites unresolved.

By depth, the needed new anchors in this particular ordered elimination are `2,0,1,2,3` at depths 11–15. Thus the count of new information does *not* decline monotonically to zero in this window.

## 6. Does the depth recursion terminate?

**Not established.** There are two different notions:

* **Local-profile closure:** a fixed simple rule correctly predicts the shape of corrections at any depth, conditional on a depth scale and front distance. For `r<=5` the evidence is strong, including 123 further conditional shape checks.
* **Global arithmetic closure:** all depth scales and all larger-r profiles are generated from finitely many initial integers by finitely many reusable recurrences. The current data do not establish this.

Post hoc, the three exact reference scale sequences `B_{α,t}`, examined at `t=0,...,15`, have Hankel matrices of maximal rank through size **8**. No exact *constant-coefficient* linear recurrence of order <=7 fits those 16 terms. We additionally tested polynomial-coefficient linear recurrences with orders 1–3 and polynomial degree <=2, and found no nonzero relation. These are only *limited model exclusions* from a finite sample; they say nothing decisive about higher-order or nonlinear/divisor-type recurrences.

In particular, the old affine grammar receives new shape defects at `r>=6`, while `B_{α,t}` contributes fresh exact values at successively larger depths. The tests do **not** show saturation in the number of independent arithmetic inputs. A finite nonlinear generation mechanism remains entirely possible, and would be a natural object to interpret post hoc in relation to modularity, but was not assumed or recovered here.

## 7. Suggested next empirical search

The most informative search is now an **exact sparse recurrence for the three scale streams** `B_{α,t}`, or a minimal residual update at `r>=6` using neighboring streams with shifted depth. Require a recurrence to predict multiple new scales exactly in forward holdouts, without importing them from the conditional reference. A useful complexity score is number of fresh integer seeds per 50 eligible sites, together with exact coverage and abstention rate. Numerical relative accuracy alone is insufficient, since several naive fits miss exact integers while achieving `10^-8` or better relative errors.

The script `validator/closure_tests.py` contains the exact finite-order nullspace tests, and `validator/profile_audit.py` and `validator/future_profile_test.py` are explicitly post-hoc reference-side diagnostics. The source code never silently treats unknown PT entries as zero.
