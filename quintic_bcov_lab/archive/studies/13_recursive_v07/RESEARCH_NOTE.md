# Beyond quadratic finite differences: localized residual pulses and a shell-locking front

*Quintic PT empirical study v0.7 — 28 October 2026*

## 1. The hypothesis

For curve degree `d=5k+s`, `s=0,…,4`, and shell depth `t`, remove the empirically fitted four-periodic sign and set

\[
T_{s,t}(k)=\sigma_s(k)\,\mathrm{PT}(5k+s,m_0(5k+s)+\delta_s+t),\qquad
\varepsilon_{s,t}(k)=T(k)-3T(k-1)+3T(k-2)-T(k-3).
\]

Here `t` is the offset from the first physical PT shell. The leading quadratic model is `eps=0`. Testing the user's suggested `linear × quadratic` or unrestricted cubic correction amounts to assuming the latest nonzero third difference persists. We test this, then make the correction **localized in k and structured across t**.

## 2. A discrete Green-function identity (exact, not a fitted result)

Fix the first position `k0` with three consecutive observations and let `Q_2` interpolate those three points. The exact Newton identity is

\[
\boxed{T_{s,t}(k)=Q^{(2)}_{s,t}(k)+
\sum_{j=k_0+3}^{k}\varepsilon_{s,t}(j)
\binom{k-j+2}{2}.}
\]

For each fixed impulse position `j`, the correction vanishes before `j`, then follows a **quadratic in k**. Therefore, when the defect turns off, the corrected sequence returns to a (different) quadratic. This is NOT a global cubic polynomial or a proof of eventual exact polynomiality.

As a working example take `s=2,t=1`, with phase-corrected PT sequence at `k=1,...,5`:

\[
(1359443250,\ 4419692500,\ 9723323500,\ 17236800750,\ 26960124250).
\]

The first-three-point quadratic predicts `Q_2(5)=27060730750`. The observed pulse at `j=4` is `eps(4)=-33535500`. Predicting that `eps(5)=0` gives

\[
27060730750 +3(-33535500)=\boxed{26960124250},
\]

exactly the held-out degree-27 integer. The original quadratic misses by 100606500; a global cubic perpetuates the wrong correction. The one-pulse assumption is a predictive hypothesis, not an exact identity before observing the next point.

## 3. Dataset and strict holdout

The dataset uses 24 `(s,t)` streams, with 3 labeled points in the degree-6–20 calibration panel, one revealed update in degree 21–25, and one withheld point in degree 26–30. The sign phase is learned from calibration first-shell labels; there is no modularity or wall-crossing input.

Every 5-point sequence was exactly reconstructed via the discrete Green-function identity. This is a consistency check, not a predictor: two pulses are observed ex post and exactly restore two held-out values. The meaningful predictive measurements are the prior-support and next-step forecasts below.

**Observed pulse types:** 11 streams have no pulse at either later step, 4 have one pulse and then heal, and 9 remain active in both steps. No `0 -> nonzero` reactivation occurred in the second step.

| Degree residue `s` | Active depths in first update | Active depths in next update |
|---:|---|---|
| 0 | 4, 5 | 5 |
| 1 | 1, 2, 3, 4 | 1, 2, 3, 4 |
| 2 | 1, 2, 3 | 2, 3 |
| 3 | 2, 3 | 3 |
| 4 | 3, 4 | 4 |

Note residue 0 has `k=5,6` at the two updates, whereas residues 1–4 have `k=4,5`.

Define the defect-generating polynomial

\[
E_s(y,k)=\sum_t\varepsilon_{s,t}(k)y^t.
\]

For the observed windows, `E_s(y,k)` is divisible by `y^{h_s(k)}` with a nondecreasing frontier `h_s`, and its nonzero support is contiguous above the frontier. The lower frontier often moves by one shell per step, but not always: `s=1` remains pinned between the two updates. **This is the candidate inter-shell recursive structure.**

## 4. A second recursive pattern: paired-residue curvature locking

Let `C=Delta^2 T`. The mirror pairs `s=1 <-> 4` and `s=2 <-> 3` exhibit exact equality of low-shell curvatures. In the first and second forward windows:

| Residue pair | Common depths, first update | Common depths, next update |
|---|---|---|
| `1 <-> 4` | `{0}` | `{0,1}` |
| `2 <-> 3` | `{0,1}` | `{0,1,2}` |

For example, at residue `s=2,t=2`, the curvature at the second update becomes exactly equal to the earlier observed curvature in its paired family `s=3,t=2`. The already existing v0.6 mirror-curvature prediction thereby repairs two extra integers relative to pure quadratic extrapolation.

An appealing hierarchy is

\[
C_{s,t}(k)=C^{\mathrm{shared}}_{\{s,5-s\},t}
       +R_{s,t}(k),
\]

where the **residual** `R` becomes zero shell by shell. This is a compression of the nontrivial curvature information; however the first two windows do not suffice to infer a universal, all-degree locking law.

## 5. Held-out model comparisons

Prediction methods were fitted on `Q<=25` and evaluated on degree `26..30`. All integer comparisons are exact; numerical errors are displayed only after the exact check.

| Model | Exact matches | Mean absolute relative error | Maximum relative error |
|---|---:|---:|---:|
| Local quadratic | 15/24 | 0.008758% | 0.199362% |
| Global cubic / constant fourth difference | 11/24 | 0.163263% | 2.83509% |
| `(A k+B) Q_leading(k)+C` | 10/19 eligible | 0.021795% | 0.393075% |
| Existing mirror-curvature rule | **17/24** | **0.0000650%** | **0.000741%** |

For `t>0`, the product model's three rational parameters were fitted exactly on the most recent three observed points. For `t=0`, predicting a shell from its own leading value would leak the target, so these five cases are excluded; **10/19 must not be compared as if it were 10/24**.

The failure of global cubic and product extrapolations is informative: adding algebraic complexity without a sparse correction gate overfits the transient.

## 6. Attenuation and modular arithmetic of the residual

Across 9 persistently active streams, `|eps_new|<=|eps_old|` holds in all 9, while the third-difference sign reverses in 8 and `v_5(eps_new)>=v_5(eps_old)` in all 9. The nine absolute attenuation factors vary enormously; **a single damping parameter cannot fit them**. Along the moving front, the first still-active shell often has a relative attenuation of order `10^-2`, whereas the pinned `s=1` family has much larger residuals in its shallower shells.

An intentionally conservative *empirical set-valued prediction* is

\[
\varepsilon_{\mathrm{next}}\in
\begin{cases}\{0\},&\varepsilon_{\mathrm{previous}}=0,\\[-|\varepsilon_{\mathrm{previous}}|,|\varepsilon_{\mathrm{previous}}|],&\text{otherwise.}
\end{cases}
\]

This covers all 24 held-out defect values (11 cases as exact singletons); it is retrospective, extremely wide at deep shells, and **does not have distribution-free coverage guarantees**. Additional arithmetic priors could restrict it to multiples of `5^{v_5(eps_prev)}`. Their empirical validity likewise needs more forward data.

## 7. Proposed recursive numerical baseline v0.7+

An economical future predictor should have three nested modules:

1. **Front/state:** identify `stable`, `transient`, or `paired-curvature locked` from `(s,t,k)`, existing neighboring-shell states, and historical non-reactivation.
2. **Shared curvature:** predict the next `Delta^2 T` from its residue partner where the locking front is likely to advance.
3. **Residual pulse:** on unresolved sites, predict `eps` in an exact integer lattice using sign, `v_5`, and features of **previous pulses in adjacent shells**. Insert it through the exact binomial Green function instead of fitting large PT coefficients directly.

Train state and integer models separately. Use chronological degree-group folds, penalize model complexity (minimum description length), preserve Python `int`/`Fraction`, and report sign, exactness, signed-log magnitude error, and abstention coverage. A neural network has no clear advantage with only 24 site labels.

## 8. Boundary of the evidence

The published PT reference integers were derived using techniques including modularity, but the **prediction algorithm** uses none of those techniques. The front model and feature family were refined after earlier inspection of Q26–30, so the result is a **retrospective empirical baseline**. Degrees 31–35, saved in `results/forecasts_Q31_Q35.csv`, remain unverified. No iteration through 36–40 is performed. In particular, the exact Newton-pulse identity is true for *all* sequences and does not itself imply future defect sparsity; the latter is the scientific conjecture to test.

## 9. Direct test of a recurrence **in shell depth**

One can go beyond separate shell pulses and ask whether the whole defect polynomial transports recursively,

\[
E_s(y,k+1)\stackrel{?}{=} y^{\ell_s}(A_s+B_sy)E_s(y,k).
\]

The observed front suggests `ell_s=1` for the advancing families and `ell_1=0` for the pinned family. As a conditional composability test, use the first **two actual next-block defect coefficients** to solve for `A_s,B_s`, and test the remaining next-block shells. Only the pinned `s=1` family has enough remaining shells for this to constitute a nontrivial check. Exact arithmetic yields `A=-45/128`, `B=-23675/16384`. The model fails badly at depths 3 and 4, including an incorrect sign at depth 3. These are **not blind holdouts** because the first two next-block labels were used in the fit; the experiment only rejects the overly simple `linear × old defect polynomial` ansatz.

Thus both kinds of global multiplication fail: `T_t(k)` as a linear-times-leading quadratic, and `E_s(y,k+1)` as a scalar/linear multiple of the older defect shell polynomial. A genuinely useful recurrence must allow **multiple transient modes or additive new shell sources**, together with the moving-front support. The code `shell_transfer.py` and its output `results/conditional_shell_transfer.csv` make this negative test reproducible.
