# The next quintic PT residual: a recurrence of sparse finite-difference pulses

**Exact-rational empirical study v1.2 — 28 October 2026**

## Scope, provenance, and the experimental distinction

This is a **purely empirical learner** of near-extremal PT integer sequences. It uses the previously learned period-four signs and quadratic envelopes; it does not use modular forms, wall-crossing equations, DT-series coefficients, or GV invariants in any prediction routine. The conditional exact PT **reference** is generated separately using the quintic rank-zero DT series and the PT wall-crossing implementation inherited from v1.1. Thus these comparisons evaluate prediction of that conditional reference, not independent mathematical truth of its assumptions.

The mature normalized scales used to extract the signals are computed **from observed mature PT integers**: one degree-labelled PT datum per residue-pair group and depth. This uses only integer arithmetic, not a DT formula, but those mature PT observations can occur at *larger curve degree than the target*. Consequently all results involving such a scale are conditional residual/shape benchmarks; they are **not blind, causal forecasts of complete PT coefficients**.

The representation and pulse-support feature engineering are retrospective across this research program. The exact forward-index tests reported below protect the fitted numerical coefficients from their holdout labels, but are not preregistered confirmations of a previously unknown architecture.

## Definitions

For degree `d=5k+s`, with residue `s∈{0,...,4}`, shell depth `t`, the earlier calibrated envelope is

\[
W_{s,t}(k)=5k^2+(2s-5)k+c_s-2t,\qquad (c_0,\dots,c_4)=(10,6,4,4,6).
\]

The known four-periodic sign is `σ_s(k)∈{−1,1}`. Put

\[
A_{s,t}(k)=\frac{\sigma_s(k)\,\operatorname{PT}(5k+s,t)}{W_{s,t}(k)},\quad
K_s(t)=\begin{cases}t-1,&s=0,\\t+3-s,&s>0,\end{cases}\quad
r=K_s(t)-k.
\]

The three mature-scale groups are `α(0)=0`, `α(1)=α(4)=1`, `α(2)=α(3)=2`. A scale `B_{α,t}` is obtained from an observed mature PT integer (`r=0`) in the same group at depth `t`: it is the ratio of PT to `σW` at that anchor. **All 628 normalized samples** generated this way agree with a separate direct scale evaluation over the rationals.

Define the *integer* (in the evaluated set) residual and its third finite difference

\[
Z_{s,r}(t)=(-1)^t\big[\sigma_s(k)\operatorname{PT}(5k+s,t)-W_{s,t}(k)B_{\alpha(s),t}\big],
\quad \varepsilon_{s,r}(t)=\Delta_t^3 Z_{s,r}(t).
\]

The original affine-in-`1/W` model is equivalent to *quadratic behavior in `t`* for `Z` at fixed `(s,r)` (more specifically, affine in `W`, itself quadratic in `t`). Therefore its first failure is exactly `ε≠0`.

## Result A. The new obstruction is a finite transient, not a permanent cubic

The `r=6` case has isolated nonzero third-difference pulses. In the observed streams, the nonzero pulses end at `t=15` for residues `0,2,3,4`; residue 1 has an earlier pulse at `t=11`. Thereafter `Δ_t^3 Z=0`, and a local two-point affine-`W` fit suffices. In a forward-*degree* holdout (`d>53`, with training `d<=53`), the two latest calibrated points predict **32/32** normalized `Z` values exactly; a three-point quadratic in `t` also gets **32/32**. A two-point fit from the very earliest calibration values gets **0/32**, because it crosses the early transient. This contrast depends on the chosen calibration window; it is not a claim that polynomial degree is the important modelling feature. What matters is *which side of the pulse front* is used to fit.

This shape score presumes the correct depth scales `B_{α,t}` are observed. A user cannot convert all 32 conditional `Z` predictions to independent PT predictions without those scales.

## Result B. Terminal front and quadratic recurrence in front-distance

For the tested streams, the final nonzero third-difference pulse is at

\[
T_s(r)=\begin{cases}2r-1,&s=1,\\ 2r+3,&s=0,2,3,4.\end{cases}
\]

Write `τ_s(r)=ε_{s,r}(T_s(r))`. Exact integer interpolation from three adjacent `r` values gives:

\[
\begin{aligned}
\tau_0(r)&=-371185562500(r-2)(5r-9),\\
\tau_1(r)&=2627390625(5r^2-27r+38),\\
\tau_2(r)&=-1523125(5r^2+5r+6),\\
\tau_3(r)&=-2627390625(5r^2-3r+2),\\
\tau_4(r)&=-371185562500(r-1)(5r-6).
\end{aligned}
\]

The shared curvatures `Δ_r²τ_0=Δ_r²τ_4=−3,711,855,625,000` and `Δ_r²τ_1=−Δ_r²τ_3=26,273,906,250` reduce the necessary numerical calibration to **13 integer pulses** rather than 15; residue 2 has its own curvature `−15,231,250`. With these 13 labels fixed, **18/18 subsequent terminal-pulse integers match exactly**. The 18 evaluated sites use `r=9..11` for residue 0, `9..12` for residue 1, `7..9` for residue 2, `6..9` for residue 3, and `7..10` for residue 4.

At each terminal site, the independent validation observation has `ε≠0` and the following `ε` is zero. At a further later-index check, the proposed terminal-front rule matches **13/13** observable terminal activity states and **33/33** observable post-terminal zeros (indexed `r>=9`). These front tests describe the finite dataset; the form of the front was developed with exploratory access to it.

Given the preceding three observed normalized values, the predicted terminal pulse predicts the **next normalized coefficient** by

\[
Z(T)=3Z(T-1)-3Z(T-2)+Z(T-3)+\tau_s(r).
\]

All **18/18** such conditional next-`Z` reconstructions are exact, with the target `Z(T)` withheld.

## Result C. The defect of the defect

At offset `ℓ≥0` behind the terminal front define

\[
E_{s,\ell}(r)=\varepsilon_{s,r}(T_s(r)-\ell).
\]

A fixed forward-index test predicts `r=10,11` in each residue/offset stream having at least three nonzero earlier observations (with training `r<=9` and `ℓ<=3`). Using the most recent three prior observations per stream, **17/19** evaluable predictions match exactly; **five** of the 24 prespecified sites lack complete conditional reference observations and are not scored. The two errors are at `(s,ℓ,r)=(3,2,10),(3,2,11)`, with `prediction - truth = −158,405,000` and `−475,215,000` respectively.

Their errors are explained by a *second-level pulse*:

\[
\Delta_r^3 E_{3,2}(10)=158405000,\qquad
\Delta_r^3 E_{3,2}(11)=0.
\]

An earlier neighboring layer `(s,ℓ)=(3,1)` has its pulse of `97,480,000` two steps earlier, giving a retrospective amplitude ratio `13/8`. This is an unvalidated possible shift/transport relation; it was **not** used to improve the benchmark.

## No finite termination claim

The experiments suggest a *recursive correction grammar*:

\[
Z\overset{\Delta_t^3}{\longrightarrow}\varepsilon
\overset{\text{terminal alignment}}{\longrightarrow}\tau_s(r)
\overset{\Delta_r^3}{\longrightarrow}0
\]

for the terminal layer, with further sparse `Δ_r³` corrections at subterminal layers. This **does not** establish all-depth finite generation. New pulse streams can appear at larger offsets, and the mature scale streams `B_{α,t}` still require independent integer inputs or an as-yet unknown recurrence. No DT modular formulas or wall-crossing relations were supplied to the empirical predictor.

## Reproducibility

`python run.py` requires only Python standard library and `data/`. Then run `python validation/score.py` to compare with the separate exact conditional PT observations. `python tests.py` reruns the predictor in an isolated temporary directory with **no validator or reference files**, verifies byte-identical outputs and reconstructs the normalization from observed mature PT integers. Optional reference regeneration uses SymPy and `validation/reference_backend/` through `python validation/generate_observations.py` and `python validation/normalize_observations.py`. A marimo notebook is included for interactive exploration; its GUI was not run in the generation environment.

## Prospective terminal-pulse integers (unverified)

The fitted terminal-pulse recurrences also give the following targets beyond the independently scored and observable terminal positions. These predictions require a future conditional PT reference of sufficient order to test; they are **not** full PT/GV predictions:

| Residue `s` | Front distance `r` | Terminal depth `T_s(r)` | Predicted `τ_s(r)` |
|---:|---:|---:|---:|
| 0 | 12 | 27 | −189304636875000 |
| 1 | 13 | 25 | +1397771812500 |
| 2 | 11 | 25 | −1014401250 |
| 3 | 11 | 25 | −1508122218750 |
| 4 | 11 | 25 | −181880925625000 |

These are saved in `results/terminal_pulse_prospective.csv` and should be checked using independent PT observations when available. Even a complete exact match on this list would not establish a finite recursion for all higher front distances, since new off-terminal pulse profiles and the mature-scale sequences remain unresolved.
