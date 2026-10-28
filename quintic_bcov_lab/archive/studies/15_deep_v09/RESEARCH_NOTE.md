# Does the empirical PT recurrence survive at deeper quintic shells?

*Numerical experiment v0.9, 28 October 2026. Depths 6–10, curve degrees 31–53.*

## Question and experimental boundaries

We previously obtained excellent forward predictions at shallow near-extremal PT shells by removing a period-four sign and modelling sequences across degrees using quadratic finite differences and sparse corrections. This study increases the **shell depth** to `t=6,...,10`. It does not increase genus from 53, and makes no claim to a full BCOV free-energy computation.

To respect the empirical-only requirement, the prediction algorithm reads only *PT integers*, shell/degree positions and exact rational arithmetic. It does not use modularity, BPS walls, or a DT formula. A completely separate conditional reference module supplies PT integers for scoring after prediction.

The reference is based on the rank-zero DT modular ansatz with nine rational basis coefficients fixed by seven polar numbers plus forced zeros, and the PT wall-crossing formula. Its depth extension requires rank-zero q-order 13 (rather than the previous 12), so the module is truncated at order 14. For any wall contribution with lower PT degree >=2 whose required value was not independently certified, **the site is withheld from scoring**. This is necessary to prevent fabricated zero contributions.

All 109 existing shallow PT reference values at degrees 31–53 were reconstructed exactly by the extended reference solver; 106 of 115 new candidate sites at `t=6,...,10`, `Q=31,...,53` are covered. The nine incomplete cases occur before Q=41; all 65 deep sites at Q=41,...,53 are valid within the conditional reference.

## A universal quadratic shape, learned from shallow data

Write degree `Q=5k+s`, `s=0,...,4`, and define the PT shell label `m=m_C(Q)+delta_s+t`, where `m_C(Q)=-floor((Q²+5Q)/10)` and `delta_s=(0,1,2,2,1)`.

The previously learned leading-sign phase `sigma_s(k)` has period four in k. Fit the sign-corrected *quadratic* polynomials for shallow depths t=0 and t=1, using only Q<=30. Their ratios of polynomial coefficients yield exactly:

\[
W_{s,t}(k)=5k^2+(2s-5)k+c_s-2t,
\qquad (c_0,c_1,c_2,c_3,c_4)=(10,6,4,4,6).
\]

The extension in t is an **empirical extrapolation** of the observed `-2t` constant term. It is not inserted from wall-crossing. Normalize

\[
 A_{s,t}(k)=\frac{\sigma_s(k)\operatorname{PT}(5k+s,t)}{W_{s,t}(k)}.
\]

## Frozen stabilization front

The prior shallow-shell experiment suggested a terminal third-difference location

\[
 k_{\mathrm{last}}^{(3)}(s,t)=\begin{cases}t+1&s=0\\t+5-s&s=1,2,3,4.\end{cases}
\]

The appropriate quotient-defect boundary is shifted by two steps:

\[
K(s,t)=\begin{cases}t-1&s=0\\t+3-s&s=1,2,3,4.\end{cases}
\]

The empirical law predicts that `A(s,t,k)-A(s,t,k-1)` is nonzero precisely for observed transitions k<=K, and zero for k>K. On the new deep-shell reference it classifies **91 of 91 transitions**, of which 47 are nonzero. Restricting to the forward-degree holdout Q>=41 gives **64/64** correct. Its third-difference counterpart classifies **41/41** (28 nonzero), including **39/39** on Q>=41.

This is stronger than just fitting a polynomial: after stabilization a **first-order** rational recurrence suffices,

\[
 \mathrm{PT}(5(k+1)+s,t)=
 \frac{\sigma_s(k+1)}{\sigma_s(k)}
 \frac{W_{s,t}(k+1)}{W_{s,t}(k)}\mathrm{PT}(5k+s,t).
\]

At these deep shells there are 66 observed post-front normalized scales, spanning 22 residue/shell streams. All 18 streams with at least two mature observations have identical scales. The scale agrees exactly across paired residues `1<->4` and `2<->3` in seven mature paired-shell comparisons.

## Blind integer prediction, not just relative errors

Calibration for the new depths is limited to labeled rows at Q<=40. The target universe comprises 65 PT integers for Q=41,...,53 and t=6,...,10. The no-oracle predictor's inputs are physically separated from the validator and target values.

| Predictor | Exact integer matches | Coverage | Meaning |
|---|---:|---:|---|
| One-anchor plus paired repeated-scale heuristic, without abstention | 28 | 65 | Every target is predicted; errors may be minuscule in relative terms but nonzero integers |
| Front-gated recurrence using a mature, previously observed same- or paired-residue scale | **32** | **32** | All 32 are exact, 33 deliberate abstentions |

Among the gated targets, 19 lie ahead of the learned front and 14 lack an independently observed mature scale at calibration cutoff; these 33 receive **no claimed prediction**. Thus 100% correctness must not be read as 100% coverage. Neither model has solved the remaining active-shell pulse magnitudes.

## Defect of the defect

The normalized quotient defect,

\[
 D_{s,t}(k)=A_{s,t}(k)-A_{s,t}(k-1),
\]

is much more localized than the original PT integers, and it is **exactly zero after the diagonal front** in every tested transition. For the 24 observed pairs where consecutive D are nonzero, every successor has lower absolute magnitude, and 19 of 24 reverse sign. However, the 5-adic valuation monotonicity previously noticed at shallower depth fails in some deeper sequences. We do **not** promote it to a universal constraint.

The more challenging next task is therefore not to fit higher polynomial degrees blindly. It is to predict the finite set of nonzero quotient pulses before stabilization, perhaps using residue-pair transfer plus a sparse arithmetic state model.

## Claims and limitations

- The PT and DT numerical reference is **conditional on a specific modular/wall-crossing construction**, not an independent unconditional mathematical calculation; its assumptions enter only the scorer.
- The quadratic envelope is inferred from shallower data, and the diagonal-front hypothesis was refined in a previous retrospective experiment. A new, deeper-shell match is a held-out test in the depth dimension, not a prespecified first-principles theorem.
- The exact-prediction rule *uses a mature calibration integer* at the target depth or its paired residue. It cannot determine a new deep-shell amplitude from shallow shells alone.
- Not all deep low-degree sites can be validated with the current lower PT inputs. Nine sites at Q31–40 are explicitly flagged as incomplete.
- The important unsolved target is the active-state integer displacement `D`, particularly near the moving front at t>=9.

## Reproduce

See `README.md`; run predictor scripts *before* the validator. The `tests.py` script checks the strict 109/109 shallow reference replay, 91/91 quotient-front and 41/41 cubic-front identities, 32/32 no-oracle gated forecasts, and a clean predictor-only copy.
