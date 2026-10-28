# Near-Castelnuovo PT/DT arithmetic without modularity as a predictive input

**Executable empirical study — 28 October 2026**

## A. Question and restriction of information

Can we predict higher-genus quintic holomorphic-ambiguity lattice integers using *only* patterns in previously observed integers? No modular ansatz, representation-theoretic decomposition, DT/PT wall crossing, or BCOV high-genus amplitude is inserted into the predictors here. For evaluation, known published integer tables are used **after fitting**, and their modularity-conditional provenance is explicit.

The evidence currently consists of one CY threefold, degree labels 1–30 and at most six rows near the genus bound. This is a discovery experiment, not a proof.

## B. Coordinate system for near-boundary data

Let

\[
g_C(d)=\left\lfloor\frac{d^2+5d+10}{10}\right\rfloor,
\qquad m_C(d)=1-g_C(d),\qquad d=5k+s.
\]

The first potentially nonzero near-extremal PT/GV coefficient occurs at offset

| residue `s` | 0 | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|---:|
| `delta_s` | 0 | 1 | 2 | 2 | 1 |

Define the shell coordinate `t ≥ 0` by

\[
P_{s,t}(k)=\mathrm{PT}(5k+s,\;m_C(5k+s)+\delta_s+t).
\]

Our input has training labels at degrees 1–20 and **48** held-out PT numbers at degrees 21–30. Exactly the same number of GV reference integers are held out. The complete set of 48 PT coefficients consists of 10 each at shell `t=0,1,2,3`, six at `t=4` and two at `t=5`. Whole-degree splits avoid leaking a high-degree near-boundary label into an adjacent training shell.

## C. Sign periodicity: a robust and cheap feature

The lowest nonzero PT signs at degrees 1–20 determine the *smallest common* periodic phase `p=4` in `k` across the five residues. In the learned convention,

\[
\sigma_s(k)=\begin{cases}
+1,&k\bmod4\in\{0,1\},\ s=0,2,4,\\
+1,&k\bmod4\in\{0,3\},\ s=1,3,\\
-1,&\text{otherwise}.
\end{cases}
\]

Train a constant orientation for each shell `t` after multiplying by this phase. On degrees 21–30 the PT sign rule is **48/48** correct. The corresponding GV rule is **43/48** correct, and begins changing signs unexpectedly at several interior sites. The PT basis is therefore more convenient empirically.

## D. Finite differences: direct quadratic vs mirror transfer

Let `T_{s,t}(k) = sigma_s(k) P_{s,t}(k)`. Fit the unique quadratic passing through the three samples at degrees 6–20 in each residue/shell; predict two later values at degrees 21–30 without updates. The PT fit is **22/48** exact; the corresponding GV direct quadratic gives **10/48** exact.

Several PT shell sequences are exactly quadratic across all five available points. For instance `s=0`, `t=0..3`; `s=4`, `t=0..2`; and `s=3`, `t=0..1` all satisfy

\[
\Delta_k^3 T_{s,t}(k)=0
\]

through the available data, where `Delta` is the forward finite difference. Higher shells frequently exhibit small early-degree defects rather than exact polynomial continuation.

Moreover, mirrored residues `s ↔ 5-s` share affine relationships between the first shell and later shells. Fitted on Q≤20:

\[
\begin{aligned}
T_{4,1}&=132T_{4,0}-1\,138\,500,\\
T_{4,2}&=437968T_{4,0}-7\,554\,948\,000,\\
T_{3,1}&=\frac{1767877}{4874}T_{3,0}-441\,969\,250.
\end{aligned}
\]

Transfer these relations to `s=1` (first two) and `s=2` (third) **only when the relation is exact throughout the source family's three training points and the target's last training relative error is ≤10^-4**. The rule thus makes only six mirror-based test predictions, all six correct. The combined predictor gives **28/48 exact PT coefficients**, 48/48 correct signs, and median relative error **0.01809% on the 20 incorrect integers**. It attains all 10/10 exact values on each of the first two shells, six on the third, two on the fourth, and none on the two deepest shell bands.

Notably, at degree 21 the three predicted PT values are:

\[
\mathrm{PT}(21,-53)=319125,\quad
\mathrm{PT}(21,-52)=40986000,\quad
\mathrm{PT}(21,-51)=132211590000.
\]

The first equals the previously discussed `n_{54,21}`. All three match the published held-out integers. The next predicted PT coefficient, at `m=-50`, is *not* exact; relative error approximately 0.00478%, giving a concrete limit to these rules. The worst held-out relative error is substantially larger, **8.50%** at degree 26 on shell 4. The very small median miss error should not obscure this depth-dependent tail risk.

## E. Learn PT → GV triangularity from integers, not a generating-function identity

Training pairs reveal that `P_0 = GV_0` and the first neighbor coefficient obeys the exceptionally simple *empirical* relation

\[
P_1-G_1=(2h-2)G_0,
\qquad h=g_C(d)-\delta_s,
\]

at 14/15 training degrees Q6–20 (exception Q=6). A maximum-consensus rational interpolation algorithm, given only the numerical data, discovers this affine factor. Subtract it and repeat: a second empirical correction is

\[
P_2-G_2=(2h-4)G_1+\binom{2h-2}{2}G_0,
\]

at 12/15 calibration degrees (exceptional Q=6,7,11). No analytic PT/GV formula was input to these fits. They are, of course, consistent with the known PT/GV generating function, which can be used for explanation *after* fitting.

On the independent degree-21–30 split, these rules and the **predicted**, rather than observed, PT coefficients recover **26/30 GV integers** over the first three shells: all ten at t=0, all ten at t=1, six at t=2. PT and GV scores are not independent pieces of evidence.

## F. A DT-like sequence emerges from finite differences

A particularly useful compression is the curvature

\[
C_{s,t}(k)=T_{s,t}(k+2)-2T_{s,t}(k+1)+T_{s,t}(k).
\]

For the three representative residues `s=0,4,3`, the curvature, after dividing by 5 and an alternating sign, yields the following integer streams *using only degrees ≤20*:

| Representative residue | Extracted sequence (successive shells) |
|---|---|
| s=0 | `5, -800, 58500, 5817125` |
| s=4 | `8625, -1138500, 3777474000` |
| s=3 | `-1218500, 441969250` |

These **nine** integers coincide exactly with the low-order rank-zero DT coefficients in equation (5.8) of the source paper, *after fitting*. The predictors did not use the DT coefficients. The comparison is not an independent confirmation of the original reference tables, which are related to DT data themselves.

For the next coefficient of each stream, using only the earliest PT curvature yields an **approximate** forecast:

| DT series slot | Train-PT estimate | Published reference | Relative error |
|---|---:|---:|---:|
| h0, next | 75,503,815,750 | 75,474,060,100 | +0.039425% |
| h1, next | 3,099,947,260,875 | 3,102,750,380,125 | −0.090343% |
| h2, next | 951,927,805,125 | 953,712,511,250 | −0.187133% |

After incorporating the withheld degree-21–30 PT data, the *last two consecutive* second differences become equal, and these three curvature estimates become exact. Thus the deeper PT shell data behave like an arithmetic channel carrying three DT coefficient streams.

The experiment does **not** establish a recurrence generating the DT coefficients along their own q-series index. The known streams grow irregularly in the first few orders and straightforward ratio extrapolation cannot predict the next coefficient exactly. The predictive operation demonstrated here is **PT degree propagation at a fixed shell**, not arbitrary-depth DT sequence generation.

## G. Conventional ML ablations

Training ordinary regressors on a one-hot-residue / degree / shell feature matrix and the log magnitude of PT coefficients yields considerably weaker held-out predictions:

| Model | Median held-out log10-magnitude error | Median multiplicative factor |
|---|---:|---:|
| Random forest | 0.354 dex | 2.26× |
| ExtraTrees | 0.267 dex | 1.85× |
| Five nearest neighbours | 2.159 dex | 144× |
| Polynomial ridge | 0.989 dex | 9.74× |
| Exact signed periodic/quadratic + mirror | 0 dex overall | 1.00× |

The exact arithmetic model has 28/48 exact matches; its median error is zero because more than half its predictions are exact. The miss-only relative error above is a better assessment of the unpredicted region. Neural networks/transformers are not used here: the integer dataset is too small to make their superiority plausible.

## H. Continuation beyond degree 30

After adding degrees 21–30 as *new calibration observations*, I looked for rows with two consecutive equal second differences, or recent exact mirrored-affine relations. Exactly **38 of 48** degree-31–40 near-extremal sites meet this recurrence-activity filter. These are *candidate forecasts, not independently verified BPS invariants*. Sites without such a pattern are blank rather than silently assumed zero. Selected candidates:

| Degree | Predicted leading PT/GV | Predicted next PT | Predicted next-next PT |
|---:|---:|---:|---:|
| 31 | −724,500 | −94,495,500 | −309,752,868,000 |
| 35 | −550 | −87,200 | −6,318,000 |
| 40 | +725 | +115,200 | +8,365,500 |

The degree-35 and degree-40 **leading terms** independently agree with the established extremal-curve formula; the subleading terms have *not* been independently verified. This is an important distinction.

## I. Interpretation, failures and safeguards

1. A **degree-period-four phase** plus **three effective DT-like integer streams** is a compact compression of many of the near-boundary PT coefficients.
2. The finite-difference polynomial law is cleanest in PT, not GV, and deteriorates with increasing shell depth. A model of the defects, perhaps via sparse change points and local neighborhood attention, is the next empirical problem.
3. The next uncovered coefficient at each deeper shell should be treated as a new integer target; empirical recurrence does not yet imply that finitely many initial DT coefficients determine all shells without additional analytic input.
4. The model family and mirror-pair idea were chosen in the context of earlier project investigations. The train/test numerical procedure is isolated, but the overall hypothesis discovery was **exploratory and retrospective**, not a registered blind competition.
5. The Q>20 published validation values are **modularity-conditional**. This paper-source dependence enters the *reference labels*, not the predictor. To establish genuine independent physical predictions, compare future candidates against an independent BCOV amplitude reconstruction or geometric curve-counting calculation.
6. The experiment concerns a near-Castelnuovo strip at degrees 1–30, not the full intermediate ambiguity lattice at arbitrary genus.

**Bottom line:** direct arithmetic discovery is genuinely productive, and the predictive structure is clearer in the finite differences of PT than in raw GV/DT q-series coefficients. But the successes currently concentrate in the first 2–3 shells; no all-genus completeness claim is justified.

### Reference table provenance

S. Alexandrov, S. Feyzbakhsh, A. Klemm, B. Pioline, T. Schimannek, *Quantum geometry, stability and modularity*, arXiv:2301.08066v2 (2024), Tables 2–3 (GV/PT) and eq. (5.8) (DT q-series): https://arxiv.org/html/2301.08066. The arithmetic formulas are **our empirical fits**, not the paper's modular reconstruction.
