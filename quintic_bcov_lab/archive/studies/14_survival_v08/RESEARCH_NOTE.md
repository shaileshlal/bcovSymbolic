# Frozen recursive-defect rules survive degrees 31–53

*Quintic PT numerical baseline v0.8 — exact-integer, forward-degree evaluation*

## Executive result

A frozen v0.7 empirical PT predictor, trained/calibrated on degree <= 30 and operating on 24 residue/shell streams, matches **22/24** PT integers at degrees 31–35, then **106/109** through degree 53 when previous true blocks are revealed only *after* scoring. With no revelation and model predictions fed back recursively, it matches **99/109** integers. Every reference comparison uses exact Python integers. The largest relative discrepancy in the autonomous run is `3.276669283956296e-05%`; the mean absolute relative error is `9.43616188155617e-07%`.

The comparator is **not an unconditional enumeration**. It is computed from the independently reconstructed quintic rank-zero DT series and an applicable wall-crossing formula; those are not present in the prediction process. The comparator first agrees on all 48 reference entries at degrees 21–30. Four leading coefficients at degrees 35, 40, 45, 50 independently match the extremal geometry theorem of Liu--Ruan.

## Definition of the empirical prediction

For degree `d=5k+s` and shell depth `t`, define the observed sign-corrected PT sequence

\[
T_{s,t}(k)=\sigma_s(k)\operatorname{PT}\bigl(5k+s,m_C(5k+s)+\delta_s+t\bigr),
\quad \varepsilon_{s,t}(k)=\Delta^3 T_{s,t}(k).
\]

Here `m_C(d)=-floor((d^2+5d)/10)`, `(delta_0,...,delta_4)=(0,1,2,2,1)`, and the sign phase was fitted solely to earlier PT integers. The v0.7 **unchanged** predictor forms the local quadratic forecast, and, if the paired residue's second difference lies between the previous two observed curvatures, substitutes that paired curvature. It never uses future PT values.

The exact Newton identity is

\[
T(k)=Q_2(k)+\sum_{j=k_0+3}^k\varepsilon(j)\binom{k-j+2}{2}.
\]

It is an exact reorganization, not by itself a prediction rule.

## Two genuinely different tests

**Rolling update (teacher forced).** For each 5-degree block, predict all shells in the block, score against reference, and *only afterward* reveal those reference values to the next block. This returns 106/109 exact integers. It is suitable for measuring one-step transfer as newly measured PT data arrive.

**Autonomous rollout.** Predict each new block, append only those predictions, and repeat to degree 53. No reference value at degree >30 is supplied anywhere in the rollout. This returns 99/109 exact integers, so it is a genuine recursive-survival test rather than an oracle-assisted bootstrap.

| Degree block | Rolling | Autonomous | Largest relative error in autonomous block (%) |
|---|---:|---:|---:|
| 31–35 | 22/24 | 22/24 | 0.0000065003 |
| 36–40 | 23/24 | 22/24 | 0.0000139342 |
| 41–45 | 24/24 | 22/24 | 0.0000209534 |
| 46–50 | 24/24 | 22/24 | 0.0000272361 |
| 51–53 | 13/13 | 11/13 | 0.0000327667 |

The model's degree-31–35 outputs match its **originally archived v0.7 forecast** entry by entry, demonstrating the tests did not change the starting prediction after inspecting labels.

## The defect front survives

With true PT reference sequences (for *evaluation only*), the counts of nonzero third differences are

\[
9\ \xrightarrow{26-30}\ 4\ \xrightarrow{31-35}\ 2\ \xrightarrow{36-40}\ 1\ \xrightarrow{41-45}\ 0\ \xrightarrow{46-50}\ 0\ \xrightarrow{51-53}.
\]

| Degree block | Exactly nonzero defect sites `(residue, depth)` |
|---|---|
| 26–30 | (0,5), (1,1..4), (2,2..3), (3,3), (4,4) |
| 31–35 | (1,2), (1,3), (1,4), (2,3) |
| 36–40 | (1,3), (1,4) |
| 41–45 | (1,4) |
| 46–50, 51–53 | None |

Thus **all** new nonzero defects are supported on the *previously active front*, and no stabilized shell reactivates. The final defect is at degree 41 in residue 1, shell 4. All 24 streams have zero third difference over their subsequent observed positions once stabilized.

A simple *post-hoc* frontier formula summarizes the observed last nonzero third differences. For residue `s=1,2,3,4`, across the active shells the last defect occurs at `k_last=t+5-s`; for residue zero's active shells 4,5, `k_last=t+1`. The previous samples already hinted at this diagonal pattern, but the exact formula was tabulated after inspecting the extension; it must not be called a prospective all-degree law.

## Curvature locking also survives

Mirrored residue pairs share second finite differences at more shell depths as k increases:

| k | Pair (1,4), locked depths | Pair (2,3), locked depths |
|---:|---|---|
| 4 | 0 | 0,1 |
| 5 | 0,1 | 0,1,2 |
| 6 | 0,1,2 | 0,1,2,3 |
| 7 | 0,1,2,3 | 0,1,2,3 |
| 8 and 9 | 0,1,2,3,4 | 0,1,2,3 |

This is a strong cross-shell, cross-residue *empirical* description that survives beyond the degree-30 calibration interval. It strongly favors shared-curvature and short-pulse models over fitting a permanent cubic.

## Where the naive models fail

| Model | Exact out of eligible | Notes |
|---|---:|---|
| Local quadratic | 102/109 | Ignores defect pulses |
| Global cubic | 93/109 | Persists transients that die out |
| Linear times leading quadratic | 79/86 | Leading shell excluded to avoid target leakage |
| Frozen mirror-curvature | **106/109** | No parameter changes since v0.7 |

All three rolling-update errors occur on the deepest residue-one shells: (31,3), (31,4), (36,4). In the autonomous run, these two early mistakes induce errors only in `(s=1,t=3)` and `(s=1,t=4)` at later degrees. The former accumulates with a linear absolute envelope and the latter with a quadratic envelope; neither spreads into the remaining 22 streams.

## Reference provenance and independence

The *prediction* code imports only `recursive_defects`, standard-library arithmetic, and observed PT labels through degree 30. Its checksum matches v0.7 exactly.

The *validation* program first reconstructs rank-zero DT series from nine exact low-order polar/vanishing conditions and then applies the AFKPS wall-crossing map at 109 sites. All sites satisfy that formula's stated applicability inequality. Only one has even a candidate extra wall `(degree=31,depth=4)`, but that term multiplies `PT(1,0)=0`, so every target equals its vacuum term. Those reference integers are conditional on the modular DT reconstruction and wall-crossing setting: comparing to them is **not** a test independent of the underlying modular-bootstrap framework.

The extremal geomerty cross-check uses the independent Liu--Ruan formula at degrees `d=5k`:

\[
n_{g_C(5k),5k} =(-1)^{\binom{k+3}{3}-\binom{k-2}{3}+3}\,5\left(\binom{k+3}{3}-\binom{k-2}{3}\right),
\]

with leading PT equal to extremal GV. It gives `(d,PT)=(35,-550),(40,725),(45,925),(50,-1150)`, all matching the reference.

## Scientific interpretation

This is now a useful **frozen empirical baseline**: it predicts not only large PT coefficients but also the support and disappearance of their residual corrections across degree blocks, using a parsimonious recurrence. However, the tested shells are extremely close to the extremal boundary; the reference itself becomes exactly quadratic in degree within each residue and shell once the non-vacuum walls are absent. We have **not** tested deeper shells with competing wall contributions, or the generic intermediate holomorphic ambiguity at arbitrary genus. The next demanding test is to increase shell depth to 6--10 rather than simply push the same shallow strip to larger degrees.

## References

- Alexandrov, Feyzbakhsh, Klemm, Pioline, Schimannek, *Quantum geometry, stability and modularity*, arXiv:2301.08066.
- Liu and Ruan, *Castelnuovo bound and higher genus Gromov--Witten invariants of quintic 3-folds*, arXiv:2210.13411.
