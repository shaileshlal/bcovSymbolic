# Recursive arithmetic in quintic PT scale pulses

*Experimental note v1.0 · 28 October 2026*

## Abstract

We treat the near-extremal quintic PT integers as an exact-rational sequence-learning problem and test recursive structure in the scale-change pulses, not the PT coefficients themselves. The input phase and quadratic envelope are inherited from the previously frozen empirical model v0.9. The primary pulse is transported diagonally from shell to shell with alternating sign at the last active position. Its off-front failure is described by a secondary residual. A stronger, recursively composable description emerges after subtracting the paired-residue stabilized scale: at fixed distance from the diagonal front, the normalized residual is *affine in the inverse quadratic envelope*. Its two coefficients can be learned from earlier data and transported between residue classes. On the 65 depth-6–10 sites at degrees 41–53, the new solver produces 59 PT integers, all matching the separately supplied conditional reference exactly, and abstains on six. This increases exact prediction coverage from the frozen v0.9 result of 32/65 to 59/65. The feature class was discovered retrospectively and requires a genuinely new held-out degree/depth region for a statistical replication.

## 1. Fixed representation and previous baseline

Let `d=5k+s`, `s in {0,...,4}`, with shell depth `t`, and use the degree-dependent sign `sigma_s(k)` learned from the first-shell PT data. The quintic's normalized polynomial envelope found earlier was

\[
W_{s,t}(k)=5k^2+(2s-5)k+c_s-2t,\quad (c_0,c_1,c_2,c_3,c_4)=(10,6,4,4,6).
\]

Define

\[
A_{s,t}(k)=\frac{\sigma_s(k)\,\operatorname{PT}(5k+s,t)}{W_{s,t}(k)}.
\]

The v0.9 empirical stabilization threshold is `K(s,t)=t-1` for `s=0` and `K(s,t)=t+3-s` otherwise. At `k >= K`, the scale is constant, `B_{[s],t}`. The phase normalization is implemented in exact rational arithmetic; because `sigma=+/-1`, division and multiplication conventions coincide.

The 65-site held-out universe is `(d,t)` with `d=41,...,53` and `t=6,...,10`. The previous frozen v0.9 model predicted 32 site labels exactly, abstaining on 33.

## 2. Diagonal transport and its exception

Put `D_{s,t}(k)=A_{s,t}(k)-A_{s,t}(k-1)` and `J_{s,t}(k)=D_{s,t}(k)+D_{s,t+1}(k+1)`.

At the terminal active pulse (`k=K(s,t)`), `J` vanishes. The earlier calibration has 4/4 such cancellations; the enlarged reference contains 15/15, with 0/17 off-front nonzero pulse pairs cancelling. For example,

\[
D_{0,7}(6)=\frac{24796875}{2},\quad
D_{0,8}(7)=-\frac{24796875}{2},\quad
D_{0,9}(8)=\frac{24796875}{2}.
\]

Using this exact empirical terminal rule as a no-oracle update gives 41/41 exact predictions and 24 abstentions, a nine-site improvement over v0.9. These tests use the conditional rank-zero DT/wall-crossing reference **only for scoring**.

Away from the terminal front, `J` is nonzero but generally much smaller than `D` (the last phase of transient corrections). A naive universal ratio `J_next = -[k/(k+4)] J` has about 0.87% mean relative error over eight conditional, true-history comparisons at degrees >=46, versus 100% for predicting zero. This is only a teacher-forced diagnostic, and it adds no exact autonomous predictions. We therefore sought a better exact representation.

## 3. The affine profile and its calibrated coefficients

Set `r=K(s,t)-k`. The observed high-depth residuals satisfy

\[
A_{s,t}(k)-B_{[s],t}=(-1)^t\left(P_{s,r}+\frac{Q_{s,r}}{W_{s,t}(k)}\right).
\tag{1}
\]

This is equivalent to a quadratic envelope multiplied by a scale shift, plus one constant. Thus the `linear × quadratic` idea becomes an *exact integer ansatz* after the sign and scale changes. At each `(s,r)`, only two rational parameters are needed.

The following representative profiles are learned from the calibration data (not from target values):

| Residue s | Distance r | P | Q | How obtained |
|---:|---:|---:|---:|---|
| 0 | 1 | 24796875/2 | 0 | exact constant repeated 3 times |
| 1 | 1 | 14375/2 | 0 | exact constant repeated 5 times |
| 1 | 2 | -1135625 | 28750 | affine fit with 5 agreeing earlier observations |
| 1 | 3 | 163630625/2 | -4513750 | affine fit with 5 agreeing earlier observations |
| 1 | 4 | 17053766875/2 | 322747500 | two-point extrapolation, lower-confidence |
| 2 | 1 | -24796875/2 | 0 | repeated 3 times |
| 2 | 2 | 1611796875 | -49593750 | two late observations, paired validation |
| 3 | 1 | 1751593750 | 0 | repeated twice |
| 3 | 2 | -631827609375 | 7006375000 | joint two-point anti-pair fit |
| 4 | 1 | -1751593750 | 0 | repeated twice |
| 4 | 2 | 631827609375 | -7006375000 | same anti-pair fit |

The `r=2` profile at residue 0 is `(-1611796875,+49593750)`, the negative of the residue-2 profile. The same sign antisymmetry, already visible at `r=1` for the residue pairs `0<->2` and `3<->4`, is then postulated at `r=2` for the second pair. One independent earlier observation from residue 3 and another from residue 4 identify the shared `P,Q` without consulting the held-out degree 41–53 reference.

This is **feature engineering followed by exact constrained fitting**, not evidence for a proven involutive symmetry of the underlying quintic. The shape was discovered during exploratory comparison with the enlarged data; independence must be tested at new depths.

### An explicit prediction

At `d=41,t=8`, `(s,k,r)=(1,8,2)`, `W=286`, and the mature scale for the residue pair `{1,4}` is

\[
B_{14,8}=1867777794427065287875.
\]

Therefore

\[
\operatorname{PT}(41,8)=286\left(B_{14,8}-1135625+\frac{28750}{286}\right)
=534184449206140347572250,
\]

exactly matching the reference.

## 4. Closed recurrence for the secondary residual

The fitted expression (1) implies an exact cancellation of the terms `P` in the diagonal residual:

\[
\boxed{J_{s,t}(k)=(-1)^t\bigg[
Q_{s,r}\left(\frac1{W_{s,t}(k)}-\frac1{W_{s,t+1}(k+1)}\right)
+Q_{s,r+1}\left(\frac1{W_{s,t+1}(k)}-\frac1{W_{s,t}(k-1)}\right)
\bigg],\quad r=K(s,t)-k.}
\]

The formula reproduces **14/14** available nonzero secondary defects at held-out degrees for which both required `Q` values have been fitted. It is an algebraic consequence of (1), hence not independent evidence beyond the reconstructed coefficients. It does explain why the defect of a defect is controlled by a still smaller set of reusable parameters.

In particular, `Q_{s,1}=0`, so the last off-front residual is a simple difference of inverse envelopes weighted by `Q_{s,2}`. The approximate `k/(k+4)` decay discovered earlier is a crude approximation of these exact rational denominators.

## 5. Held-out integer benchmark

The learned profiles and inferred stabilized scales use **only** PT integers through curve degree 40. The reference file is not imported by the predictor. All calculations use `fractions.Fraction`, and any candidate not an exact integer is rejected rather than rounded.

| Model | Exact predictions | Abstentions | Wrong exact predictions |
|---|---:|---:|---:|
| Frozen v0.9 stabilized-scale baseline | 32/65 | 33 | 0 |
| Terminal pulse transport (single-pass recursive) | 41/65 | 24 | 0 |
| Conservative fitted profiles (without cross-pair extrapolation) | 48/65 | 17 | 0 |
| Affine profiles plus residue anti-pairs | **59/65** | **6** | **0** |

Final coverage by shell: `t=6:13/13`, `t=7:13/13`, `t=8:13/13`, `t=9:13/13`, and `t=10:7/13`.

The six remaining missing sites are `(42,10),(43,10),(47,10),(48,10),(52,10),(53,10)`. They require the still-unanchored mature scale of residue pair `{2,3}` at shell depth 10, and one deeper `r=3` profile input.

In an intentionally oracle-assisted sensitivity test, supplying the single withheld integer `PT(43,10)` fixes the mature scale and predicts four additional integers exactly, leaving `(42,10)` unresolved. This makes explicit how one new integer can unlock an entire small family, but is **not** counted as an AI prediction.

## 6. Caveats and best-practice next test

1. **Retrospective feature search:** the affine-in-`1/W` family and cross-residue pairing were chosen while examining the broader data, despite the numerical fit using only earlier integers. The 59/59 score is a research baseline, not a statistically clean test of a previously frozen hypothesis.
2. **Conditional reference:** the integer targets come from the established rank-zero DT/modular plus PT wall-crossing computation. No such formula enters the predictor, but the benchmark cannot certify those external assumptions independently.
3. **Two-point profiles:** some new scales depend on hypotheses fitted with just two observations, potentially fragile beyond depth 10. Full learned provenance is saved for every profile.
4. **Not a full topological-string calculation:** the experiment is near-extremal PT arithmetic, not the all-genus YY ambiguity itself. Converting genuinely new PT predictions to BCOV ambiguity data still requires the exact GV/PT transform and boundary solver.
5. **Replicate in a new direction:** freeze the present hypothesis and test `t=11..15`, or test a new degree window after extending the conditional reference without using those values in model selection. More forward degrees at previously stabilized depths are less discriminating than deeper transient shells.

A further useful target is to learn the **minimum number of independent new integer anchors per shell**, rather than every coefficient. In this experiment, the six unresolved depth-10 sites already reduce to an unanchored mature scale plus one new depth-3 profile value.

## Reproducibility

Run `python predictor/pulse_profile.py`, then `python validator/validate.py`, then `python tests.py`. The test suite creates an isolated temporary predictor-only directory, with all reference/validator files absent, and verifies byte-for-byte regeneration of the prediction CSV and learned rational parameters. `python validator/one_seed_audit.py` performs the explicitly oracle-assisted sensitivity test separately. Optional marimo and Matplotlib materials are included.
