# A controlled two-stream GV–YY experiment for the high-degree quintic ambiguity

**Computational note, 24 October 2026**  
**Status:** exact boundary reconstruction and genuine small-data supervised extrapolation, with negative controls. **Not** a trained large neural transformer, discovery of an all-genus formula, or new rational genus-13 reconstruction.

## Summary

1. The exact GV/MUM response matrix is *strictly triangular* in the ambiguity index. In the canonical YY convention, `X=x/(x-1)`, `q=x/3125+...`, so the response of the q^d amplitude to a pure-X coefficient a_{g,j} is identically zero for `j>d`. This is verified over exact fractions, through d=4, for genera 6, 9, 10, 11, 12. The diagonal equals

   \[R_{g,d,d}=(-1)^{g-1+d}5^{5d-4g+4}.\]

   At genus 12 this leaves all 29 columns `j=5,...,33` in the exact degree-4 GV-projection kernel.

2. Nonetheless, the **full exact B-model boundary map** (conifold gap + orbifold regularity + MUM constant map), applied to the *known physical lower-genus BCOV anomaly source*, fixes every high-X coefficient of genera 10, 11, and 12. All fixed entries match their reference **rational numbers exactly**, without positive-degree same-genus GV input.

3. A small-data supervised surrogate learns to extrapolate the high-X *shape* from the mixed `v_i X^j` (antiholomorphic) terms, with modest improvement over inverse-genus polynomial regression on a broad training window. Its additional physically aligned GV channel reduces extreme coefficient errors for that window, but **does not consistently improve the median error**. RoPE-like GV phase features do not improve it.

4. A known stronger inverse-genus extrapolation fitted to a short, more recent window (g=6..9) remains hard to beat. This window sensitivity makes it inappropriate to declare victory for the more complex models.

## 1. Exact physical inputs and constraints

The archive is the corrected canonical YY reference for genera 2–12:
\[
P_g=P_g^{\mathrm{an}}(v_1,v_2,v_3,X)+a_g(X),\quad a_g(X)=\sum_{j=0}^{3g-3}a_{g,j}X^j.
\]

The genus-four entry is the GV-corrected polynomial, not the discrepant printed YY Appendix value. All computations are over Python `Fraction` for the physics; floating point is used only in learned regressors and Fourier-feature diagnostics.

We use the exact quintic period frame, mirror map, and GV/GW multicover kernel:
\[
F(\lambda,t)=\sum_{h,d,m\ge1}\frac{n_{h,d}}m Q^{md}\left(2\sin\frac{m\lambda}2\right)^{2h-2},
\]
with the genus-zero sector incorporated by its Laurent expansion, not by an inappropriate finite-character approximation. At degrees 1–4, the BPS inputs for h=0,1 are computed from exact periods and genus-one data; h=2,3,4 are calculated from the **training-side known low-genus physical amplitudes**, and genus h≥5 vanishes by the established quintic geometric genus bound. The nonzero degree-four spin data include n_{2,4}=534750 and n_{3,4}=8625.

All five `q^0,...,q^4` coefficients of the full physical amplitude reproduce the multi-cover prediction *exactly* at every genus 4–12 (nine checks). No amplitude at genus 10–12 is used to construct the GV feature dictionary; those reference amplitudes are consulted only for the exact correctness audit and final scoring. Although degree four is more MUM information than is strictly required at some earlier genera, its vanishing and lower-genus contributions are independently available.

### MUM rank-zero result for high X

Since `X=-x+O(x^2)` and `q=x/3125+O(x^2)`, while the frame factor is nonzero at the MUM point,
\[
\frac{\partial [q^d]F_g}{\partial a_{g,j}}=[q^d]\big(k(x)^{g-1}X^j\big)=0\quad(j>d).
\]
The diagonal is the exact signed power of five above. The q^0..q^4 Jacobian has rank **5**, with columns j>4 exactly zero, for all verified genera 6,9,10,11,12. This is independent of machine learning and is not an artefact of numerical precision. Its significance is that the multi-cover GV stream cannot supply **direct** training labels for the distant high-degree coefficients. It can influence them only indirectly through lower-genus BCOV recursion or by supplementing a boundary map.

### Exact high-X determination by conifold and orbifold

| genus | linear equations rank | free coefficient indices | exactly fixed high-X range | verification |
|---:|---:|---|---|---|
| 10 | 25 of 28 | 7, 8, 9 | 10–27 | all entries exact |
| 11 | 27 of 31 | 7, 8, 9, 10 | 11–30 | all entries exact |
| 12 | 30 of 34 | 8, 9, 10, 11 | 12–33 | all entries exact |

The exact rational equations impose **conifold gap, orbifold regularity and MUM constant map only**, with the previously physically selected lower-genus amplitudes used by BCOV. Unknown same-genus positive-degree GV values remain free. The fixed high-X range contains no components of any exact kernel vector. The predictions match the physical reference precisely; this is a constraint solution, **not neural extrapolation**.

## 2. What was learned

For K=3g−3, consider the reversed high-X polynomial
\[
B_g(y)=\frac{y^Ka_g(1/y)}{a_{g,K}},\qquad
C_g(y)=\log B_g(y)=\sum_{r\ge1}c_{g,r}y^r.
\]
We predict r=1,...,6 by regressing `c_{g,r}/g`, and reconstruct `b_{g,r}=[y^r]B_g` with the exact triangular exponential recurrence. All models use the physical **top coefficient `a_{g,K}` as an oracle** to normalize the target; they are tests of shape, not independent absolute ambiguity reconstructions.

**Inputs:**

- **Inverse-genus baseline:** `1, 1/g, 1/g^2, 1/g^3`.
- **YY stream:** the corresponding reverse-logarithmic coefficient in the `v1 X^{K−1−r}` mixed channel of the known BCOV anomaly polynomial. This does not use any held-out ambiguity coefficient.
- **Physical GV/MUM stream:** exactly calculated GV-to-GW `q^1,...,q^4` data and the physical LCS residual `target_GW - [q^d] F_g^{an}`. The supervised one-dimensional residual feature is the signed decimal log of its degree-four component. No ad hoc identification of BPS curve degree with YY exponent is made.
- **GV phase stream:** Fourier sine-character values with indices `(h,d,m)` exactly those of the GV multi-cover formula, sampled at several inverse-square-root-genus couplings. It is a **small character-feature dictionary**, not a trained RoPE attention operator. The h=0 channel uses its rational Laurent kernel elsewhere in the exact GW transform. Positive h=1 gives a large genus-independent constant which dominates this particular Fourier sampling; this ablation is therefore weak evidence against *other*, spin-normalized GV representations.

The learned residuals are ridge-regularized after projecting out the inverse-genus basis, so that additional channels compete on improvements **beyond** the baseline. Hyperparameters are selected exclusively on forward g=8,9 validation; g=10,11,12 are withheld until evaluation. The training data consist of only six independent genera (g=4–9), 36 `(genus,tail-position)` observations; results must be interpreted with this sample-size limitation.

## 3. Strict forward-genus heldout results

Train g=4..9; tune only using past forward validations; evaluate 18 unseen `(g,r)` samples at g=10..12. Lower is better.

| model | median coefficient relative error | p90 coefficient relative error |
|---|---:|---:|
| inverse-genus, quadratic | 1.72968% | 91.6748% |
| inverse-genus, cubic | 0.33355% | 82.5729% |
| cubic + YY mixed | **0.12469%** | 55.9520% |
| cubic + GV/MUM | 0.32940% | 81.7882% |
| cubic + YY + GV/MUM | 0.28247% | **12.3200%** |
| cubic + YY + GV phase | 0.12469% | 55.9520% |
| cubic + YY + GV/MUM + GV phase | 0.28247% | 12.3200% |

The broad-window YY model improves the median by 2.7x, but the YY+GV model does **not** improve that median. The combined model notably improves the upper error tail. Its extra phase inputs have no observable incremental effect at the displayed precision.

The particularly large p90 errors originate mostly from the sixth (farthest) high-X coefficient, which behaves very differently at genus four: c_{4,6}/4≈−97.24, versus 0.47 at genus five. This makes the broad-window extrapolation unusually difficult.

### Window sensitivity, the critical control

A different, more physically natural window `g=5..9` substantially improves every method. The median / p90 errors are:

| train window; model | median | p90 |
|---|---:|---:|
| g=5..9; cubic baseline | 0.07605% | 2.45812% |
| g=5..9; cubic + YY | 0.04617% | 1.00630% |
| g=5..9; cubic + YY + GV | **0.02593%** | **0.15513%** |
| g=6..9; cubic baseline | **0.02288%** | 0.51405% |

The g=6..9 cubic baseline interpolates four samples per position and is the *pre-existing, very strong* high-X benchmark. Adding learned residual streams to that exactly interpolated four-sample fit supplies no additional improvement without re-regularizing the genus basis. The `g=5..9` combined model improves the p90 versus that strong benchmark, but not its median. Because multiple windows have been inspected, this comparison is exploratory, not a pre-registered superiority claim.

### Permutation control

To test whether YY inputs contain information rather than being an arbitrary correlated regressor, we randomly permuted the training-genus YY signals 199 times (seed 13579), retaining the true held-out features. The unshuffled g=4..9 YY model has a median 0.12469% relative error; **about 95%** of the 199 shuffled controls do worse on median. The resulting approximate one-sided randomization p-value is 0.055. This is suggestive rather than conclusive, and multiple investigated models further caution against statistical overinterpretation.

## 4. GV/RoPE mathematical identity and the negative result

For a BPS sector h≥1, with s=h−1,
\[
\left(2\sin\frac{m\lambda}2\right)^{2s}
=\sum_{r=0}^{2s}(-1)^{s+r}\binom{2s}{r}e^{im(s-r)\lambda}.
\]
This equality is a direct finite character decomposition, with weighted Fourier rotations exactly analogous to the complex unit rotations in RoPE. The test verifies it independently to machine precision across h=1,..,4, m=1,..,4 and three angles; the earlier project's exact integer Laurent-polynomial identity provides the stronger exact version. But it does **not** identify the GV spin/multicover mode with the YY power index X^j.

In this initial physically aligned *feature-bank* experiment, charge-phase features were dominated by the h=1 constant channel and contributed no incremental performance. One must not infer that a trained character-valued attention model, a subtracted-spin character basis, or a digit-place positional embedding cannot improve with more data. Such models have **not** been trained here.

## 5. Conclusions and most productive next test

**Exact mathematical:** The high-X ambiguity of g10–12 is already fixed by the B-model gap and orbifold equations once the previous physically selected BCOV amplitudes are known. The low-degree GV->GW stream has identically zero direct Jacobian on that high-X range.

**Empirical AI:** YY anomaly mixed-term features help median high-X extrapolation in a six-genus broad-window comparison. A GV/MUM feature sometimes suppresses worst-case error, but window-sensitive results do not establish robust superiority over the existing inverse-genus baseline. Character phases by themselves add no measured predictive value in this pilot.

**Recommended actual next architecture:** train an explicitly *charge-space* sequence model on much larger algebraically generated, physically consistent BPS/GW data, with (i) exact multicover coupling of `(d,m,h)` labels, (ii) digit-place/prime-valuation rational representations, (iii) spin-subtracted GV characters, (iv) exact sparse conifold/orbifold projection or a learned surrogate to those matrices, and (v) a **fully separate heldout-genera/charge** evaluation. Use the current finite reference as calibration and exact verification, not sufficient training data for a deep model. Never credit neural predictions for coefficients forced by the exact boundary projector.

## Reproduction

See `README.md`; the principal commands are `python run_experiment.py`, `python statistical_controls.py`, `python exact_response_audit.py`, `python verify_high_boundary.py`. No internet required. Source and exact genus-2..12 reference, all per-genus predictions, boundary ranks, permutation summaries, and numerical logs are included in the bundle.
