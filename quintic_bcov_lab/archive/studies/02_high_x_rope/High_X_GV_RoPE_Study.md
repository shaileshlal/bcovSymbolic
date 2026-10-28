# Upper-degree quintic holomorphic ambiguity and GV-inspired positional embeddings

*Computational pilot, 24 October 2026. Exact reference: canonical Yamaguchi–Yau gauge with our audited genus-four correction. Calculations are reproducible with the attached code.*

## Executive findings

1. **Exact boundary projection already fixes the upper YY ambiguity conditional on the lower-genus BCOV solution.** With conifold gap, orbifold regularity, and only the large-complex-structure (MUM) constant-map condition, the exact affine solution spaces for genera 7, 9 and 10 have respectively free indices `(5,6)`, `(6,7,8)` and `(7,8,9)`. All coefficients beyond the largest free index agree with their exact reference values. Neither upper-X coefficient machine learning nor a positive-degree GV oracle contributes to these exact statements. The genus-10 boundary calculation is moderately costly and is available through `verify_boundary.py --all`.
2. A **reversed polynomial logarithm** gives smoothly varying high-degree coefficients, permitting substantial short-range genus extrapolation.
3. The GV multi-cover sine factor is an **exact finite sum of rotary complex phases** at each fixed positive BPS genus; this is a real representation-theoretic analogy with RoPE. However, identifying its curve degree and spin frequencies with YY *X-coefficient index* is not mathematically justified.
4. Our compact GV-character and compact RoPE-like feature banks **do not outperform** index-aligned inverse-genus regression in an honest (small) future-genus holdout. Full-rank phase feature banks exactly reproduce the aligned baseline: a useful change-of-basis sanity check.

## 1. Conventions and exact upper-degree determinacy

Use

\[
X=\frac{x}{x-1},\quad x=\psi^{-5},\quad
P_g=P_g^{\mathrm{an}}+a_g(X),\quad
 a_g(X)=\sum_{j=0}^{K} a_{g,j}X^j,\quad K=3g-3.
\]

The particular solution is normalized by `P_g^an(0,0,0,X)=0`. The exact physical boundary rows depend linearly on the full unknown vector `a_g`. They require the *already physically chosen* lower-genus amplitudes to compute `P_g^an`. We impose the conifold gap, orbifold regularity and the MUM constant map; no positive-degree GV input or MSW constraints are used. `a_g,0` is exactly supplied by the Bernoulli constant-map formula.

| Genus | Number of slots | Exact rank | Free indices | Upper fixed indices |
|---|---:|---:|---|---|
| 7 | 19 | 17 | 5,6 | 7–18 |
| 9 | 25 | 22 | 6,7,8 | 9–24 |
| 10 | 28 | 25 | 7,8,9 | 10–27 |

For example, the physical genus-nine reference has

\[
a_{9,24}=\frac{2614676952362060546875}{229824},\quad
 a_{9,23}=-\frac{768059253871440887451171875}{11732745024}.
\]

Both entries are reproduced *exactly* from the boundary linear system with all nonconstant MUM data withheld, as are all other upper entries in the range specified in the table. This is a condition-on-lower-genus result, not an independent all-genus closed formula.

## 2. Reverse logarithmic cumulants: a useful gauge of high-degree structure

Introduce the reversed polynomial and its formal logarithm:

\[
 B_g(y)=\frac{y^{3g-3}a_g(1/y)}{a_{g,3g-3}}
 =1+\sum_{r\geq1}b_{g,r}y^r,
 \qquad C_g(y)=\log B_g(y)=\sum_{r\geq1}c_{g,r}y^r.
\]

Then `b_g,r=a_g,K-r/a_g,K`. The transformation is *exact, triangular and invertible*: for `r>=1`,

\[
 c_{g,r}=b_{g,r}-\frac1r\sum_{k=1}^{r-1}k\,c_{g,k}\,b_{g,r-k},\qquad
 b_{g,r}=\frac1r\sum_{k=1}^{r}k\,c_{g,k}\,b_{g,r-k}.
\]

The first six coefficients vary smoothly across higher available genera. We tested, **with no genus 10–12 training data**, the model

\[
 \frac{c_{g,r}}{g}=\beta_r+\gamma_r g^{-1}+\delta_r g^{-2}+\epsilon_r g^{-3},
 \qquad r=1,\ldots,6.
\]

Four training genera (6–9) give four points for each four-parameter curve: this is interpolation, *not* an overdetermined statistical fit. Three unseen genera (10–12), comprising 18 tail coefficients, provide the actual evaluation. All reconstructed tail values are compared against exact rational data, but prediction itself is numerical.

### Predeclared holdout: train g=6–9, evaluate g=10–12

**Target:** relative tail `b_g,1` through `b_g,6`; the true top coefficient `a_g,K` is supplied for this shape-only test, which is explicitly an oracle normalization. Relative errors are computed on `b` after cumulant-to-coefficient reconstruction, not on the easier `c` target.

| Index / genus representation | Median relative error | 90th-percentile relative error |
|---|---:|---:|
| Position-aligned, linear in g | 40.9448% | 125.9922% |
| Position-aligned, `1,1/g,1/g²` | 0.1692% | 3.1169% |
| **Position-aligned, `1,1/g,1/g²,1/g³`** | **0.02288%** | **0.5141%** |
| Polynomial in index + inverse-g² | 15.6872% | 51.3885% |
| Compact RoPE-like Fourier features + inverse-g² | 3.9733% | 21.1619% |
| Compact GV sine-character features + inverse-g² | 16.4179% | 64.4383% |
| Full-rank RoPE phase bank + inverse-g³ | 0.02288% | 0.5141% |
| Full-rank GV character bank + inverse-g³ | 0.02288% | 0.5141% |

**Interpretation:** Full-rank phase representations are invertible coordinate changes on the six discrete index positions and therefore carry *no extra information*. Compressed phase dictionaries performed worse, which may partly reflect unequal capacity/regularization; this is not a negative theorem about attention with RoPE. The strong result in this small experiment is the inverse-genus structure after reverse-cumulant normalization. Model choices were explored on the available data, so these numbers are exploratory rather than an unbiased final estimator of generalization.

A second evaluation predicts the top scale `a_g,K/a_g,0` by **affine regression on its base-10 logarithm**, using only g=6–9, then predicts all seven high-X ratios including `k=0`. Only the exact constant-map `a_g,0` is supplied; the top coefficient is predicted. Across g=10–12 (21 coefficients), the median relative error is **0.08113%** and 90th-percentile error is **0.40453%**. The genus-12 median is 0.12117%. These are still far from exact rational recovery.

## 3. GV as a finite phase representation: exact relationship to RoPE

The GV multi-cover formula (suppressing classical and degree-zero terms) is

\[
 \mathcal F(\lambda,t)=\sum_{h\geq0}\sum_{d\geq1} n_{h,d}\sum_{m\geq1}\frac{Q^{md}}{m}
 \left(2\sin\frac{m\lambda}{2}\right)^{2h-2}.
\]

For `h>=1` let `s=h-1`. Then **exactly**

\[
 \left(2\sin\frac{m\lambda}{2}\right)^{2s}
 =\sum_{r=0}^{2s}(-1)^{s+r}\binom{2s}{r}
 e^{\mathrm i m(s-r)\lambda}.
\]

This identity is verified in integer Laurent-polynomial arithmetic (`verify_phases_and_digits.py`), and separately numerically for genus-channel exponent `s=1..5` at multiple angles. In RoPE, a two-dimensional query/key plane is rotated as `e^{i p omega}` and the query-key product generates `e^{i(q-p)omega}`. GV uses the *same unitary Fourier characters* `e^{i m\ell\lambda}`, but with a specific alternating binomial character weight and an additional multi-cover arithmetic map `(d,m) -> md`, `Q^{md}/m`.

It is important not to overstate the analogy: the ordinary genus-zero (`h=0`) GV sector has negative exponent and is not a finite character polynomial; its sine inverse requires a different Laurent/rational treatment. More importantly, a Fourier mode of string coupling or curve degree is **not** the index of an `X` coefficient in a YY polynomial. Fourier features should be tested on the actual GV-to-GW multicover operator and BPS charge lattice, where their representation theoretic meaning is exact, rather than identified with the ambiguity's X powers without a mirror-frame transfer map.

### Proposed physically aligned dual-position embedding

Use the following separate streams in the eventual trainable architecture:

- **Physical charge stream:** `(d,m, h)` with Fourier modes `exp(i m ell lambda)` and integer BPS spin-character weights. Include divisor/multicover incidence `m*d`, the `1/m` weight, and genus/Taylor order. This can be implemented as structured relative-phase attention or a sparse convolution. The Fourier/character frequencies arise directly from GV.
- **YY polynomial stream:** `(g,j,K-j, v_1,v_2,v_3)` exponents with index-aligned or learned embeddings, conditioned on computed BCOV mixed terms and boundary rows. Its map to GV charges is the *exact mirror-map and GV coefficient-extraction layer*, not an assumed identification of index `j` with degree `d`.
- **Integer representation stream:** signed numerator/denominator limbs in base 1000, place significance aligned from the least significant end (Abacus-style), plus valuations at small primes and modular residues. This represents enormous rationals without float truncation. Apply exact `Fraction` and congruence checks to every decoded candidate.

An important distinction: Abacus embeddings address *digit-place arithmetic*; RoPE addresses *relative phase*, and GV characters encode *specific physical Fourier amplitudes*. The three notions of position are complementary, not interchangeable. The original Abacus embeddings and position coupling results make the case for explicit place significance; they do not imply that a GV character bank necessarily improves YY coefficient extrapolation.

## 4. Explicit prospective genus-13 forecast (not certified)

The table reports **dimensionless numerical ratios** `a_13,j/a_13,0`, where `a_13,0` has a separate Bernoulli formula, using the inverse-g³ cumulant extrapolation and linear top-scale extrapolation fitted over **g=8–12**. `K=36` at genus 13. These are neither rational reconstructions nor physically boundary-verified coefficients.

| j | Forecast `a_13,j/a_13,0` |
|---:|---:|
| 36 | +2.93533 × 10^8 |
| 35 | −3.48804 × 10^9 |
| 34 | +1.94492 × 10^10 |
| 33 | −6.74864 × 10^10 |
| 32 | +1.62966 × 10^11 |
| 31 | −2.90156 × 10^11 |
| 30 | +3.93648 × 10^11 |

The archive contains forecasts for four nested windows, g=6–12 through g=9–12, to assess sensitivity. Their mutual agreement is *not* a statistical confidence interval; it does not account for shared model bias or extrapolation breakdown. The decisive next test is full genus-13 BCOV + exact conifold/orbifold projection (or an independent physical reference), not more rounds of fitting these same eleven genera.

## 5. What has and has not been learned

**Established with exact checks:** anomaly-conditioned upper-X boundary determinacy for the checked genera; the rational cumulant inversion; the integer Laurent-character identity; lossless digit tokenization; no positive-degree GV or MSW oracle in the partial boundary projection.

**Empirical, genuinely tested:** inverse-genus cumulant holdout on genera 10–12; forecast of the unknown top scale from earlier genera; Fourier/GV-character dictionary ablations against full-rank controls.

**Not demonstrated:** exact new high-degree rational coefficients at g>=13; a general analytic inverse-genus expansion; superior performance from RoPE or GV positional embeddings; determination of the full genus-g ambiguity without lower-genus physical inputs; a new modular or Airy hierarchy identity. The known leading Airy coefficient controls a different, propagator-high-power sector and is not implicitly equated to `a_g,K`.

The next useful neural experiment is to train **on charge-space GV/PT data and exact mirror-map projections**, using the character kernel to handle the arithmetic multi-cover map. The high-X surrogate should then consume the physically projected BPS variables and the lower-genus BCOV anomaly tensor. This tests the analogy in the correct representation and makes possible cross-genus/charge holdouts rather than short interpolation of only eleven polynomials.

## Reproducibility

- `python high_degree_phase_study.py`: pure NumPy, reads the exact published/reconstructed reference, tests models, writes `results_high_degree_phase.json`.
- `python verify_phases_and_digits.py`: pure Python exact Laurent-identity check and lossless rational digit round trips.
- `python verify_boundary.py` (add `--all` for g10): exact rational B-model constraints, verifies the reference against the rank and kernel statement, independently of ML.
- `forecast_genus13.json`: prospective coefficients and four training-window sensitivity checks.
- `exact_high_tail_projection.json`: cached exact genus 7/9/10 rank and verification results.

All raw genus 2–12 amplitudes and source engine files are bundled. The YY genus-four ambiguity is the GV-corrected one; do not replace it with the discrepant printed YY Appendix vector.

## Sources

- S. Yamaguchi and S.-T. Yau, *Topological String Partition Functions as Polynomials*, https://arxiv.org/abs/hep-th/0406078 .
- M. Alim, S.-T. Yau and J. Zhou, *Airy Equation for the Topological String Partition Function*, https://arxiv.org/abs/1506.01375 .
- Gopakumar–Vafa multi-cover formula as reviewed in *Mathematical Structures of Non-perturbative Topological String Theory: From GW to DT Invariants*, https://link.springer.com/article/10.1007/s00220-022-04571-y .
- J. Su et al., *RoFormer: Enhanced Transformer with Rotary Position Embedding*, https://arxiv.org/abs/2104.09864 .
- S. McLeish et al., *Transformers Can Do Arithmetic with the Right Embeddings*, NeurIPS 2024, https://proceedings.neurips.cc/paper_files/paper/2024/hash/c35986bc1ee29b31c1011481b77fe540-Abstract-Conference.html .
- H. Cho et al., *Position Coupling: Improving Length Generalization of Arithmetic Transformers Using Task Structure*, NeurIPS 2024, https://proceedings.neurips.cc/paper_files/paper/2024/hash/27aa3a0e6d63db269977bb2df5607cb8-Abstract-Conference.html .
