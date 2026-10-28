# Quintic holomorphic ambiguity: AI-assisted predictor (experimental v0.1)

**Status:** Runnable research prototype with exact B-model constraint projection,
integer-coset decoding, a learned low-X long-genus predictor, and a conditional
antiholomorphic-feature benchmark. Exact YY reference through genus 12 is included;
low-X coefficients through genus 1000 come from the separately audited BCOV quotient
computation. No unpublished genus-1001 rational coefficient is claimed to be exact.

## Mathematical conventions and provenance

- The **canonical** original YY convention here is `x=psi^-5`, `X=x/(x-1)` and
  `P_g = P_g^an + sum_j a[g,j] X**j`, `P_g^an(0,0,0,X)=0`.
- **Do not mix** the HKQ, `z`, and YY `psi` conventions without the complete
  transformation of the polynomial and generators.
- The genus-four reference uses our previously audited GV-corrected polynomial,
  `P4_corrected = P4_YY_printed + (6746/5**12)*X**5*(1-X)**2`,
  rather than the three discrepant pure-X terms of the original printed appendix.
- Mixed terms determine an anomaly particular solution, **not** the same-genus
  ambiguity. The B-model conifold gap, orbifold regularity, and MUM give Q-linear
  constraints on the entire ambiguity. GV input is separate and explicitly tagged.
- A previous exact BCOV computation fixes `a[g,0..4]` through genus 1000 using
  lower-genus seeds, constant maps and known degree<=4 GV vanishings; it is not
  a proof of uniqueness from integrality or an AI prediction.
- The known **leading Airy** number `alpha_g` is an independent feature. Higher
  speculative Airy descendants and graded Ward relations have not been
  silently added as exact training labels; they require explicit verified files.

## Quick start

Requires Python 3.10+, NumPy, scikit-learn (for anomaly benchmark); exact
boundary and lattice modules themselves use standard-library `Fraction`.
Run from this package's root:

```bash
python -m hybrid.tests
python -m hybrid.model --data results/gmp/direct_J4_exact.json.gz --output hybrid/results
python -m hybrid.predict --genus 1001 --train-end 1000
python - <<'PY'
from hybrid.anomaly_learning import benchmark, experimental_full_projection
print(benchmark('reference_amplitudes.json'))
print(benchmark('reference_amplitudes.json','results/all_constant_and_leading.json.gz'))
print(experimental_full_projection('reference_amplitudes.json',9))
PY
```

The complete exact dataset may take around 30–60 seconds to parse because some
rational numerators have many thousands of digits. The separately cached
`hybrid/results/float_dataset.npz` makes model experiments and smoke tests
fast; **never** treat this floating cache as an exact arithmetic certificate.

## Pipeline and implemented methods

1. **Exact teacher:** the original BCOV/YY sparse polynomial engine computes
   the particular solution; the reference dataset and degree-4 quotient data
   supply checks. Store the source of every boundary condition.
2. **Physical constraint matrix:** `hybrid.boundary.boundary_system` builds
   exact conifold, orbifold and selected MUM rows and returns a particular
   solution plus a certified rational kernel. `mum_degrees=[0]` withholds all
   positive-degree enumerative information, leaving a true test of the learner.
3. **Arithmetic encoding:** `hybrid.digits` losslessly encodes sign, base-1000
   numerator/denominator limbs and small-prime valuations. These are suitable
   for token-based sequence models, unlike direct floats of gigantic rational
   coefficients. `hybrid.lattice` solves exact integrality congruences for a
   proposed denominator cap, then selects the coset member closest to a
   predicted BPS index; the cap and selection are extra assumptions.
4. **Learned continuous prior:** `hybrid.model` models the first four
   coefficients of `log(a_g(X)/a[g,0])` using a multioutput fit on
   `1,g^-1,g^-2` (after factoring an overall genus). This is an arithmetic-aware
   shallow regression, not a transformer or exact all-genus identity.
5. **Mixed-term learner:** `hybrid.anomaly_learning` extracts local slices of
   the anomaly tensor indexed by the X degree, generator exponents, and
   relative coefficient scale, then fits an ExtraTrees model with strict
   held-out future genera. It is considerably weaker than physics-first
   low-X fitting, but captures nontrivial predictive signs.
6. **Certification:** only the exact BCOV residual, exact boundary equations,
   and independently justified integer/GV constraints can certify the result.
   A decimal coefficient forecast is never automatically rounded to a large
   rational. The output type keeps forecasts distinct from certificates.

## Quantitative tests (run in this prototype)

- **Genus 54–1000 holdout**, trained ONLY on g40..53 (14 genera), four
  predicted ratios j=1..4: median relative error `0.000455783554` versus
  `0.158567` for the simple binomial model.
- **With exact a1 supplied from the independent Bernoulli identity**, g54..1000
  predicted j=2,3,4: median relative error `9.76641072e-7`. a0, a1 are *not*
  scored as predictions. The all-genus formula for a1 was independently checked
  at g54, 200, 500, 1000 against the cached table.
- **Antiholomorphic-feature baseline**, train g2..8, test g9..12:
  sign accuracy `0.8070` (constant baseline `0.5965`); signed-log MAE `3.0148`
  with independent leading Airy features (without `3.0310`). This is not
  numerically precise enough to reconstruct rational coefficients by itself.
- **Genus-9 partial physics projection**, conifold + orbifold + MUM constant:
  rank 22/25, free slots (6,7,8). It enforces every supplied boundary
  constraint exactly. The mixed-term learner predicts the free slots poorly;
  some remaining coefficients have enormous relative error. The observed
  17/25 exact matches after projection arise from the *physics constraints*,
  not successful AI extrapolation.
- **Exact arithmetic control**: the genus-3 integral coset calculation
  correctly recovers a3,4 = `235625/2016`, a3,5 = `33125/1008`, and a3,6 =
  `15625/1008` with previously supplied empirical denominator `D3=36288`.
- **Exact boundary controls**: full physical reconstructions through genus 10
  match the archived corrected ambiguity vectors; supplied light tests check
  g2, g3, g4 and the partial rank at g7.

## Why g54 is a practical new target

With the usual counting, the ambiguity has `3g-2` rational slots. For the
quintic the conifold gap imposes `2g-2` equations, orbifold regularity another
`ceil((3g-3)/5)`, and the remaining MUM slots correspond to degrees
`0..floor((2g-2)/5)`. Geometric Castelnuovo vanishing establishes
`n_{g,d}=0` when `10g > d*d+5d+10`.

At g54, degree d=21 is the first unprotected primitive-GV input after the
known special g51,d20 entry; at g55 the same degree remains unprotected. In
contrast, the g54 **first five X coefficients** are already determined by
low-degree vanishings and the exact BCOV quotient solver. The next research
problem is therefore to learn the *small missing BPS index*, not all 160
coefficients in an arbitrary basis. At g1000 there are 302 potentially
nonvanishing positive-degree MUM inputs among d<=399 (before any further
MSW/modular/DT/PT constraints). These are counts of possible free enumerative
inputs, not a proof of linear independence for the entire boundary system.

## Scaling and next AI experiments

- Build a dataset of `N_{g,d}` near the Castelnuovo edge with **degree and
  genus held out together**, from independently certified GV/PT tables.
- Model in triangular LCS/BPS coordinates, with a **classification head for
  exact zero**, signed-log-magnitude head, and a discrete **integer-coset
  offset head**; never regress the final thousands-digit numerator directly.
- Start with CatBoost/LightGBM-style tabular models and a 1D sequence model;
  only consider a transformer after there are enough nonzero edge labels.
- For a digit decoder use signed base-1000 limbs, explicit place embeddings,
  denominators/prime valuations, and grammar-constrained output. Employ exact
  CRT/modular reconstruction and `Fraction` verification after decoding.
- Incorporate leading Airy normalization now. Feed higher Airy-dressed layers
  and Ward/graded relations as **separate ablations**, once their normalization
  and certified rank are audited. Include an optional MSW conditional source
  but do not train on it in the nominal pure-B-model experiment.
- Measure (a) exact-coset match, (b) GV exact match, (c) uncertainty calibration,
  (d) genuses/curve degrees withheld, (e) sensitivity to YY gauge changes, and
  (f) violations of conifold/orbifold/MUM equations. Avoid row-level random
  splits that leak the same genus into both train and test.

## Literature / independent context

YY, *Topological String Partition Functions as Polynomials* (2004),
https://arxiv.org/abs/hep-th/0406078 ;
Alim–Yau–Zhou, *Airy Equation for the Topological String Partition Function*
(2015), https://arxiv.org/abs/1506.01375 ;
McLeish et al., *Transformers Can Do Arithmetic with the Right Embeddings*
(NeurIPS 2024), https://proceedings.neurips.cc/paper_files/paper/2024/hash/c35986bc1ee29b31c1011481b77fe540-Abstract-Conference.html ;
Schwartz et al., *NumeroLogic* (EMNLP 2024),
https://aclanthology.org/2024.emnlp-main.12/ .
