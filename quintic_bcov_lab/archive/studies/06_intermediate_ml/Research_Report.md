# A small-data AI benchmark for the intermediate quintic holomorphic ambiguity

**Date:** 27 October 2026  
**Status:** Executed experiments; preliminary and conditional predictions, not new exact ambiguity reconstructions.

## Question

Does concentrating on the physical, *intermediate* YY ambiguity directions improve numerical prediction, compared with the first version of our unrestricted holomorphic-ambiguity learner?

We use the exact decomposition

\[
\delta a_g(X)=X(1-X)^{m_g}\sum_{r=0}^{D_g-1}c_{g,r}X^r,
\quad m_g=\left\lceil\frac{3g-3}{5}\right\rceil,
\quad D_g=\left\lfloor\frac{2g-2}{5}\right\rfloor.
\]

In the section frame \(X=x/(x-1)\), the same data are the coefficients of

\[
\delta f_g(x)=\sum_{r=0}^{D_g-1}b_{g,r}x^{r+1}(1-x)^{D_g-1-r},
\quad b_{g,r}=(-1)^{g-2-r}c_{g,r}.
\]

This is a precise change of *homogeneous boundary-kernel basis*, not an \(x\mapsto1-x\) physical symmetry.

## Data and leakage policy

- Reference corrected YY amplitudes in genera **2–12**; **22** free coefficient labels in genera **4–12** (nine genera, not 22 independent statistical observations).
- Exact rational labels are stored as strings in `data/intermediate_training.csv` and `data/intermediate_exact.json`.
- At each free slot \(r\), the only mixed-YY data used by the simple models are summaries of the BCOV **particular/anomaly** polynomial at the pivot degree \(j=m_g+1+r\).
- The target genus true YY coefficients, its reference-derived zero-pivot anchor, and true midpoint values are excluded from features. The source polynomial at the target genus *does*, however, depend on the physically selected lower-genus amplitudes. All predictions are conditional on those inputs and **do not constitute a sequential bootstrap**.
- The canonical zero-pivot anchor for the released labels is reference-derived, not an independently computed higher-genus affine section. It is used for labeling, not for fitting the test observations.
- Data at genera 2,3 have no residual direction under the chosen boundary equations and are not artificially included in the learner.

## Model comparison: fixed forward holdout

Train on **genera 4–8** (eight coefficient observations), predict entire **genera 9–12** (14 coefficients, four genus-level observations). This is the primary holdout; no individual coefficients from any test genus enter the fit.

The magnitude task predicts \(L_{g,r}=\log_{10}|c_{g,r}|\). Magnitude accuracy is scored **without using the true sign**. The sign task is reported separately.

| Model | Magnitude MAE (dex) | Median absolute magnitude error (dex) | Genus-sign accuracy |
|---|---:|---:|---:|
| Linear in genus and residual index | 0.602 | 0.559 | 2/4 with majority sign |
| Linear in mixed-YY maximum \(M\) and index | 0.280 | 0.268 | 2/4 with majority sign |
| **Linear in \(M\), genus, index** | **0.251** | **0.231** | **2/4** with majority sign |
| Fixed slope 1 in \(M\), linear residual correction | 0.251 | 0.232 | 2/4 with majority sign |
| Original v0.1 ExtraTrees on signed logarithm | 5.994 | 6.303 | 2/4 |

Here \(M_{g,r}=\log_{10}\max |\text{mixed nonholomorphic coefficient at }X^{m_g+1+r}|\). The fixed-slope model is

\[
\widehat{L}_{g,r}=M_{g,r}+\alpha+\beta g+\gamma r.
\]

For training through genus 8, the fitted values are

\[
\alpha=0.2266305194,\qquad\beta=-0.1859568352,\qquad\gamma=-0.3136230034.
\]

The magnitude median error of 0.231 dex corresponds to a multiplicative factor \(10^{0.231}\approx1.70\). This is **not** 70% relative error in the signed coefficient; it measures *orders of magnitude*. On unseen residual positions \(r=2,3\) (six observations), the mean magnitude error of the flexible \(M,g,r\) model is 0.215 dex.

The unrestricted v0.1 ExtraTrees regression was designed for descriptive signs/magnitudes over many more coefficients. It extrapolates poorly here, because tree regression does not extrapolate beyond the training magnitude range and because each genus supplies so few truly free directions.

## Chronological robustness

| Data protocol | Model | Mean absolute magnitude error (dex) | Median magnitude factor |
|---|---|---:|---:|
| Train through g=8, test g=9–12 | \(M,g,r\) | 0.251 | 1.70 |
| Train through g=9, test g=10–12 | \(M,g,r\) | 0.120 | 1.32 |
| Expanding origin: train g<target, target g=9,10,11,12 | \(M,g,r\) | **0.116** | **1.36** |
| Expanding origin: same targets | Genus and index only | 0.254 | 1.72 |
| Expanding origin: same targets | v0.1 signed-log ExtraTrees | 2.251 | 195.35 |

The median magnitude factor is \(10^{\operatorname{median}(|\widehat L-L|)}\). Refitting from earlier samples in the expanding-origin test uses **true earlier labels** (teacher forcing): it is a forward holdout, not a self-consistent generated sequence.

## Within-genus shape

All **nine** available nontrivial genera have \(\operatorname{sgn}c_{g,r}\) independent of \(r\). Their absolute magnitudes also decrease approximately geometrically in \(r\).

Fitting only the three \(r=1\) ratios at genera 6–8 gives

\[
\log_{10}\left|\frac{c_{g,r}}{c_{g,0}}\right|\approx -0.466308\,r.
\]

Against the ten higher-index observations in genera 9–12, the median absolute shape error is **0.0856 dex** (mean 0.0969 dex). This **uses the true first coefficient** \(c_{g,0}\) at each test genus and is therefore **an oracle-conditioned structural diagnostic, not blind inference**. It is encouraging evidence for predicting one scale, one sign, and a simple within-genus shape.

## Sign: the hard unsolved component

A single genus sign determines all \(D_g\) coefficient signs in the current dataset, but none of the standard local mixed-anomaly sign summaries predicts that genus sign reliably. A majority-sign rule and an ExtraTrees sign classifier score **50% by coefficient** on the primary held-out split (two of four genera).

The short sequence of observed genus signs is

```
 genus:  4  5  6  7  8  9 10 11 12
 sign:   -  +  +  -  +  +  -  -  +
```

A period-three sign hypothesis suggested by the initial five training genera correctly classifies **three of four** held-out genera (10 of 14 coefficients), but **fails at genus 11**. This pattern was noticed retrospectively; it is **not** evidence of an established periodicity or a valid blind discovery. The four test genera are too few to judge statistical significance.

Sign errors dominate fully signed-log metrics (the model with accurate magnitudes but majority signs has signed-log MAE 19.13). The arithmetic lattice/GV boundary decoder therefore remains the plausible route to discrete sign and exact rational choices.

## Positional embeddings and rational representation

Exact rational targets are retained without float truncation; only fitted logarithms are floats. A reasonable next architecture would predict a genus-wide sign, \(\log|c_{g,0}|\) using the YY source scale, and a low-rank shape \(\log|c_{g,r}/c_{g,0}|\).

The current first training split has only indices \(r=0,1\). The \(r=2,3\) targets occur only in future genera. A GV/RoPE feature dictionary over these positions is thus **underdetermined** without stronger analytic priors. We did not claim that training a rotary-attention network on 22 numbers would give a meaningful independent test.

Likewise, the exact integer/GV congruence decoder cannot be scored blindly at higher genus without independent denominator/charge-lattice bounds. Supplying true denominators would leak information.

## Verification

- `python run_experiment.py` regenerates `results/benchmark.json` and `results/predictions.csv` and fits only chronological training data.
- `python tests.py` checks the exact \(c\leftrightarrow b\) sign/basis identity, target count, no forbidden feature columns, and out-of-genus splits.
- `cd source && python tests_intermediate.py` independently checks the original reference section round-trip identities (11 genera and 22 targets). Optional `--boundary --max-genus 9` additionally invokes the full exact boundary-system solver.
- `marimo edit intermediate_ml_marimo.py` launches an interactive version of the fixed-split comparisons and shape diagnostic. The marimo GUI was not run in this container; Python syntax and the numerical backend were tested.

## Conclusion

**Magnitude extrapolation becomes significantly better when we factor out the local BCOV mixed-term scale and fit the residual index**, but the sign and arithmetic lattice representative are not learned. The best next target is a two-stage model: genus-wide sign/charge-lattice offset followed by cross-index shape, with exact verification against the independently computed boundary system.
