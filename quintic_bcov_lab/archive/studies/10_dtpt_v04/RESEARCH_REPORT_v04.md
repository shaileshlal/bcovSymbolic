# How far can empirical DT/PT arithmetic predict the quintic BCOV ambiguity?

**Executed 27 October 2026.** Exact rational/int computations, forward-degree holdouts, rank-zero DT modular reconstruction, and explicit DT→PT→GV wall-crossing tests. The published genus-54 target is **excluded** from parameter fitting. The mathematical feature family was chosen retrospectively, so this is a non-oracular *within-model* holdout, not a preregistered historical forecast.

## Executive result

The first genuinely free residual lattice charge beyond the standard genus-53 bootstrap is

\[\boxed{n_{54,21}=319125.}\]

There are **three reproducible routes** to it:

1. **Empirical PT charge normalization.** Divide each leading PT integer by the signed wall-crossing charge multiplier. Degree-6–10 examples yield just the three integers \(\Omega_0=5,\ \Omega_1=8625,\ \Omega_2=-1218500\). With no more fit parameters, all **20** degree-11–30 test values match exactly, including 319125.
2. **Rank-zero DT from seven polar coefficients.** Construct the standard finite-dimensional modular ansatz from exact \(E_4,E_6,\eta^{-70},\vartheta^{(5)},D\vartheta,D^2\vartheta\) series. Seven input polar coefficients plus the two forced leading zeros produce an exact 9×9 rational matrix of rank 9. The resulting h0,h1,h2 nonpolar coefficients match the literature.
3. **DT→PT solver.** The genus-54 site corresponds to \(Q=21,m=-53\). Its wall-crossing charge is \(\chi=m+Q-5=-37\). The only admissible summand is the vacuum, so \(\mathrm{PT}(21,-53)=(-1)^{-37}(-37)\Omega_{1,1}(83/120)=37\cdot8625=319125\). This is the leading GV invariant because the previous top-genus positions vanish.

The same DT→PT solver computes the **first five PT coefficients** at degree 21 exactly, without loading the target PT table. It correctly includes degree-one PT wall contributions at the last two orders.

## Input provenance and holdout hygiene

- Empirical calibration: published leading PT/GV integers at degrees 6,7,8,9,10 (genera at most 16). Model-family validation: degrees 11–20. Hidden comparison set: degrees 21–30, all read only after fitting.
- Modular input: seven polar entries \(h_0:5,-800,58500;\ h_1:0,8625;\ h_2:0,0\). In the theta-series implementation, each of h1 and h2 has an additional forbidden leading power set to zero, yielding nine independent linear conditions. No nonpolar h_mu, degree-21 PT, or degree-21 GV target enters this matrix.
- Lower-degree wall states: \(\mathrm{PT}(1,0)=0,\ \mathrm{PT}(1,1)=2875,\ \mathrm{PT}(1,2)=-5750\). The fixed genus-51 degree-20 extremal input \(n_{51,20}=175\) is used **only** for the disconnected subtraction at the fifth PT/GV coefficient.
- Validation-only values come from Alexandrov et al., *Quantum geometry, stability and modularity*, arXiv:2301.08066v2, quintic Tables 2–3 and eq. (5.8). The code puts reference values in assertions after computing the predictions; they do not enter the solvers.

## 1. Empirical feature extraction and degree holdouts

Let \(m_C(Q)=-\lfloor Q(Q+5)/10\rfloor\), and for \(s=Q\bmod5\) take the first potentially nonzero row \(m=m_C(Q)+\delta_s\), where

| s | 0 | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|---:|
| delta_s | 0 | 1 | 2 | 2 | 1 |
| mu_s | 0 | 1 | 2 | 2 | 1 |

Define \(\chi=m+Q-5\). The empirically useful transformation is

\[\boxed{\frac{\mathrm{PT}(Q,m)}{(-1)^\chi\chi}=\Omega_{\mu_s}}.\]

The five degrees Q=6,…,10 give \(\Omega_0=5,\Omega_1=8625,\Omega_2=-1218500\), without reference to any higher-degree target. Predicting Q=11,…,30 then gives **20/20 exact integer matches**. Ten are genuinely beyond the g=53 frontier in the previously published modularity-conditional table. A rolling training-only scan over periods 1–15 selects 5 as the smallest exact period after charge normalization; period 10 is an equivalent nonminimal refinement.

On the **Q=21–30** evaluation: the zero, same-residue nearest-neighbor, signed quadratic, random forest, ExtraTrees, and KNN baselines each score **0/10 exact**. The normalized arithmetic model scores **10/10 exact**. These are small-sample outcomes; no claim is made that the transformation was discovered without the wall-crossing insight.

### Why the old phase/quadratic models fit

At \(Q=5k+1\), \(m=-\frac{5k^2+7k}{2}+1\), so

\[\chi= -\frac{5k^2-3k+6}{2}.\]

Therefore \(n_{g,Q}=8625(-1)^\chi\chi\) is precisely the prior quadratic envelope times the four-periodic Fourier sign. The other three nontrivial residues work similarly, with \(\Omega_2=-1218500\). The apparent empirically learned Fourier law is therefore a change of coordinates on the exact wall-crossing multiplier. This is a **post hoc explanation**; it is not an independent second prediction.

## 2. Reconstructing the rank-zero DT series from the polar data

The modular ansatz is

\[h_\mu=\eta^{-70}\left[A(E_4,E_6)\vartheta_\mu^{(5)}+B(E_4,E_6)D\vartheta_\mu^{(5)}+C(E_4,E_6)D^2\vartheta_\mu^{(5)}\right],\]

with the nine weight-correct basis monomials appearing in eq. (5.7) of the source paper. The code constructs the q-series using exact Bernoulli-divisor-sum Eisenstein series, the eta product, theta definitions and iterated Serre derivatives—**it does not import the published nine fitted coefficients**.

Solving the 9×9 matrix gives full rank, all rational coefficients, and the following independently computed terms (exponents shown relative to the natural theta base):

- \(h_0=q^{-55/24}(5-800q+58500q^2+5817125q^3+\cdots)\)
- \(h_1=q^{-323/120}(0+0q+8625q^2-1138500q^3+3777474000q^4+3102750380125q^5+\cdots)\)
- \(h_2=q^{-347/120}(0+0q+0q^2-1218500q^3+441969250q^4+\cdots)\)

The shifted \(h_1\) exponent is equivalent to the more conventional \(q^{-55/24+3/5}(0+8625q-1138500q^2+\cdots)\). All coefficients are exact integers after summation. Independent published nonpolar coefficients agree.

## 3. Solving degree 21 by exact wall crossing

Use the rank-zero DT/ PT formula

\[\mathrm{PT}(Q,m)=\sum_{Q',m'}(-1)^{\chi'}\chi'\mathrm{PT}(Q',m')\Omega_{1,Q-Q'}(\hat q_0'),\]

with \(\chi'=m-m'+Q+Q'-5\) and

\[\hat q_0'=m'-m-\frac{(Q'-Q)^2}{10}-\frac{Q+Q'}2+\frac{55}{24}.\]

The code enumerates *all admissible* \((Q',m')\) with the exact rational inequalities in eqs. (4.14)–(4.15). All five target sites satisfy \(\alpha>f(Q/5)\), the theorem's applicability condition. At \((Q,m)=(21,-53)\), the allowed Q' candidates are 1–4 but their m' ranges are empty; the result is literally a one-term calculation:

\[\boxed{319125=37(8625).}\]

The complete five-coefficient calculation gives:

| m | Predicted PT | Extra nonvacuum correction |
|---:|---:|---:|
| -53 | 319125 | 0 |
| -52 | 40986000 | 0 |
| -51 | 132211590000 | 0 (PT(1,0)=0) |
| -50 | -105493513413000 | -488750 |
| -49 | 19064998024136500 | -74922500 |

The final correction decomposes as \(-75,900,000+977,500\), from the two degree-one PT states. Every row matches published Table 3 *after* execution.

## 4. Converting the PT neighborhood to GV and checking disconnected pairs

At degree 21 the first four PT coefficients are triangular binomial combinations of the genus-54,…,51 GV invariants. For the fifth coefficient, the product of two lower-degree PT states also contributes:

\[\mathrm{PT}(20,-50)\,\mathrm{PT}(1,1)=175\times2875=503125.\]

Subtracting this term is essential. The reconstructed GV indices are

| genus | predicted GV(21) |
|---:|---:|
| 54 | 319125 |
| 53 | 7158750 |
| 52 | 129691149375 |
| 51 | -118821918509250 |
| 50 | 30276261813046500 |

Every entry matches published Table 2. Without the \(503125\) subtraction the fifth entry fails exactly by that amount. This is a valuable neighborhood consistency check, but PT and GV are related algebraically; it is not five independent physical measurements.

## 5. Exact BCOV ambiguity response and what it can/cannot verify

The predicted remaining genus-54 integer gives the homogeneous change

\[\boxed{\Delta a_{54}(X)=319125\,5^{107}X^{21}(1-X)^{32}.}\]

Restoring the one-dimensional section factor gives

\[\boxed{\Delta f_{54}(x)=319125\,5^{107}x^{21}.}\]

The exact MUM response is

\[\Delta F_{54}(q)=n_{54,21}q^{21}\left[1-3450q+4790625q^2-3685910000q^3+1751393550000q^4+\cdots\right].\]

This expansion has been independently calculated using the quintic Picard–Fuchs periods and mirror map over exact rationals. It verifies that the integer becomes the desired GV boundary value and obeys the expected triangular charge transport.

**Important:** the BCOV anomaly and ordinary conifold/orbifold boundary conditions alone do not select \(319125\). They admit this entire integer lattice. The independent selection here comes from the rank-zero DT input and the wall-crossing relation. We have *not* generated the complete physical genus-54 BCOV amplitude (which requires the genus-53 input polynomial); the BCOV check is a homogeneous-response certificate, not an independent full-amplitude verification. Alexandrov et al. report comparing modular predictions with direct integration through genus 64 in their 2024 paper.

## 6. Further out-of-table predictions

The exact DT-normalized leading-coefficient formula gives the following **unverified against external GV/PT tables** integers; the code separately verifies that the one-term theorem is applicable and that no nonvacuum pair is admissible in each case:

| Q | genus | predicted GV/PT |
|---:|---:|---:|
| 31 | 111 | -724500 |
| 32 | 117 | -108446500 |
| 33 | 124 | -115757500 |
| 34 | 132 | -879750 |
| 35 | 141 | -550 |
| 36 | 147 | 991875 |
| 37 | 154 | -147438500 |
| 38 | 162 | 155968000 |
| 39 | 171 | -1173000 |
| 40 | 181 | 725 |

These are conditional consequences of the three DT coefficients and wall-crossing theorem rather than claims of previously unknown mathematical invariants. Their provenance is distinct from the published Q≤30 holdout table.

## 7. How to reproduce

```bash
python -m pip install sympy numpy scikit-learn
python modular_core.py
python empirical_benchmark.py
python period_scan.py
python model_ablation.py
python pt_from_dt_solver.py
python pt_to_gv.py
python bcov_exact_response.py
python tests.py
```

All inputs, comparison-only labels, and tests are in plain-text Python, CSV, and JSON. No float is used in modular reconstruction, wall crossing, PT→GV, or BCOV response. Floats are used solely by the deliberately naive ML baselines.

## Sources and status

- S. Alexandrov, S. Feyzbakhsh, A. Klemm, B. Pioline, T. Schimannek, *Quantum geometry, stability and modularity*, arXiv:2301.08066v2 (2024), §§4–5, especially eqs. (4.12)–(4.15), (5.5)–(5.8), Tables 2–3: https://arxiv.org/html/2301.08066.
- M. Liu, Y. Ruan, *Castelnuovo bound and higher genus Gromov–Witten invariants of quintic 3-folds*, https://arxiv.org/abs/2210.13411.
- The wall-crossing theorem uses stability hypotheses, as discussed in the source; the extrapolation is conditional on those assumptions. The modular ansatz is physics-motivated, with independent confirmations in the source. The residue-class interpolation is retrospective and not a theorem of all-genus physics by itself.

## Supplement: stronger no-oracle protocol

For convenience, `prepare_blind_splits.py` physically writes three files: the degree-6–20 training values, the unlabeled list Q=21–40, and a separate reference-only Q=21–30 file. The stand-alone `blind_predict.py` reads **only the first two files** and writes the predictions with SHA-256 input hashes. `validate_blind.py` independently opens the third file and scores the held-out results. This achieves **10/10** exact integer predictions at Q=21–30, including 319125, with no target file opened by the predictor. The choice of charge-normalized feature family is retrospective and disclosed.

**Independent geometric spot checks at Q>30.** The d=5k extremal curve formula gives -550 at (Q,g)=(35,141) and +725 at (40,181), exactly matching the two out-of-table DT-based forecasts. These geometric checks do not yet validate the other eight Q31–40 forecasts.
