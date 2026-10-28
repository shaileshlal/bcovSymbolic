# A minimal-integer quintic bootstrap and periodic GV/PT arithmetic

**Status:** executed proof of concept with certified linearized BCOV propagation, and a forward-held-out arithmetic predictor. 27 October 2026.

## 1. Analytic baseline: exactly one exceptional input at genus 45–53

In the canonical YY gauge, after conifold, orbifold, and MUM constant-map constraints, the remaining genus-g ambiguity has `D_g=floor((2g-2)/5)` rational degrees of freedom. Represent these in the GV-charge lattice `n[g,d]`, `1<=d<=D_g`, and choose minimum height (zero) as the initial point. The proven quintic vanishing `n[g,d]=0` for `10*g > d^2 + 5*d + 10` fixes 168 of the 169 lattice sites at g=45..53. The sole exception is (51,20).

For `d=5k`, the Liu–Ruan projective-bundle formula yields the extremal value

`N_k = binomial(k+3,3) - binomial(k-2,3)`,

`n[ 1+5*k*(k+1)/2, 5*k ] = (-1)**(N_k+3) * 5*N_k`, for k>=2.

At k=4, `N_4=35`, so `n[51,20]=+175` exactly. Prior examples `n[16,10]=-50`, `n[31,15]=-100` are consistent. Simple linear extrapolation from those prior two values would incorrectly predict `-150`, illustrating why the geometric/parity feature matters.

The full integer charge vector at g45..53 is thus `0,0,...,0` with one nonzero integer +175. This is an **analytic/geometric no-target-oracle construction**, not a neural predictor. Its available published charge benchmark is exact at 169/169 sites, but the complete physical YY polynomials through genus 53 are not independently loaded or recomputed in this release.

## 2. Exact BCOV sensitivity g51 → g52 → g53

Set `x=psi^-5` and `X=x/(x-1)`. The reduced section `f_g=(x-1)^(g-1)*a_g(X)` changes under the exceptional integer by

`delta f_51=175*5**100*x**20`.

With `C=175*5**100` the pure-X changes are

`delta a_51 = C*X**20*(1-X)**30`,

`delta a_52 = C*X**20*(1-X)**30*(15625/6 - 46675*X/2)`,

`delta a_53 = C*X**20*(1-X)**30*(162109375/48 - 4573065625*X/72 + 40994061325*X**2/144 - 46375*X**3/18 - 3125*X**4/3)`.

The exact nonholomorphic change is also propagated. At genus 52 the source difference is

`dR52 = -1/2 * (D51^2 dP51 + 2*P1^(1)*D51 dP51)`.

At genus 53 it is

`dR53 = -1/2 * (D52^2 dP52 + 2*P1^(1)*D52 dP52 + 2*P2^(1)*D51 dP51)`.

Here `Dg` denotes the YY covariant insertion operator, not the floor-index `D_g` above. The integration in `v1,v2,v3` uses the canonical zero-pure-X primitive. The additional pure-X completions are fixed from the homogeneous conifold/orbifold and MUM/GV conditions for these *differences*. The script verifies exact BCOV anomaly residuals and all three boundary conditions at genera 52 and 53; it generates exact serializable coefficient lists with numerator and denominator strings.

Terms in the nonholomorphic primitive: 157 at genus 52, 679 at genus 53. These are exact linearized (in the exceptional input) changes, not physical full P_g amplitudes.

## 3. Near-Castelnuovo EDA and held-out integer ML

Read `g_C(d)=floor((d*d+5*d+10)/10)`. For `d=5k+s`, s=1,2,3,4, the first possibly nonzero charge occurs at `g=g_C(d)-delta_s` where `delta_1=delta_4=1` and `delta_2=delta_3=2` in the studied table.

The explicit integer-phase laws recovered from training on *only k=1,2,3* in each class are:

| Degree class | Genus index | Integer-valued extrapolation |
|---|---|---|
| `d=5k+1` | `g_C(d)-1` | `8625/2 * (5k^2-3k+6) * (cos(pi k/2)-sin(pi k/2))` |
| `d=5k+2` | `g_C(d)-2` | `609250 * (5k^2-k+4) * (cos(pi k/2)+sin(pi k/2))` |
| `d=5k+3` | `g_C(d)-2` | `609250 * (5k^2+k+4) * (cos(pi k/2)-sin(pi k/2))` |
| `d=5k+4` | `g_C(d)-1` | `8625/2 * (5k^2+3k+6) * (cos(pi k/2)+sin(pi k/2))` |

Because k is an integer, each trigonometric character takes values 0 or ±1, and the predicted results are integers. In the implementation, these phases are computed with a modular lookup, never floating point. The quadratic magnitude coefficients are inferred from the three training magnitudes by exact integer finite differences, and the phase sine/cosine coefficients from the signs of the first two training points. Their agreement with the third training sign is explicitly tested.

Crucially all 12 training labels lie at genus ≤53 (e.g. d=6,g=6; d=11,g=17; d=16,g=33), while all eight *held-out* labels are at genus ≥54. The holdouts are at d=21,22,23,24 and 26,27,28,29, respectively g=54,58,63,69 and 80,85,91,98. All eight are recovered exactly from the first three earlier points per class (no polynomial is fitted to a heldout integer).

Baselines: minimum height zero, nearest-genus label, and unconstrained quadratic interpolation of the *signed* values each give zero exact matches. The phase-normalized model scores 8/8. This is an arithmetic model with a physically motivated **preselected period-four Fourier dictionary**, not a claim that period four was uniquely discovered using so few examples.

The eight high-genus reference numbers originate from Table 2 of arXiv:2301.08066v2 and are **modularity conditional**. The matching formulas may have an exact geometric explanation involving near-extremal curve moduli, but this has not been proved in this package. Nor should agreement between published PT and GV entries be counted as independent evidence, because the PT/GV relation is an identity.

The frontier result is a candidate `n[54,21]=319125` and hence **the residual** (not the complete F54)

`delta f_54 = 319125 * 5**107 * x**21`.

## 4. PT/DT: which periodicity is real and where new constraints enter

For a primitive GV charge at genus g and degree d, its contribution to the rank-one stable-pair partition function at degree d begins with

`delta[Q**d]Z_PT(q) = n[g,d]*q**(1-g)*(1+q)**(2g-2)`.

At d=20, the known `n[51,20]=175` contributes `175*q^-50*(1+q)^100`. The next PT coefficient additionally includes `n[50,20]=9700`, so

`PT(20,-49) = 100*175 + 9700 = 27200`, matching the published row.

At d=21, the near-extremal candidate `n[54,21]=319125` gives leading `PT(21,-53)=319125`. The next entry includes `n[53,21]=7158750`:

`PT(21,-52) = 106*319125 + 7158750 = 40986000`, again matching the published modularity-conditional PT Table 3.

The rank-one ideal-sheaf DT function is related by the universal degree-zero MacMahon factor

`Z_DT_rank1(Q,q) = M(-q)^chi * Z_PT(Q,q)` with `chi=-200`,

`M(-q)^(-200) = 1 + 200*q + 19500*q**2 + ...`.

For the isolated n51 change this gives `delta DT_rank1` coefficient sequence `175,52500,7778750` at powers q^-50,q^-49,q^-48. These integer relations do **not** independently determine the GV lattice point.

The different **rank-zero** DT (D4–D2–D0) spectral-flow charges `mu in Z/5Z`, modular weight -3/2 and three independent components `h0,h1,h2` are the genuine additional structure. The 2024 source constructs the associated vector-valued modular form from its seven polar coefficients and uses wall crossing to constrain/predict the near-boundary GV/PT tables. The `d mod 5` and Fourier `k mod 4` pattern discovered above should not be confused with an equality `mu=d mod 5`: the mapping passes through charge and spectral-flow conventions.

## 5. Reproducibility and nonclaims

- Every inferred model coefficient and GV integer is calculated with Python integers; no float rounding is used anywhere in fitting.
- The conifold/orbifold/MUM projection and the high-genus sensitivity polynomial are exact rational calculations (run `python code/bcov_sensitivity_51_53.py`).
- `python code/tests.py` reloads the serialized YY polynomials and checks both anomaly equations from the sources independently.
- The holdout references are published *modularity-conditional* numbers; this is a powerful forecasting and consistency test, but does not prove our near-extremal formula for arbitrary degree.
- Full physical P_g for g=45..53 and an independent BCOV P54 have not been recomputed or coefficient-by-coefficient compared.
- There is no evidence here that rank-one PT and DT integrality add new constraints beyond GV integrality.

## Sources

1. Z. Liu and Y. Ruan, *Castelnuovo bound and higher genus Gromov-Witten invariants of quintic 3-folds*, https://arxiv.org/abs/2210.13411.
2. S. Alexandrov, S. Feyzbakhsh, A. Klemm, B. Pioline, T. Schimannek, *Quantum geometry, stability and modularity*, https://arxiv.org/html/2301.08066 (version 2, 24 January 2024; §§2.6,3,5.2, Tables 2–3).
