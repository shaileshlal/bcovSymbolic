# Quintic GV-lattice displacement prediction: exploratory study v0.3

**October 2026.** Exact BCOV calibration, geometry-gated dataset to genus 1000, arithmetic-aware model search and integer neighbor transport.

## Executive conclusion

This study targets **integer GV-charge displacements** from the unconstrained-lattice minimum-height origin, **not the enormous rational YY ambiguity coefficients**. The first exact reference dataset (g=4–12) contains only zeros. A genuine high-genus training study therefore requires independently known nonzero enumerative inputs. We extracted a provable extremal family from the theorem of Liu–Ruan. A compact phase-plus-quadratic symbolic model, fitted on the first eight extremal cases, predicts the next eight **exactly**, whereas naive black-box regressors fail to produce any exact integers. This is recovery of a known formula, not a computation of unknown physics. The g=54,d=21 charge remains **unpredicted**.

An additional **exact new organizational identity** simplifies the first nonzero lattice displacement and its next-degree consequence; it was verified against all nine genus-4–12 transport calculations.

## 1. Define the target without data leakage

For the mirror quintic in the original YY gauge, set

\[
 x=\psi^{-5},\quad X={x\over x-1},\quad
 D_g=\lfloor(2g-2)/5\rfloor,\quad m_g=\lceil(3g-3)/5\rceil.
\]

After the conifold gap, orbifold regularity, and MUM constant map,

\[
 \delta a_g(X)=X(1-X)^{m_g}\sum_{r=0}^{D_g-1}c_{g,r}X^r,
\qquad
 \delta f_g(x)=\sum_{r=0}^{D_g-1}(-1)^{g-2-r}c_{g,r}
 x^{r+1}(1-x)^{D_g-r-1}.
\]

Let \(n_{g,d}\), \(d=1,\dots,D_g\), be the GV charges obtained after the MUM mirror transform. In the **unrestricted integer GV lattice**, the minimum-height integer representative is 0, so our integer target is

\[
 \Delta n_{g,d}=n_{g,d}-0=n_{g,d}.
\]

This equivalence depends on how we define the lattice. If a *nontrivial denominator cap or congruence sublattice* is imposed, its minimal representative may instead be nonzero. Such a cap must be specified without using the target denominators. We do not use the retrospectively fitted denominator caps in prediction.

The historical reference audit has 22 residual coefficients across g=4–12. All 22 have \(\Delta n=0\) **by established geometric vanishing**, not by ML. Fourteen were tested with an independently solved boundary particular; the remaining eight used a reference-derived particular. Their rational and lattice checks are bundled for reproducibility.

## 2. Build a three-way, censored training set

Use the quintic genus bound \(G(d)=(d^2+5d+10)/10\). A point with \(g>G(d)\) is labelled **certified_zero**. The extremal case \(d=5k\) at \(g=1+5k(k+1)/2\) has a proven formula for k≥2 and is labelled **extremal_nonzero** when it lies in the residual space. Every other site remains **unresolved**. Unknown values are represented by empty cells, never as 0.

| Genus scope | Certified zeros | Proven extremal nonzeros | Unresolved |
|---|---:|---:|---:|
| 4–12 | 22 | 0 | 0 |
| 4–53 | 529 | 1 | 0 |
| 4–100 | 1,683 | 2 | 255 |
| 4–1000 | **63,544** | **16** | **135,840** |

This is **199,400** sites in total, of which **68.1% remain unresolved**. Thus it would be an error to train a generic classifier by filling the missing labels with zero. In particular, a 99.9% apparent accuracy on proven zeros would not establish ability to select nonzero displacements.

The first unprotected cell is **(g,d,r)=(54,21,20)**, followed at g55 by the same degree. The previously known triplet is \(n_{51,20}=175\), \(n_{52,20}=n_{53,20}=0\). The next source-complete extremal point in the residual sector is \((g,d)=(76,25)\), with \(n_{76,25}=275\).

### Geometry-aware features

The useful purely structural features are \(g\), \(d=r+1\), \(D_g\), \(g\bmod5\), \(r\bmod5\), \(d\bmod5\), and, importantly, the exact signed margin

\[
 M_{g,d}=d^2+5d+10-10g.
\]

A hard gate imposes \(M<0\Rightarrow\Delta n=0\). On the known top family, \(M=0\), \(g\equiv1\pmod5\), \(r\equiv4\pmod5\), and \(k=d/5\) becomes a meaningful discrete position. **These conditions do not predict generic unresolved cells with M>0.**

## 3. Discover the extremal periodic law using a training holdout

For \(d=5k\), the known extremal theorem (not a fitted assumption) supplies

\[
  n_{1+5k(k+1)/2,\,5k}
  =(-1)^{k(k-1)/2}\,25\left(1+{k(k-1)\over2}\right),\quad k\ge2.
\]

Within our residual sector, the examples begin at k=4:

| k | genus | d | exact displacement | Split |
|---:|---:|---:|---:|---|
| 4 | 51 | 20 | 175 | train |
| 5 | 76 | 25 | 275 | train |
| 6 | 106 | 30 | −400 | train |
| 7 | 141 | 35 | −550 | train |
| 8 | 181 | 40 | 725 | train |
| 9 | 226 | 45 | 925 | train |
| 10 | 276 | 50 | −1150 | train |
| 11 | 331 | 55 | −1400 | train |
| 12–19 | 391–951 | 60–95 | **held out** | test |

Train only k=4,…,11; predict k=12,…,19. A small finite search over sign periods and exact rational polynomials selects the *shortest repeated sign period* p=4 and the *lowest fitting magnitude degree* 2. Both patterns repeat over training examples; the period is not an arbitrary high-period interpolation. The learned formula is

\[
 \widehat{\Delta n}_k
  =25\left(1+{k(k-1)\over 2}\right)s(k),
 \quad s(k)=\begin{cases}+1 & k\equiv0,1\pmod4\\-1 & k\equiv2,3\pmod4.\end{cases}
\]

This admits an **exact finite-Fourier / RoPE-like realization on integer positions**:

\[
 s(k)=\cos(\pi k/2)+\sin(\pi k/2).
\]

The phase is not merely a cosmetic embedding here; it is precisely the sign pattern of this *special extremal family*. No phase advantage is established for general interior GV invariants.

### Holdout benchmark, eight withheld integers

| Model | MAE in integer charge | Exact integers / 8 |
|---|---:|---:|
| Minimum-height zero | 2900.00 | 0 |
| Raw ridge quadratic (no phase) | 5561.90 | 0 |
| Raw 3-nearest-neighbors | 2900.00 | 0 |
| Raw random forest | 2900.00 | 0 |
| Raw ExtraTrees | 2900.00 | 0 |
| **Period-four sign + degree-two magnitude (exact rational solver)** | **0** | **8** |
| **Two-dimensional rotary phase × quadratic envelope** | ~5×10⁻¹² (float) | **8 after integer rounding** |

The phase period was also chosen in a **pre-holdout** search, training on k=4–9 and validating on k=10–11. Period 4 with quadratic envelope had validation MAE below 10⁻¹¹, beating alternate tested periods. Then the training set was expanded to k=4–11, with the eight late examples untouched until scoring.

Control: among 70 balanced random sign sequences of length eight, six have an exact period-four repetition. Among all 40,320 permutations of the eight observed positive magnitudes, only two are exactly quadratic. These are descriptive permutation controls, **not blind hypothesis-significance tests**—we are rediscovering an already known analytic formula.

## 4. A stronger *neighborhood regulator* from the exact BCOV map

The previous exact BCOV lattice map is triangular:

\[
 n^{\rm phys}_{g,d}=n^{\rm anchor}_{g,d}+\sum_r H_{g;d,r}c_{g,r},
 \qquad H_{g;d,d-1}=(-1)^{g-1+d}5^{5d-4g+4}.
\]

If all charges before degree d vanish, the leading affected Taylor coefficient in the reduced section is

\[
 \boxed{[x^d]\Delta f_g(x)=5^{4g-4-5d}\,\Delta n_{g,d}.}
\]

This statement holds for the first affected coefficient in the section, even if d is not the last free index. For the *last free degree* \(d=D_g\), it is the entire displacement:

\[
 \boxed{\Delta f_g(x)=5^{4g-4-5D_g}\Delta n_{g,D_g}\,x^{D_g}.}
\]

The exact historical H matrices verify these formulas in every one of their 22 unit-charge directions.

The MUM amplitude response for a last-charge perturbation simplifies dramatically:

\[
 \boxed{\Delta F_g^{\rm MUM}(q)
  =\Delta n_{g,D_g}\,\omega_0(z(q))^{2g-2}z(q)^{D_g},\qquad z=x/3125.}
\]

At low order, \(\omega_0(z)=1+120z+\cdots\), and \(q=z+770z^2+\cdots\). Therefore the next-degree coefficient obeys a universal **linear nearest-neighbor transport law** for D≥2:

\[
 \boxed{\Delta n_{g,D_g+1}
  =[240(g-1)-770D_g] \Delta n_{g,D_g}.}
\]

If D=1, subtract the degree-one twice-covered GV contribution \(2^{2g-3}\Delta n_{g,1}\). We reproduced the previously recorded exact transport entries in **9/9 genera (4–12)**, including the D=1 exceptions. This is an exact structural relation between the effect of two different ambiguity choices, not a new equation selecting the physical charge.

### Concrete lattice transitions

At genus 51, the only free charge not fixed to vanish is n_(51,20)=175. Therefore

\[
 \boxed{\Delta f_{51}(x)=175\,5^{100}x^{20},\quad
    \Delta n_{51,21}=-3400\cdot175=-595000.}
\]

At genus 54, lower charges n_(54,1…20) vanish and the remaining integer is unknown:

\[
 \boxed{\Delta f_{54}(x)=5^{107}\,n_{54,21}\,x^{21},\quad
    \Delta n_{54,22}=-3450\,n_{54,21}.}
\]

These formulas are *exact* even though the value of n_(54,21) is not available from the labelled training data. The naive minimum-height candidate n=0 is a baseline, **not a validated physical forecast**.

## 5. What is genuinely learned, and what is not

1. **Established by exact rational computation:** the 22 genus-4–12 calibration directions, the charge-to-section leading coefficient, the last-direction simplification, and the nine archived neighbor responses.
2. **Recovered from analytic labels by ML/symbolic inference:** the special extremal family’s period-four sign, quadratic amplitude, and matching Fourier embedding (8/8 exact heldout labels).
3. **Unresolved:** the integers at g=54,d=21 and the 135,839 other unknown sites up to genus 1000; there is no supervised evidence for a genus/degree neighborhood law inside that region. The exact periodic extremal formula should *not* be extrapolated by substituting fractional k=d/5.

### Recommended continuation

Obtain independent near-boundary nonzero GV labels beyond the extremal d=5k line (MSW polar data, DT/PT wall crossing, or published BCOV reconstructions). Then train a **gated integer model** on \((g,d,M_{g,d},d\bmod5,g\bmod5,\text{nearby known charges},\text{YY-source features})\), predicting sign, \(\log(1+|n|)\), and a small integer correction separately. Use an exact solver to translate the proposed integer vector by \(H_g^{-1}\), and only score **forward-genus, completely hidden enumerative inputs**. Unknown cells must remain censored; calibration and abstention are preferable to false precision.

For g=54,d=21, the next meaningful computation is a genuinely independent boundary value or modular/BPS constraint selecting n_(54,21). The structural formulas above then give the complete polynomial displacement and a nontrivial next-degree consistency check automatically.

## Source provenance

- Zhiyu Liu and Yongbin Ruan, *Castelnuovo bound and higher genus Gromov–Witten invariants of quintic 3-folds*, https://arxiv.org/abs/2210.13411 (geometric zero and extremal labels).
- *Topological String Theory on Compact Calabi–Yau: Modularity and Boundary Conditions*, https://arxiv.org/abs/hep-th/0612125 (BCOV boundary framework).
- Previous exact YY/BCOV residual audit and integer transport: `source/exact_lattice_audit.json` and `source/integral_transport_degree26.json` (bundled from the earlier conversation's released experiment).

**Reproducibility:** `python build_data.py && python eda.py && python learn_top.py && python permutation_control.py && python exact_solver.py && python neighbor_regulator.py && python frontier_solver.py && python tests.py`. `displacement_marimo.py` provides an optional notebook UI. All physical labels have an explicit provenance string; missing labels are never imputed as zero.
