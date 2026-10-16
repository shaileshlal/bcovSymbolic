# Towards Polyakov-type blocks for BCOV on the mirror quintic

**Status:** proven finite-dimensional homogeneous gap/LCS results, with an exact genus-two and genus-three computation. **Not yet:** the full conifold connection coefficients, metaplectic symplectic transport of the BCOV partition function, or an independent modular D4-D2-D0 constraint engine.

## 1. The crossing analogy and the relevant group

Put `s=5^5 z`. The Picard–Fuchs operator is

\[
\mathcal L=\theta^4-s(\theta+1/5)(\theta+2/5)(\theta+3/5)(\theta+4/5),\quad \theta=s\partial_s.
\]

At `s=0,1,infty` the respective exponents are `(0,0,0,0)`, `(0,1,1,2)`, `(1/5,2/5,3/5,4/5)`. Continuation between MUM, conifold, and Gepner frames is analogous to changing OPE channels, but the involution `s -> 1-s` is **not** a symmetry of this differential equation. One must compute connection matrices between local Frobenius bases, and for higher genus also the symplectic polarization change of the full BCOV partition function.

Related groundwork: Movasati, *Modular-type functions attached to mirror quintic Calabi–Yau varieties*, arXiv:1111.0357 (a Ramanujan-type differential system and monodromy functional equations); Alim–Scheidegger–Yau–Zhou, *Special Polynomial Rings...*, arXiv:1306.0002; Huang–Klemm–Quackenbush, arXiv:hep-th/0612125.

## 2. Homogeneous conifold gap = a distinguished subspace

Use original YY variable `X=1/(1-psi^5)`, discriminant `Delta=1-3125 z`, and `X=1-Delta^-1`. The allowed pure-X ambiguity space is

\[
 \mathcal A_g=\mathbb Q[X]_{\le 3g-3},\quad \dim \mathcal A_g=3g-2.
\]

In a holomorphic conifold frame the ambiguity contribution has form

\[
 \delta F_g^c=\Delta^{g-1}U_g(\Delta)f_g(X),\quad U_g(0)\ne 0.
\]

Then the pole-free **homogeneous** conifold kernel is exactly

\[
\mathcal K_g=\operatorname{span}\{(1-X)^r:0\le r<g\},\quad \dim \mathcal K_g=g.
\]

This is an exact order-of-pole statement and is independent of the coefficients of analytic nonvanishing `U_g`. A choice of flat conifold coordinate does not change which deformations are pole-free, because `t_c=alpha Delta+O(Delta^2)` is a local uniformizer.

Conifold boundary constraints have rank `2g-2`. The *inhomogeneous* prescribed gap determines a coset of this kernel after adding the BCOV particular solution; it does **not** assert that the gap-polar coefficient of the YY polynomial's pure-X projection vanishes.

Genus 2: `dim A_2=4`, rank 2, kernel 2.

Genus 3: `dim A_3=7`, rank 4, kernel 3.

The trivial single-locus LCS filtration has one-dimensional quotients because `X=-3125q+O(q^2)` is already a uniformizer. The genuinely joint filtration is `K_g^(d)= {f in K_g: [q^e] E_g f=0 (0<=e<d)}`, where

\[
 E_g=\varpi_0(z)^{2g-2}\Delta^{g-1},\qquad q=z+O(z^2).
\]

It has dimensions `g,g-1,...,0`. Its successive quotients are one-dimensional as well, **but only after the nontrivial conifold restriction**.

## 3. Exact MUM + conifold dual bases

Define `T[g]_{dr}=[q^d]E_g(1-X)^r` for `d,r=0,...,g-1`. This matrix is invertible, with determinant `3125^(g(g-1)/2)` for the present normalized LCS coordinate (verified at g=2 and g=3; follows in general by leading-order Vandermonde/binomial evaluation). Define

\[
K_{g,d}(X)=\sum_{r=0}^{g-1}(T[g]^{-1})_{rd}(1-X)^r
=\left[\frac{q(X)^d}{E_g(X)}\right]_{X^{<g}}.
\]

Then `K[g,d]` lies in the homogeneous conifold gap kernel and `[q^e] E_g K[g,d] = delta_{ed}` for `0<=e,d<g`.

### Genus two

`T[2]=[[1,1],[-2885,240]]`, determinant `3125`.

\[
 K_{2,0}=1-\frac{577}{625}X,\quad K_{2,1}=-\frac{X}{3125}.
\]

### Genus three

\[
T[3]=\begin{pmatrix}
1&1&1\\ -5770&-2645&480\\11748525&1076650&170400
\end{pmatrix},\quad\det T[3]=30517578125=5^{15}.
\]

\[
\begin{aligned}
K_{3,0}&=1-\frac{1154}{625}X+\frac{318241}{390625}X^2,\\
K_{3,1}&=-\frac{X}{3125}+\frac{683}{1953125}X^2,\\
K_{3,2}&=\frac{X^2}{9765625}.
\end{aligned}
\]

## 4. GV-dual and modular-compatible residual basis

The primitive GV response at genus g is obtained from the Fourier coefficients by the multiple-cover Möbius transform. Within the gap kernel, define

\[
 C_{g,0}=K_{g,0},\qquad
 C_{g,d}=\sum_{\substack{r\ge 1\\rd<g}}r^{2g-3}K_{g,rd},\quad 1\le d<g.
\]

Then a shift of `P_g` by `C_{g,d}` changes the constant-map Fourier coefficient by zero (`d>0`) and changes `n_e^g` by `N_g delta_{ed}` for `1<=e<g`, holding lower genera fixed. Here `N_g=(-1/625)^(g-1)` is the normalization carried by the YY convention.

The concrete bases are

\[
C_{2,0}=1-\frac{577}{625}X,\qquad C_{2,1}=-\frac{X}{3125};
\]

\[
\begin{aligned}
C_{3,0}&=1-\frac{1154}{625}X+\frac{318241}{390625}X^2,\\
C_{3,1}&=-\frac{X}{3125}+\frac{3423}{9765625}X^2,\\
C_{3,2}&=\frac{X^2}{9765625}.
\end{aligned}
\]

These are simultaneous *homogeneous gap-preserving and GV-dual* modes. They are not presently known to obey a scalar automorphy law; the name "crossing-adapted" refers specifically to their diagonal boundary data.

## 5. The complementary conifold-polar modes

The leading-discriminant basis is `(1-X)^(g-1+j)`, `j=1,...,2g-2`, because its conifold amplitude is `Delta^-j` when `U_g=1`, `t_c=Delta`. To diagonalize the **exact** gap in physical flat coordinate, replace these by

\[
 B_{g,j}^c=\operatorname{PP}_\Delta\left[\Delta^{-(g-1)}U_g(\Delta)^{-1}t_c(\Delta)^{-j}\right],
\]

where `PP_Delta` keeps negative powers of `Delta`. Then

\[
 \Delta^{g-1}U_g B_{g,j}^c=t_c^{-j}+O(t_c^0).
\]

A genus-g ambiguity is uniquely decomposed into `2g-2` gap-dual polar modes and `g` gap-regular `C_{g,d}`. This decomposition is **exact once** `U_g` and `t_c` have been obtained from period continuation.

For genus two, if `t_c=Delta(1+c1 Delta+c2 Delta^2+...)` and `U_2=1+u1 Delta+u2 Delta^2+...`, the two pole-dual modes are

\[
\begin{aligned}
B_{2,1}^c&=\Delta^{-2}-(u_1+c_1)\Delta^{-1},\\
B_{2,2}^c&=\Delta^{-3}-(u_1+2c_1)\Delta^{-2}\\
&\quad +(u_1^2-u_2+2u_1c_1+3c_1^2-2c_2)\Delta^{-1}.
\end{aligned}
\]

This is a concrete recipe for the physically correct "gap diagonalization", not the (in general incorrect) identification of `Delta` with the flat conifold coordinate.

## 6. Application to YY genus two and three polynomials

For a fixed **chosen** particular polynomial, a simple bookkeeping convention is its pure-X projection `f_g^(0)(X)=P_g(0,0,0,X)`. Its conifold Laurent-basis decomposition is informative but must NOT be mistaken for a decomposition of the full physical conifold free energy, because the `v_i` pieces also contribute at the conifold.

Genus two:

\[
f_2^{(0)}=\sum_{r=0}^3 c_r(1-X)^r,
\qquad(c_0,c_1,c_2,c_3)=\left(-\frac{12695}{72},\frac{1760}{9},\frac{215}{72},-\frac{25}{48}\right).
\]

Its gap-regular part is

\[
 f_{2,\mathrm{reg}}^{(0)}=-\frac{5(2816X-277)}{72}
 = \frac{1385}{72}C_{2,0}+\frac{40004275}{72}C_{2,1}.
\]

Genus three:

\[
(c_0,\ldots,c_6)=\left(-\frac{7020115}{1296},\frac{10561435}{864},-\frac{211210255}{36288},-\frac{757655}{567},\frac{1035625}{2016},-\frac{18125}{144},\frac{15625}{1008}\right).
\]

Its gap-regular part is

\[
 f_{3,\mathrm{reg}}^{(0)}=-\frac{5(42242051X^2+4231952X-7161359)}{36288},
\]

and in the GV-dual `C` coordinates it has coefficients

\[
\left(\frac{35806795}{36288},-\frac{70240478575}{18144},-\frac{233326698293725}{4536}\right).
\]

The coefficients are still large because this is a convention-dependent projection, not the physical fully gap-subtracted ambiguity. The true objective is to use full period transport to determine `B^c` and the BCOV particular solution and **then** test whether the remaining `C` coefficients become simpler.

## 7. Pullback of modular/PT constraints

Write the residual homogeneous ambiguity as `delta P_2 = b20*C20+b21*C21`, `delta P_3 = b30*C30+b31*C31+b32*C32`. In these coordinates the exact degree-four GV response is

\[
\begin{aligned}
\delta n_4^2&=-\frac{1}{625}\left(42804082800 b_{20}-50613940b_{21}\right),\\
\delta n_4^3&=\frac{1}{390625}\left(-15499620000 b_{30}+72523636b_{31}+367142b_{32}\right).
\end{aligned}
\]

The quintic modular/PT boundary data discussed previously yield `PT(4,-2)=n_4^3` and `PT(4,-1)=n_4^2+4n_4^3`. Their *homogeneous variations* are hence two independent linear constraints: one from `delta n_4^3`, one from `delta n_4^2+4 delta n_4^3`. Constant-map conditions fix `b20,b30`. The remaining three unknowns `b21,b31,b32` face two degree-four equations, leaving at least one degree of freedom until an additional PT/GV, orbifold or modular constraint is supplied. This is a rank statement, not proof that these inputs alone determine the entire ambiguity.

Inserting numerical PT target values requires subtracting the BCOV particular solution and the conifold gap-fixed part first. This subtraction has NOT been executed in the present two-boundary script.

## 8. Next test for genuine automorphy

The full meaningful symmetry is a gluing problem of local BCOV amplitudes: conifold, LCS and orbifold local frames are related by period connection matrices and (for the partition function) a symplectic change of polarization. A universal `Polyakov-BCOV` basis would have a simple connection law under those transitions, with only genus-lower anomaly data as prescribed inhomogeneities.

Movasati's seven differential generators / Ramanujan-type system, arXiv:1111.0357, give an independently motivated candidate automorphic enhancement of the YY generators. The actual conifold period matrix also has arithmetic connections to periods and quasiperiods of a weight-four form of level 25 (see work discussed by Scheidegger/Klemm). Neither fact implies that the `C[g,d]` themselves are classical modular forms.

A decisive experiment is to transport `C[g,d]` to the conifold and Gepner frames via the full fourth-order PF system and BCOV symplectic wave-function action. If the transformed blocks have uniformly sparse decomposition coefficients at g=2,3,4,..., with a genus-independent recursion, the Polyakov analogy would have acquired content beyond a finite-dimensional basis choice.

## Reproducibility

Run `python yy_crossing_gap_basis.py` next to `yy_p2_p3_exact.py`. SymPy only. `yy_crossing_gap_results.txt` contains the complete exact matrices and checks. All reported genus-two and genus-three identities and Möbius tests are exact rational computations.
