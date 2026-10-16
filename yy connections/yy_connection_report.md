# Mirror-quintic connection action on YY generators — exact matrices, local conifold connection, and P2/P3 cocycles

**Date:** 2026-10-16. **Status:** Exact genus-zero period monodromy and holomorphic-limit YY-ring action; numerical LCS-to-conifold Frobenius connection independently checked; induced substitution on the full published genus-two and genus-three YY polynomials exactly tested. This is **not** yet the full higher-genus wave-function / metaplectic polarization operator.

## 0. Coordinates and conventions

Use `s=5^5 z=3125 z`, `theta=s d/ds`, `D=psi d/dpsi=-5 theta`, `Delta=1-s`, and original YY `X=1-Delta^{-1}`. The Picard–Fuchs equation is

```
[theta^4 - s (theta+1/5)(theta+2/5)(theta+3/5)(theta+4/5)] f = 0.
```

The companion Gauss–Manin matrix acts on the matrix `A_ij=theta^{j-1} x_i` by

```
theta A = A * G(s),
G(s) = [[0,0,0,(24/625)*s/(1-s)],
        [1,0,0,( 2/  5)*s/(1-s)],
        [0,1,0,( 7/  5)*s/(1-s)],
        [0,0,1,2*s/(1-s)]].
```

The period basis `x=(x1,x2,x3,x4)` is Movasati's integral symplectic basis, rescaled by a **common constant** so `x2=omega0` and `x1/x2=psi1/(2 pi i psi0)` is the large-radius flat coordinate. This rescaling does not change monodromy matrices. Movasati, *Modular-type functions attached to mirror quintic Calabi–Yau varieties*, Math. Z. (2015), equations (1), period base change near eq. (1), and Theorems 1–2; https://arxiv.org/abs/1111.0357 ; https://archive.ymsc.tsinghua.edu.cn/pacm_download/188/8236-Modulartypefunctions.pdf.

## 1. Exact monodromies and symplectic pairing

```
M0 = [[1, 1, 0, 0],
      [0, 1, 0, 0],
      [5, 5, 1, 0],
      [0,-5,-1, 1]]
Mc = [[1, 0, 0, 0],
      [0, 1, 0, 1],
      [0, 0, 1, 0],
      [0, 0, 0, 1]]
Omega = [[ 0, 0, 1, 0],
         [ 0, 0, 0, 1],
         [-1, 0, 0, 0],
         [ 0,-1, 0, 0]].
```

All were checked exactly: `M0.T Omega M0 = Mc.T Omega Mc = Omega`, `(M0 Mc)^5=I`, `(Mc-I)^2=0`.

`Mc` sends `(x1,x2,x3,x4)` to `(x1,x2+x4,x3,x4)`. Therefore in special coordinates

```
tau0=x1/x2, tau1=x3/x2, tau2=x4/x2, tau3=d tau1/d tau0
```

we have `tau0' = tau0/(1+tau2)`, `tau1'=tau1/(1+tau2)`, `tau2'=tau2/(1+tau2)`. The symplectic identity `<x,Dx>=0` implies

```
D tau2 = tau1 D tau0 - tau0 D tau1.
```

## 2. Exact action on the exterior square

Define `Wij = xi D xj - xj D xi`; order `(W12,W13,W14,W23,W24,W34)`.

```
wedge2(M0) =
[[  1,  0, 0, 0, 0, 0],
 [  0,  1, 0, 1, 0, 0],
 [ -5, -1, 1,-1, 1, 0],
 [ -5,  0, 0, 1, 0, 0],
 [  0,  0, 0,-1, 1, 0],
 [-25, -5, 5, 0, 5, 1]]

wedge2(Mc) =
[[1,0,1,0,0, 0],
 [0,1,0,0,0, 0],
 [0,0,1,0,0, 0],
 [0,0,0,1,0,-1],
 [0,0,0,0,1, 0],
 [0,0,0,0,0, 1]].
```

The simultaneous fixed subspace of the six-dimensional representation is exactly `span(W13+W24)`. **But it vanishes on our periods**: Griffiths transversality gives `W13+W24=0`. Hence the five-dimensional physical Wronskian space, basis `(W12,W13,W14,W23,W34)`, carries

```
N0 = [[  1,  0, 0, 0, 0],
      [  0,  1, 0, 1, 0],
      [ -5, -2, 1,-1, 0],
      [ -5,  0, 0, 1, 0],
      [-25,-10, 5, 0, 1]]
Nc = [[1,0,1,0, 0],
      [0,1,0,0, 0],
      [0,0,1,0, 0],
      [0,0,0,1,-1],
      [0,0,0,0, 1]].
```

It has no nonzero common fixed vector. Both matrices preserve the quadric

```
Q = W12*W34 + W13^2 + W14*W23,
```

and the Plücker identity makes `Q=0` for actual period Wronskians. This provides an obstruction to finding a *linear scalar invariant* among the five Wronskians (not a no-go theorem for vector-valued blocks or the entire higher-genus amplitude).

## 3. Exact action on the holomorphic YY differential generators

Define `W=x2 D x1-x1 D x2 = -W12` and `Jn=D^n W / W`. The conifold period continuation multiplies W by

```
S = 1 + W14/W12
  = 1+tau2-tau0*tau1+tau0^2*tau3.
```

With `sigma_k=D^k(S)/S` and `sigma0=1`,

```
J1' = J1 + sigma1
J2' = J2 + 2 J1 sigma1 + sigma2
J3' = J3 + 3 J2 sigma1 + 3 J1 sigma2 + sigma3.
```

For every n, `Jn' = sum_{k=0}^n binomial(n,k) sigma_{n-k} Jk`. `X' = X`. The LCS monodromy fixes W and hence fixes every Jn.

The original YY generators are obtained by the inverse triangular change of variables

```
v1 = J1
v2 = (5 X J1 - 10 X - J2)/2
v3 = (-J3 + 5 X J2 + (25 X^2-35 X) J1 + 50 X -50 X^2)/2.
```

It follows that

```
delta v1 = sigma1

delta v2 = ((5 X -2 J1) sigma1 - sigma2)/2

delta v3 = [-3 J2 sigma1 -3 J1 sigma2 -sigma3
            +10 X J1 sigma1 +5 X sigma2 +(25 X^2-35 X) sigma1]/2.
```

The basis `Hi` defined earlier by `H1=J1`, `H2=J2+10X`, `H3=J3-50X+100X^2` transforms by the same substitution (the X shifts are fixed). Mixed monomials `X^k Jn` acquire only the displayed binomial connection terms.

The Gaussian/Bell-jet nature of the action is especially valuable: *unnormalized* `W Jn=D^n W` transform **linearly** under the constant exterior-square monodromy, while normalized `Jn` transform by the cocycle above.

## 4. Full P2 and P3 polynomials: calculated connection cocycles

Using the exactly transcribed polynomials of Yamaguchi–Yau (hep-th/0406078), the script `yy_connection_action.py` explicitly substituted the connection action into both full polynomials. Independent substitution after reconstructing `v1',v2',v3'` agreed exactly. The change is `delta Pg := Pg(J')-Pg(J)`.

For genus two the full result is

```
delta P2 =
 - (5/24)*sigma1^3
 + sigma1^2*(-115*J1/24+45*X/4+125/24)
 + sigma1*(-115*J1^2/24+14935*J1*X/36+875*J1/9
            -725*J2/12-2825*X^2/288-133085*X/144
            -25*sigma2/12-15625/288)
 + sigma2*(-725*J1/12+14125*X/72+3125/72)
 - (175/9)*sigma3.
```

The genus-three cocycle contains **22** nonzero monomials in the sigma variables (total sigma degree **6**). Its `sigma3^2` coefficient is `-8225/108`, while the coefficient linear in `sigma3` at zero other sigmas equals `partial P3 / partial J3`. Every term is written exactly in `yy_connection_p2_p3_cocycles.txt`, avoiding an unwieldy expanded expression in the report. The special-geometry continuation *does not* leave the raw P2/P3 polynomials invariant; it produces highly structured, explicitly computable, inhomogeneous terms.

The pure-X polynomial `f_g(X)` itself is unchanged as an algebraic polynomial, since `X'=X`. However, the LCS holomorphic normalization `omega0^{2g-2} Delta^{g-1}` changes by `(1+tau2)^{2g-2}` when the defining period changes from x2 to x2+x4. This must NOT be mistaken for the full higher-genus physical transformation law: BCOV wave-function polarization changes also generate Feynman/heat-operator corrections from lower genus. Those are not implemented here.

## 5. Numerical **LCS-to-conifold** Frobenius connection matrix

To compute the actual connection between local *solution bases*, not just their monodromies, choose `t=1-s` and solve the conifold Frobenius recurrences **over exact rationals**. The normalized basis is

```
f0 = 1 + (2/625)*t^3 + ...,
f1 = t + (7/10)*t^2 + (41/75)*t^3 + ...,
f2 = t^2 + (37/30)*t^3 + ...,
fL = (f1*log(t)+g(t))/(2*pi*i),
     g(t) = -(23/360)*t^3 + ... .
```

The log solution satisfies `fL -> fL+f1` on a positive conifold circuit. Solve for the constant matrix `C` in

```
(x1,x2,x3,x4)^T = C (f0,f1,f2,fL)^T.
```

Matching all four theta jets at `s=1/2` using 250 terms at 65-digit working precision gives (18 displayed significant digits)

```
C =
[[+1.29357398450410864 i,   +0.150766695123547301 i,  -0.0277921800168652449 i, 0],
 [ 1.07072586843015580,    -0.0247076138044847181,    0.0057845115995744471,  +0.355881271708588529 i],
 [-6.19501627714957489+3.23393496126027159 i,
  -1.01660471670258207+0.376916737808868252 i,
   0.140889979448830909-0.0694804500421631122 i, 0],
 [0, +0.355881271708588529 i, 0, 0]].
```

`T_c` is the identity except `(T_c)[3,1]=1` (zero-indexed), so `fL' = fL+f1`. The independently calculated connection satisfies

```
maxabs(C*T_c*C^{-1}-Mc) = 3.87781348e-64.
```

The sole log coefficient in x2 agrees with the vanishing-period coefficient in x4 to `2.37e-66`; all other forbidden log coefficients are < `1.6e-65`. A second 43-digit, 175-term run agreed, with monodromy error < `5.5e-42`. These tests pin down the actual analytic-connection normalization and the direction/sign of conifold monodromy.

## 6. Does the YY differential ring close after continuation?

This was tested independently of `P2`, `P3`. At `s=1/2`, derive period derivatives up to `D^5` directly from the Picard–Fuchs equation, compute `J1,...,J4` from both Wronskians, and test the nontrivial fourth-jet relation implied by the YY equation for `Dv3`:

```
J4 = -(25*J1^2*X^2 -70*J1^2*X +10*J1*J2*X
       -2*J1*J3 -100*J1*X^2-350*J1*X+J2^2
       +190*J2*X -30*J3*X+100*X^2+404*X)/2.
```

Errors: `8.35e-53` on the LCS branch and `3.76e-53` after conifold monodromy. This is a nontrivial numerical test that the augmented connection action is compatible with the original differential-ring closure; it is much stronger than merely checking the Leibniz rule for J2/J3.

## 7. Interpretation and next mathematical step

The no-scalar-invariant calculation suggests a useful correction to a naive Polyakov-block analogy: **we should seek equivariant, vector-valued blocks (or gluing data in a bundle), not scalar invariants in the bare Wronskians.** The promising objects are `X^k D^n W_ij`, transforming by constant matrices. The holomorphic YY polynomials are expressible using *ratios* of these sections, and the induced cocycles at genus 2 and genus 3 have just been computed.

What is still missing is the full BCOV change-of-polarization action on the *partition function*, including the quantum/propagator anomaly. In general `P_g(J') != P_g(J)` is **not** an obstruction to a crossing-symmetric partition function: the missing lower-genus graph terms must be included before interpreting the residual ambiguity. The next substantive calculation is to compute the symplectic-polarization Gaussian kernel for the same `Mc`, apply its expansion through `g=3`, and ask whether the genus-two and genus-three cocycles above are cancelled by the predicted graph terms. Only then can a residual, possibly simpler, holomorphic ambiguity be identified in a genuinely crossing-covariant basis.

## Reproducibility

- `yy_connection_action.py`: exact symplectic/wedge² matrices, YY ring substitution, P2/P3 action and assertions.
- `yy_connection_numeric.py`: exact rational conifold series, LCS Frobenius epsilon-series, 65-digit connection match, monodromy check.
- `yy_connection_ring_check.py`: independent fourth Bell-jet/ring-closure verification using PF-derived fifth derivatives.
- `yy_connection_p2_p3_cocycles.txt`: entire exact genus-three sigma-cocycle coefficient list (22 terms).
- Earlier `yy_p2_p3_exact.py` is the original full genus-two/genus-three polynomial module and is required beside the exact-action script.
