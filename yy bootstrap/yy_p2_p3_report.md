# Exact YY basis transformations for the quintic: genus 2 and genus 3

## Input and conventions

This is an exact rational computation from Yamaguchi--Yau, *Topological String Partition Functions as Polynomials* (2004), eqs. (4.1), (4.8), (A.1); see https://arxiv.org/pdf/hep-th/0406078 . No modular-bootstrap constraints are assumed.

\[
D=\psi\partial_\psi=-5z\partial_z,\quad z=(5\psi)^{-5},\quad
X=(1-\psi^5)^{-1}=-\frac{3125z}{1-3125z},\quad P=1-X=\frac{1}{1-3125z}.
\]
At large complex structure (MUM), \(X=0\); at the conifold, \(X\) has a pole. The earlier proposal to use \((X-1)^k\) as a MUM-triangular basis applied to the *inverse-discriminant* generator \(P\), not the original YY variable \(X\). For the original YY variable, \(X^k=O(q^k)\).

The ring is
\[
Dv_1=-v_1^2-2v_2-10X+5v_1X,\quad
Dv_2=-v_1v_2+v_3,\quad
Dv_3=v_2^2-24X-25v_2X-5v_1v_2X+10v_3X,\quad
DX=5X(X-1).
\]

## Bell/Wronskian generators and centered generators

\[
\begin{aligned}
J_1&=v_1,\\
J_2&=v_1^2+Dv_1=-2v_2+5Xv_1-10X,\\
J_3&=(D+v_1)J_2=-2v_3-10Xv_2+(50X^2-35X)v_1+50X-100X^2.
\end{aligned}
\]

These raw Bell generators contain pure functions of \(X\). We therefore center them:
\[
\begin{aligned}
H_1&=J_1=v_1,\\
H_2&=J_2+10X=-2v_2+5Xv_1,\\
H_3&=J_3-50X+100X^2=-2v_3-10Xv_2+(50X^2-35X)v_1.
\end{aligned}
\]
The inverse triangular transformation is
\[
v_1=H_1,\qquad v_2=\tfrac12(5XH_1-H_2),\qquad
v_3=\tfrac12(-H_3+5XH_2+(25X^2-35X)H_1).
\]
At fixed \(X\), all \(H_i=0\) iff all \(v_i=0\). Consequently, the pure-\(X\) projection is preserved by centering.

## Genus two: full centered expression

\[
\begin{aligned}
P_2={}&\underbrace{\frac{3125}{144}-\frac{28795}{144}X
+\frac{205}{144}X^2+\frac{25}{48}X^3}_{f_2^{(0)}(X)}
-\frac5{24}H_1^3
+\left(\frac{125}{24}+\frac{45}{4}X\right)H_1^2\\
&-\left(\frac{15625}{288}+\frac{130085}{144}X
+\frac{2825}{288}X^2\right)H_1
+\left(\frac{3125}{72}+\frac{14125}{72}X
-\frac{25}{12}H_1\right)H_2-\frac{175}{9}H_3.
\end{aligned}
\]

The number of monomials is **14** in the original \(v\), raw \(J\), and centered \(H\) bases.

## Genus three: ambiguity and higher generator dependence

The full \(P_3\) is 62 terms in each of the three polynomial bases (all coefficients are encoded in the companion Python file). Its pure-\(X\) projection, unchanged in the centered \(H\) basis, is
\[
\boxed{f_3^{(0)}(X)=\frac{5}{72576}\bigl(
781250+27683000X-112828006X^2+3339968X^3+1696500X^4+477000X^5+225000X^6
\bigr).}
\]
The \(H_3\)-dependence is at most quadratic. In particular,
\[
[H_3^2]P_3=-\frac{8225}{108},
\]
while the complete linear coefficient is
\[
\begin{aligned}
[H_3]P_3=-\frac{235}{1728}\bigl(&12H_1^3-648H_1^2X-300H_1^2+120H_1H_2+565H_1X^2\\
&+52034H_1X+3125H_1-11300H_2X-2500H_2\\
&-30X^3-82X^2+11518X-1250\bigr).
\end{aligned}
\]
The other 41-plus terms in \(P_3\) are generated and checked exactly by the code, as are all coefficients of \(X^kH_1^aH_2^bH_3^c\).

## What changes in the holomorphic ambiguity?

A *holomorphic ambiguity* is the freedom to add \(f_g(X)\) in the chosen YY gauge. The full \(P_g\) has no unique split into \(\mathrm{BCOV\ particular}+f_g\) until one specifies how the particular solution is normalized.

One possible convention is to set the nonholomorphic part to zero on the section \(v_1=v_2=v_3=0\); then the polynomial \(f_g^{(0)}(X)=P_g(0,0,0,X)\) is that convention's holomorphic term. This section is not a physical holomorphic limit and is **not** the conifold gap condition. A change to the centered \(H_i\) preserves it exactly. A change to the raw \(J_i\) does not, and shifts the projection:

\[
\begin{aligned}
f_2^{(J)}(X)-f_2^{(0)}(X)&=\frac{625}{36}X(X+81),\\
f_3^{(J)}(X)-f_3^{(0)}(X)&=\frac{125X}{2592}\left(3375X^4+25521X^3+40930X^2-5737125X-223125\right).
\end{aligned}
\]

**This shift is a change of the chosen particular/holomorphic splitting, not a change in the physical free energy or its intrinsic ambiguity freedom.** The vector space of allowed ambiguity polynomials remains \(\mathbb Q[X]_{\le 3g-3}\).

## Conifold power basis

Setting \(P=1-X\), the projected genus-two ambiguity becomes
\[
f_2^{(0)}(1-P)=-\frac{12695}{72}+\frac{1760}{9}P+\frac{215}{72}P^2-\frac{25}{48}P^3.
\]
For genus three,
\[
\begin{aligned}
f_3^{(0)}(1-P)={}&-\frac{7020115}{1296}+\frac{10561435}{864}P
-\frac{211210255}{36288}P^2-\frac{757655}{567}P^3\\
&+\frac{1035625}{2016}P^4-\frac{18125}{144}P^5+\frac{15625}{1008}P^6.
\end{aligned}
\]
The gap condition constrains the *complete* holomorphic-limit free energy (the particular solution plus the ambiguity, with frame factors). It is not a requirement that the coefficients of the displayed projected polynomial separately vanish.

## Fourier-dual basis for boundary and modular data

Let \(q=e^{-t}=z\exp\bigl(5\sum_{n\ge1}a_n(H_{5n}-H_n)z^n/\omega_0(z)\bigr)\),
\(a_n=(5n)!/(n!)^5\), \(\omega_0(z)=\sum_{n\ge0}a_nz^n\), and let
\[
\xi=-X/3125=z/(1-3125z),\quad
E_g(\xi)=\omega_0(z(\xi))^{2g-2}(1-3125z(\xi))^{g-1}.
\]
In YY conventions, \(F_g^A =(-1/625)^{g-1}E_g P_g^{\mathrm{hol}}\). For the ambiguity subspace, define the exact truncated Fourier-dual polynomials
\[
\Psi_{g,k}(\xi)=\left[\frac{q(\xi)^k}{E_g(\xi)}\right]_{\le 3g-3}.
\]
Thus \(E_g(\xi)\Psi_{g,k}(\xi)=q^k+O(q^{3g-2})\), and the coordinates in the \(\Psi\) basis couple *diagonally* to the Fourier coefficients of the ambiguity contribution through the target degree.

For genus two, an explicit calculation yields
\[
\begin{aligned}
\Psi_{2,0}&=1+2885\xi-183600\xi^2+312150000\xi^3,\\
\Psi_{2,1}&=\xi+530\xi^2-1010375\xi^3,\\
\Psi_{2,2}&=\xi^2-1825\xi^3,\\
\Psi_{2,3}&=\xi^3.
\end{aligned}
\]
For genus three, e.g.,
\[
\begin{aligned}
\Psi_{3,0}&=1+5770\xi+7956025\xi^2-435072000\xi^3+\cdots,\\
\Psi_{3,1}&=\xi+3415\xi^2+335075\xi^3+\cdots,\\
\Psi_{3,4}&=\xi^4-3650\xi^5+10748375\xi^6,\\
\Psi_{3,5}&=\xi^5-6005\xi^6,\quad\Psi_{3,6}=\xi^6.
\end{aligned}
\]
Notice that, *in these calculations*, all basis-change coefficients are integers. This reflects the integral power series of the quintic mirror map and the fact that \(E_g(0)=1\). The rational numbers reside in the physical polynomial's coordinates, while the reorganization itself is integral.

The genus-two projected ambiguity has the Fourier-dual coordinates
\[
(b_{2,0},b_{2,1},b_{2,2},b_{2,3})=
\left(\frac{3125}{144},\frac{40484375}{72},-\frac{40337734375}{144},\frac{1232037109375}{36}\right).
\]
The seven exact genus-three coordinates and every \(\Psi_{3,k}\) are printed by the program.

## Independent enumerative sanity checks

After applying the original YY holomorphic limit, the corresponding period factors, and the exact mirror map, the code computes
\[
F_2^A=-\frac5{144}+\frac{575}{48}q+\frac{5125}{2}q^2+\frac{7930375}{6}q^3+\cdots,
\]
\[
F_3^A=\frac5{36288}+\frac{2875}{6048}q+\frac{35125}{336}q^2+\frac{9915125}{189}q^3+\cdots.
\]
These match the known genus-zero \(n_1^0=2875\), \(n_2^0=609250\) low-degree GV multiple-cover contributions in the sectors where higher-genus invariants vanish.

## Interpretation, limitations, and next test

The centered basis achieves an *exact ambiguity-preserving factorization*, but neither \(H_i\) nor \(X^kH_i\) has yet been proved to transform as a scalar modular form. The Fourier-dual basis is compatible with bootstrap information, but modular predictions are not imposed by the change of basis itself. A full bootstrap test must translate D4--D2--D0 modular coefficients into GV/PT constraints and apply those linear or nonlinear conditions to the Fourier-dual ambiguity coordinates.

The calculations suggest testing \(X^kH_i\) and Wronskian-covariant combinations **alongside** (not instead of) the Fourier-dual basis. The files provide a reproducible small-genus benchmark for such experiments.
