# Exact bootstrap-adapted YY polynomial combinations: quintic, genera two and three

**Status.** This is a reproducible **pullback** of published modular-bootstrap/near-Castelnuovo GV and PT data to the Yamaguchi–Yau polynomial ring. It verifies the published numbers and constructs exact Fourier/GV-dual polynomial combinations. It **does not** reconstruct the weight −3/2 vector-valued modular form or implement the full D4–D2–D0/PT wall-crossing machinery. In particular, the resulting polynomials are not claimed to be scalar modular forms.

References:

- Yamaguchi–Yau, *Topological String Partition Functions as Polynomials* (hep-th/0406078), equations (4.1), (4.8) and (A.1): https://arxiv.org/abs/hep-th/0406078 .
- Alexandrov–Feyzbakhsh–Klemm–Pioline–Schimannek, *Quantum geometry, stability and modularity* (arXiv:2301.08066v2), equations (3.3)–(3.12), (4.12), (5.7)–(5.8), Tables 2 and 3: https://arxiv.org/html/2301.08066 .

## 1. Variables and frame

Write

\[
 z=(5\psi)^{-5},\qquad X=\frac1{1-\psi^5}=-\frac{3125z}{1-3125z},\qquad
 \xi=-\frac X{3125}=\frac{z}{1-3125z}.
\]

Let

\[
\omega_0(z)=\sum_{n\ge0}\frac{(5n)!}{(n!)^5}z^n,\qquad
q(z)=z\exp\left[\frac{5}{\omega_0(z)}\sum_{n\ge1}\frac{(5n)!}{(n!)^5}(H_{5n}-H_n)z^n\right].
\]

The exact period-factor normalization used in the provided original YY calculation is

\[
F_g^{A}(q)=N_g E_g(\xi) P_g^{\rm hol}(\xi),\qquad
N_g=\left(-\frac1{625}\right)^{g-1},\qquad
E_g=\omega_0(z(\xi))^{2g-2}(1-3125z(\xi))^{g-1}.
\]

Here `hol` denotes the holomorphic-limit substitution in the YY polynomial, **not** evaluation at `v1=v2=v3=0`. The latter is a convenient but noncanonical projection that defines a convention for splitting a particular BCOV solution and its pure-\(X\) remainder.

## 2. GV-dual basis — exact formula

A shift in genus-\(g\) holomorphic ambiguity, with all lower genera held fixed, induces

\[
\delta F_g=N_g E_g(\xi)\,\delta P_g(X).
\]

The GV multicovers at genus \(g\) carry the factor \(r^{2g-3}\). Consequently, with \(m=3g-3\), define

\[
\boxed{\mathcal B_{g,d}(\xi)
=\left[E_g(\xi)^{-1}\sum_{1\le r\le m/d}r^{2g-3}q(\xi)^{rd}\right]_{\xi^{\le m}},\qquad 1\le d\le m,}
\]

and set \(\mathcal B_{g,0}=[E_g^{-1}]_{\le m}\) for the constant-map coordinate. The truncation is **polynomial** and is carried out with exact rational arithmetic. Since \(q(\xi)=\xi+O(\xi^2)\), the basis is triangular, with unit diagonal.

For any \(1\le d,k\le m\),

\[
\sum_{r\mid d}\mu_{\rm Mob}(r)r^{2g-3}
[q^{d/r}]\left(E_g\mathcal B_{g,k}\right)=\delta_{dk}.
\]

Thus, writing

\[
\delta P_g=\sum_{d=0}^{m}a_{g,d}\mathcal B_{g,d}
\]

gives the **homogeneous** exact response

\[
\boxed{\delta n_d^g=N_g a_{g,d},\qquad1\le d\le m,}
\]

where lower-genus data are held fixed. This is *not* a claim about degrees \(d>m\): truncation produces nonzero higher-order tails.

### Genus 2 (\(m=3\))

\[
\begin{aligned}
\mathcal B_{2,0}&=1+2885\xi-183600\xi^2+312150000\xi^3,\\
\mathcal B_{2,1}&=\xi+532\xi^2-1014022\xi^3,\\
\mathcal B_{2,2}&=\xi^2-1825\xi^3,\\
\mathcal B_{2,3}&=\xi^3.
\end{aligned}
\]

### Genus 3 (\(m=6\))

\[
\begin{aligned}
\mathcal B_{3,1}&=\xi+3423\xi^2+343582\xi^3-517462326\xi^4+940317805275\xi^5-1906438887994735\xi^6,\\
\mathcal B_{3,2}&=\xi^2+1060\xi^3-1739842\xi^4+3322135800\xi^5-6943558209848\xi^6,\\
\mathcal B_{3,3}&=\xi^3-1295\xi^4+1731250\xi^5-2018939992\xi^6,\\
\mathcal B_{3,4}&=\xi^4-3650\xi^5+10748375\xi^6,\\
\mathcal B_{3,5}&=\xi^5-6005\xi^6,\\
\mathcal B_{3,6}&=\xi^6.
\end{aligned}
\]

The expression for \(\mathcal B_{3,0}\), the full coefficient matrices, and all coefficients for the original \(P_2,P_3\) projections are in `yy_modular_bootstrap_results.txt`.

## 3. A modular PT constraint genuinely mixing genera two and three

From Tables 2–3 of the modular-bootstrap paper,

\[
\operatorname{PT}(4,-2)=8625,\qquad
\operatorname{PT}(4,-1)=569250.
\]

The near-Castelnuovo PT/GV formulas imply

\[
\operatorname{PT}(4,-2)=n_4^3,\qquad
\operatorname{PT}(4,-1)=n_4^2+4n_4^3,
\]

since \(g_{\rm max}(4)=3\) and \(2g_{\rm max}(4)-2=4\). The exact original YY data give \(n_4^3=8625\), \(n_4^2=534750\), hence \(534750+4\cdot8625=569250\).

Choose the noncanonical but *fixed* BCOV splitting \(P_g=P_g^{\rm part}+f_g^{(0)}\), where \(f_g^{(0)}=P_g(0,0,0,X)\). From the exact holomorphic calculations the responses of the non-pure-\(X\) parts at \(d=4\) are

\[
(n_4^2)_{\rm part}=-\frac{9095513274625}{144},\qquad
(n_4^3)_{\rm part}=\frac{861590899312225}{18144}.
\]

Write the genus-two ambiguity \(f_2^{(0)}=\sum_{d=0}^3a_{2,d}\mathcal B_{2,d}\). Because the genus-two basis is truncated at \(\xi^3\), its response at degree four is non-diagonal:

\[
\boxed{\delta n_4^2=-\frac1{625}
\bigl(627948300000 a_{2,0}-2203986464 a_{2,1}
-3708877 a_{2,2}+4180 a_{2,3}\bigr).}
\]

At genus three, however, \(d=4\le m=6\), so with \(f_3^{(0)}=\sum_{d=0}^6 a_{3,d}\mathcal B_{3,d}\),

\[
\boxed{\delta n_4^3=\frac{a_{3,4}}{390625}.}
\]

Therefore the pulled-back modular/PT constraints are *explicitly*

\[
\boxed{\frac{861590899312225}{18144}+\frac{a_{3,4}}{390625}=8625,}
\]

\[
\boxed{-\frac{9095513274625}{144}-\frac{627948300000a_{2,0}-2203986464a_{2,1}-3708877a_{2,2}+4180a_{2,3}}{625}
+4\left(\frac{861590899312225}{18144}+\frac{a_{3,4}}{390625}\right)=569250.}
\]

Using the YY polynomial solutions,

\[
 a_{2,0}=\frac{3125}{144},\quad a_{2,1}=\frac{40484375}{72},\quad
 a_{2,2}=-\frac{13499890625}{48},\quad a_{2,3}=\frac{2463952765625}{72},
\]

and

\[
a_{3,4}=-\frac{336558883914150390625}{18144},
\]

both relations are verified exactly. The enormous cancellation between the particular solution and the pure-\(X\) remainder explains why looking at either component separately does not reveal modular-looking coefficients.

## 4. Mixed \(X^kJ_i\) combinations and a nontrivial filtered generator

Define the centered Bell/YY generators

\[
H_1=v_1,\qquad H_2=-2v_2+5Xv_1,\qquad
H_3=-2v_3-10Xv_2+(50X^2-35X)v_1.
\]

These are related to the period-Wronskian Bell polynomials \(J_i\) by \(H_2=J_2+10X\), \(H_3=J_3-50X+100X^2\).

Their exact holomorphic jets, expressed in \(\xi=-X/3125\), begin

\[
\begin{aligned}
H_2^{\rm hol}&=-6000\xi+81756250\xi^2-40027881250\xi^3+\dots,\\
H_3^{\rm hol}&=30000\xi-363043750\xi^2-2310136093750\xi^3+\dots.
\end{aligned}
\]

Consequently \(H_3+5H_2\) has no \(\xi\) term, and its \(\xi^2\) coefficient is \(45737500\). One can remove both the \(\xi^2\) and \(\xi^3\) terms using only powers of \(X\). Define

\[
\boxed{\Theta=H_3+5H_2
-\frac{14636}{3125}X^2
-\frac{160657632}{1953125}X^3.}
\]

The resulting **mixed YY combination** obeys

\[
\boxed{\Theta^{\rm hol}=O(q^4).}
\]

In particular it is invisible to *all* genus-two or genus-three GV bootstrap tests involving only degrees 1–3, with lower genera fixed. The first nonzero responses in the two amplitude normalizations are

\[
[q^4]\bigl(N_2E_2\Theta^{\rm hol}\bigr)=-1908599241000,\qquad
[q^4]\bigl(N_3E_3\Theta^{\rm hol}\bigr)=\frac{15268793928}{5}.
\]

Likewise \(X^k\Theta\) starts at order \(q^{4+k}\), suggesting an explicit filtration on the allowed YY monomials. This construction takes the `X^power J_i` proposal seriously, but **it is a Fourier/GV filtration, not a proof that \(\Theta\) transforms modularly**. Also, although such combinations can be used as a *basis of the full fixed BCOV polynomial*, a mixed \(H\)-dependent polynomial cannot be freely added as a holomorphic ambiguity without potentially violating the anomaly equation.

## 5. Exact GV checks from the original \(P_2,P_3\)

The full A-model series, after all period factors, begin

\[
\begin{aligned}
F_2^A={}&-\frac5{144}+\frac{575}{48}q+\frac{5125}{2}q^2+\frac{7930375}{6}q^3+1010821250q^4+\frac{3092760574975}{3}q^5+\frac{45745969626745975}{24}q^6+\ldots,\\
F_3^A={}&\frac5{36288}+\frac{2875}{6048}q+\frac{35125}{336}q^2+\frac{9915125}{189}q^3+\frac{5046984125}{126}q^4+\frac{23896218123175}{756}q^5-\frac{86000485655933125}{3024}q^6+\ldots.
\end{aligned}
\]

At genus 2, subtract \(\sum_{r\mid d}r n_{d/r}^0/240\), then invert the Dirichlet convolution with \(r\). At genus 3, subtract \(\sum_{r\mid d}r^3(n^0_{d/r}/6048-n^2_{d/r}/12)\), then invert the Dirichlet convolution with \(r^3\). The exact results are

| d | \(n_d^2\) | \(n_d^3\) |
|---:|---:|---:|
| 1 | 0 | 0 |
| 2 | 0 | 0 |
| 3 | 0 | 0 |
| 4 | 534750 | 8625 |
| 5 | 75478987900 | -15663750 |
| 6 | 871708139638250 | 3156446162875 |

These coincide exactly with the cited modular-bootstrap paper, Table 2. Because these low-genus values were already accessible by direct integration, this is a consistency check rather than an independent modular determination of new genus-two/three invariants.

## 6. Computational interpretation and extension

The finite dual basis is defined entirely in \(\mathbb Q[\xi]\), computed using rational Taylor operations, the normalized period, and the mirror map. It is neither genus-universal nor a modular representation, and its uniqueness applies only modulo \(\xi^{m+1}\). The *space* of allowable holomorphic ambiguities remains \(\mathbb Q[X]_{\le3g-3}\). Its basis has changed but its dimension has not.

The next genuinely new mathematical step is to start with the vector-valued weight \(-3/2\) modular form \(h_\mu(\tau)\), compute new stable-pair invariants through the precise wall-crossing map (4.12), convert these by plethystic logarithm to \(n_d^g\), and apply the GV-response matrix shown here to ambiguities at large genus. The available modular form controls *near-Castelnuovo* data at high genus; it should not be mistaken for an ordinary modular transformation law of \(P_2\) or \(P_3\).

## 7. Reproduction

Keep both files in one directory, install SymPy, then run:

```bash
python yy_p2_p3_exact.py
python yy_modular_bootstrap.py
```

The new script tests the GV-duality identities, verifies Table 2 at \(d\le6\), checks the degree-four PT relations, computes mixed-generator pullback matrices, and verifies that the constructed \(\Theta\) begins at degree four. All calculations use exact integers and fractions.
