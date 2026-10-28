# Mirror quintic: the midpoint involution and the residual ambiguity

**24 October 2026 — exact rational calculation in the canonical Yamaguchi–Yau gauge**

## Executive result

The coordinate involution exchanging the MUM and conifold coordinates is
`x -> 1-x` for `x=psi^(-5)=3125*z_PF`. It sends

`X=x/(x-1) -> 1/X`; the fixed point is `x=1/2`, `X=-1`.

This is **not an automorphism of the mirror-quintic variation of Hodge
structure**. The MUM and conifold singularities have different local
monodromy. Therefore a crossing relation between their local expressions
cannot be imposed as an extra boundary condition.

Nevertheless, the exact boundary kernel admits a concise all-genus
algebraic description, and its midpoint jets form an invertible set of
coordinates *if midpoint observables can be independently computed*.

## 1. Transformation of the Picard–Fuchs system

The standard operator in the normalized conifold coordinate `x` is

    L_x = theta_x^4 - x prod_{k=1}^4(theta_x + k/5),
    theta_x=x*d/dx.

With `w=1-x`,

    theta_x=-(1-w)*d/dw,
    L_x -> [-(1-w)*d/dw]^4
           -(1-w)*prod_{k=1}^4([-(1-w)*d/dw]+k/5).

It is not the standard quintic operator with `theta_w=w*d/dw`.
The local exponent multisets at the two singularities are `(0,0,0,0)`
(MUM, at x=0) and `(0,1,1,2)` (conifold, at x=1). In particular the
monodromy logarithms have different Jordan type. Therefore `x -> 1-x`
is *only a coordinate interchange* and not a physical self-duality.

`X(1-x)=1/X(x)` is an exact rational identity; it does not imply
`a_g(X)=X^(3g-3)*a_g(1/X)` or any palindromicity.

## 2. The section coefficient and gluing at the midpoint

The canonical polynomial ambiguity is

    a_g(X)=sum_{j=0}^{3g-3} a_{g,j} X^j.

The natural *reduced section coefficient* (after suppressing the common
period normalization `omega(x)^(2g-2)/625^(g-1)` in a fixed holomorphic
frame) is

    f_g(x)=(x-1)^(g-1) a_g(x/(x-1)).

In the conifold coordinate `w=1-x`, the same coefficient is

    f_g(1-w)=(-w)^(g-1) a_g((w-1)/w).

At `w=x=1/2` these expressions coincide **identically for every set
of ambiguity coefficients**. Derivatives satisfy

    d_w^n f_g(1-w) |_{w=1/2} = (-1)^n d_x^n f_g(x)|_{x=1/2}.

A further transition function for the Hodge line-bundle frame multiplies
the local coefficient, but still introduces no constraints on a global
section. Changing the holomorphic limit/polarization requires the
appropriate propagator/special-geometry transformation, not merely
substituting periods.

**Counterexample to a fictitious x <-> 1-x symmetry:**

    a_2(X)=3125/144-(28795/144)X+(205/144)X^2+(25/48)X^3.
    a_2(-1)=16025/72 != 0.
    f_2(1/3)=-140515/1728,
    f_2(2/3)=-60935/432.

Neither simple reciprocal palindromicity nor `f_2(x)=f_2(1-x)` holds.

## 3. Closed form of the *homogeneous* boundary kernel

Fix the lower-genus BCOV anomaly and impose on two ambiguity choices the
same conifold gap, orbifold regularity, and MUM constant-map term. Let
`delta a_g` be their difference. Define

    m=ceil((3g-3)/5),   D=floor((2g-2)/5)=g-1-m.

Then:

1. At the conifold, the holomorphic frame factor behaves as
   `(1-x)^(g-1)`; removing every negative Laurent power forces
   `deg_X(delta a_g)<=g-1`.
2. At the orbifold, the frame carries `psi^(-3(g-1))`, while
   `1-X ~ -psi^5`; regularity forces a zero of order at least `m`
   at `X=1`.
3. MUM constant-map equality forces `delta a_g(0)=0`.

It follows that

    delta a_g(X) = X*(1-X)^m*R_g(X),
    deg R_g <= D-1.

Since `1-X=-1/(x-1)`,

    delta f_g(x)=(-1)^m*x*(x-1)^(D-1)*R_g(x/(x-1))
                =x*Q_g(x),
    deg Q_g <=D-1.

These are vector spaces of dimension `D`. No new midpoint condition was
used in obtaining this result. For `D=0`, the homogeneous space is zero.

A particularly convenient exact basis is

    delta a_g^(r)(X)=X^(r+1)*(1-X)^m, r=0,...,D-1,
    delta f_g^(r)(x)=(-1)^m*x^(r+1)*(x-1)^(D-1-r).

The standard boundary equations annihilate all these vectors; this was
verified by exact rational arithmetic for genera 4,5,6,7,8,9. Their
count equals the kernel dimension in each genus. Existing independently
computed boundary ranks for genera 10,11,12 agree with the formula.

## 4. Midpoint collocation: what it can and cannot do

If *independent* physical midpoint data for a chosen global holomorphic
frame were available, the D values

    J_r = (d/dx)^r (f_g - f_g^particular)(1/2), r=0,...,D-1

would determine all D ambiguity parameters. This is because the Wronskian
of the simple polynomial basis `x,x^2,...,x^D` has determinant

    det[d^r(x^k)/dx^r |_{x=1/2}]
       = (1/2)^D prod_{j=1}^{D-1} j! != 0,

with rows `r=0,...,D-1` and columns `k=1,...,D`.
The basis change to `x^(r+1)*(x-1)^(D-1-r)` has determinant +/-1.

In genus 7, `m=4`, `D=2`, and

    delta a_7(X)=X*(1-X)^4*(c0+c1*X),
    delta f_7(x)=c0*x*(x-1)+c1*x^2.

Therefore

    J0=delta f_7(1/2)=(-c0+c1)/4,
    J1=delta f_7'(1/2)=c1,
    c1=J1, c0=J1-4*J0.

At genus 9, `m=5`, `D=3`, and all three coefficients are determined by
value, first derivative, and second derivative at the midpoint.

**Critical distinction:** evaluating the MUM and conifold coordinate
*descriptions* at the midpoint produces no values for the `J_r` and hence
no new exact equations. Those values need an independent observable
(e.g. an independently established enumerative/modular condition or a
separate non-BCOV numerical calculation). A constraint made by equating
independent approximate neural predictors is a consistency regularizer,
not a new mathematical boundary condition.

## 5. Exact computation and oracle experiment

The bundled standard-library Python program `midpoint_test.py` uses the
original corrected exact genus-2..12 reference and the unmodified BCOV
boundary engine. For each genus 4..9, it checks that the explicit
factorized basis lies in the exact conifold/orbifold/MUM-constant kernel,
and that its cardinality is the exact nullity. It also verifies
nonvanishing midpoint response for undetermined directions.

For genera 7 and 9, the script builds midpoint derivative values from
the **held-out known reference ambiguity** and then reconstructs every
ambiguity coefficient using exact rational collocation. Both experiments
pass exactly. This is a proof of reconstructability *with an oracle*, not
a discovery of new physical data or a blind coefficient prediction.

Run:

    python midpoint_test.py --genera 4 5 6 7 8 9

Results and exact numerical rational vectors are in `results.json`.

## 6. Implications for the proposed AI program

There is a clear basis change worth adopting immediately. Instead of
predicting the intermediate YY coefficients directly, predict `D`
coefficients of `Q_g` in `delta f_g=x Q_g(x)` (or equivalently `R_g(X)`),
project them through the explicit exact triangular basis, and reconstruct
the full ambiguity. The MUM/GV response on this space is low-degree and
triangular after expanding the known mirror map. Use exact arithmetic
for this map and integer-lattice constraints on the unknown BPS indices.

Two independently trained local-amplitude heads, expanded about MUM
and conifold and analytically continued to the midpoint with correct
line-bundle/polarization transitions, can be assessed for consistency.
Such a two-head loss can detect approximation error, but if both heads
are generated from the same rational ambiguity vector it is identically
zero and cannot disambiguate candidates.

At genus 54, `D=21` before positive-degree MUM/GV constraints. Known
low-degree geometric vanishings leave the one degree-21 BPS input
unfixed. A separately computed midpoint amplitude value could in
principle determine that remaining parameter; applying `x -> 1-x`
alone cannot.

**Reference conventions:** Yamaguchi–Yau polynomial formalism;
standard quintic PF operator; and BCOV genus-g amplitudes as sections
of the appropriate Hodge line bundle.
