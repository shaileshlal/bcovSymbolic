# Exact BCOV algebra in Yamaguchi–Yau generators

A runnable preliminary quintic implementation, using Python 3.10+ and the standard
library. It computes normalized B-model anomaly polynomials and keeps boundary
ambiguities separate. No SymPy, Sage, or floating point arithmetic is needed. The physical driver
adds exact period and mirror-map series, without explicitly inverting z(q).

## Physical boundaries and the genus-53 target (version 0.2)

See [BOUNDARIES.md](BOUNDARIES.md) for the new exact boundary solver, local
monodromy frames, symbolic coordinate scale mu, Castelnuovo/GV constraints,
and checkpointed physical driver:

```bash
python run_physical.py --plan --genus 53
python run_physical.py --genus 6 --gv-degree 4 --output my_run
python run_physical.py --genus 53 --mu mu --output my_run
```

All 43 tests pass (25 algebra tests and 18 boundary/driver tests). Physical
recursion through genus eight, with GV checks through degree four, is included
in `validated_physical/`. The physical driver fixes each ambiguity before proceeding. Genus 53 is an
unbenchmarked target; resource limits retain the last completed checkpoint.
The original `bcov_yy.py` demo below remains a **formal** recursion demo.

## Original algebra demo

Unzip, enter `bcov_yy`, and run:

```bash
python worked_example.py
python -m unittest discover -s tests -v
python bcov_yy.py --genus 3 --checkpoint my_checkpoint.json
```

`bcov_yy.py` is also a standalone module and executable. The other files provide
tests, a step-by-step example, an optional benchmark, and installation metadata.
The bundled `example_checkpoint.json` contains the reference genus-two amplitude
and the genus-three anomaly with its ambiguity unresolved.

## Exact representation and atomic expressions

`Poly` stores a sparse dictionary

```text
(exponent_u, exponent_v1, exponent_v2, exponent_v3, exponent_X) -> Fraction
```

Only nonzero terms survive. All arithmetic produces canonical sparse polynomials;
float input is rejected. Weights are `(1,1,2,3,1)`. Multiplication, powers, exact
scalar division, partial derivatives, primitives, coefficient extraction, simultaneous
substitution, exact evaluation, and weighted components are implemented.

The normalized quintic ring is **QQ[u,v1,v2,v3,X]**. Modulus dependence is carried
by X, so this layer needs rational numbers as coefficients, rather than a generic
QQ(z) rational-function engine. This is a deliberate specialization of our earlier
sparse-propagator design. Polynomial division and arbitrary rational functions
are outside the present algebra.

`RecipeArena` adds an interned arithmetic DAG above the sparse layer. Lower-genus
amplitudes, derivative vertices, and repeated products are shared. Construction
does not expand the DAG; `recipe.expand()` explicitly evaluates it, with memoization.
Genus and insertion count label each covariant derivative. `repr(poly)` reports
size instead of printing the expression; `str(poly)` requests the full polynomial.

The gradient solver currently requires expanded sparse polynomials. The DAG
preserves the computation recipe and reduces repeated work; it is not a claim
that all intermediate expansion has been eliminated. The monomial derivative
cache is bounded. `solver.compact_recipes()` discards temporary DAGs and caches
while retaining amplitudes. Release any externally held old recipe handles too.

## Mathematical conventions

We use the **original quintic psi convention**, not rescaled z-coordinate YY
generators. In [Yamaguchi–Yau, hep-th/0406078v2](https://arxiv.org/abs/hep-th/0406078),
see equations (3.1), (3.7), (3.16), (3.22), (3.28)–(3.30), and (4.1).

Let `theta = psi*d/dpsi`, `z=(5*psi)^(-5)`, and define

```text
X = 1/(1-psi^5) = 3125*z/(3125*z-1)
C = psi^3 C_(psi psi psi)
P_g = C^(g-1) F_g
B_p = theta^p(exp(-K))/exp(-K),   A = theta(G)/G
u = B_1
v1 = A+1+2u
v2 = B_2-u*v1
v3 = B_3+u*v2-5uX(v1-2)
```

The factors of 5 in v3 are convention dependent. Do not mix this implementation
with differently rescaled v-generators. `theta_z()` implements `z*d/dz=-theta/5`,
not `d/dz`. The closure table is

```text
theta(u)  = v2+u*v1-u^2
theta(v1) = -v1^2-2v2-10X+5v1X
theta(v2) = -v1*v2+v3
theta(v3) = v2^2-24X-25v2X-5v1*v2X+10v3X
theta(X)  = 5X(X-1)
```

The algebra is commutative; theta acts by the Leibniz rule. There is no need to
represent arbitrary noncommutative products of differential operators.
For normalized vertices `P_g^(n)=C^(g-1) psi^n D_psi^n F_g`, use

```text
P_g^(n+1) = theta(P_g^(n))
            + [-n(v1-2u)+(2g-2)(u-5X/2)] P_g^(n).
```

The stable seeds are `P_0^(3)=1` and
`P_1^(1)=125/12-v1/2-28u/3+5X/12`. The genus-one seed is a one-point function,
not a polynomial representation of F_1 itself. Tensor rank increases after each
covariant derivative.

## The complete polynomial anomaly recursion

For g>=2, form the exact source

```text
R_g = -1/2 [P_(g-1)^(2) + sum_(h=1)^(g-1) P_h^(1) P_(g-h)^(1)].
```

The two anomaly equations become u independence and

```text
(d_v1 - u*d_v2 - u*(u+5X)*d_v3) P_g = R_g.
```

The source is at most quadratic in u. With `r_k=[u^k]R_g`, recover all three
derivatives as `(d_v1 P_g,d_v2 P_g,d_v3 P_g)=(r_0,-r_1+5X*r_2,-r_2)`.
The solver checks mixed-partial consistency, integrates the gradient, and verifies
the result. Integrating only the first derivative would leave spurious freedom
in v2 and v3. u is retained in the intermediate vertices even though P_g is
independent of u.

`Amplitude.anomaly` uses the integration convention `P_g(0,0,0,X)=0`.
The separate ambiguity vector has slot j multiplying X^j, through degree 3g-3.
`None` means unresolved and blocks continuation to the next genus. An explicit
empty tuple selects zero coefficients. This is an algebraic choice; physical
boundary conditions are supplied by `run_physical.py`. The full polynomial ambiguity space is
an enlarged ansatz; geometric regularity can reduce its allowed subspace.

```python
from bcov_yy import BCOV, YY_GENUS2_AMBIGUITY

s = BCOV()
a2 = s.integrate(2)                     # unresolved ambiguity
s.set_ambiguity(2, YY_GENUS2_AMBIGUITY)  # supplied published boundary data
a3 = s.integrate(3)                     # exact anomaly, seven free X slots
assert s.residual(3, a3.anomaly) == 0
s.save("checkpoint.json")
s = BCOV.load("checkpoint.json")        # verifies conventions and anomaly
# For a formal continuation only:
s.set_ambiguity(3, ())
a4 = s.integrate(4)
```

At genus two the recursion gives `d_v3 P_2=350/9` and
`d_v2 P_2=-3125/36+25v1/6-2375X/12`. Supplying the bundled pure-X vector
reproduces the complete published equation (4.1), coefficient by coefficient.
The next anomaly has 55 terms of weight at most six; the independent equation
(4.7) gives `d_v3 P_3=-47 P_2/3`, which is checked exactly.

## Validation and higher-genus extension points

The 25 focused tests cover exact arithmetic, ring identities, original AB/Picard–Fuchs
closure, covariant tensor rank, DAG sharing, published genus two, genus-three and
genus-four consequences of equation (4.7), ambiguity propagation, and checkpoint
validation/resume. A formal run through genus eight has 2,548 anomaly terms and
zero BCOV residuals at every step. It uses explicit zero ambiguities after genus
two; it is not a physical genus-eight quintic free energy. Measure locally with
`python benchmark.py --genus 8`; timings are machine dependent.

Practical next extensions:

- Substitute GMP/FLINT rational and sparse-polynomial kernels behind `Poly`;
  profile multiplication and rational reduction before changing algorithms.
- Group monomials by weight or v3 degree, reuse symmetric products, and add
  modular computation with certified rational reconstruction for large coefficients.
- Preserve more derivative computations at DAG level and stream gradient blocks
  instead of materializing the whole source. Check memory as well as run time.
- Physical boundary solvers now live in `physical_boundaries.py`. Generalizing
  the geometry or leaving genuinely unfixed lower-genus symbols still requires
  a larger coefficient algebra.
- Exact period-series/mirror-map evaluation and frame conversion now live in
  `period_frames.py`. Global integral connection matrices remain an extension.
- For another geometry or gauge, replace the closure table, normalization,
  genus-one seed, anomaly vector fields, and ambiguity basis together. Merely
  changing the Euler characteristic does not generalize this quintic engine.

The versioned JSON checkpoint uses integer strings for exact coefficients and
atomic replacement; it stores amplitude data, not executable pickle or DAG caches.
No genus-64 performance or boundary-fixing capability is claimed for this prototype.
