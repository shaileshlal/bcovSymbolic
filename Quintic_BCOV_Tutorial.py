import marimo

__generated_with = "0.25.1"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Exact BCOV on the quintic: a guided laboratory

    **A reproducible tutorial for the `bcov-symbolic` package**
    Canonical recursion: **original Yamaguchi–Yau (YY, 2004)**. Other frames: **Huang–Klemm–Quackenbush (HKQ)** and **Gu/Klemm small-$z$**.

    ### What we will learn

    1. Build and differentiate sparse polynomials in $\mathbb Q[u,v_1,v_2,v_3,X]$.
    2. Track the *coordinate*, the *Yukawa normalization*, and the *propagator gauge* as separate choices.
    3. Solve the BCOV anomaly through genus 2 and then genus 3, supplying—not inventing—the holomorphic ambiguity.
    4. Extract $[X^k]P_g$ and distinguish it from the **pure-$X$ coefficient** $[X^k]P_g|_{v_1=v_2=v_3=0}$.
    5. Extract the latter efficiently from the factored Klemm free energies, without constructing enormous polynomials.
    6. Test Bernoulli recurrences, denominator invariance, and the physical GW/GV $q$-expansion.

    **Status of this tutorial.** The repository's existing physical data are explicit through genus 15, and the built-in $q$-series checks are to degree 6. This notebook does **not** claim to fix the ambiguity at arbitrary genus. Every displayed rational result is exact.

    **How to run:** From the repository root (the directory containing this notebook and `bcov/`), run `marimo edit Quintic_BCOV_Tutorial.py`, or `python Quintic_BCOV_Tutorial.py` to execute it as a script. Python 3.10+; all core calculations use only the standard library. marimo is needed to display the notebook.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 0. Locate the package and read its convention ledger

    Install the package with `python -m pip install -e ".[notebooks]"`. Imports use the installed package; generated files go into `notebook_output/` in the current directory.
    """)
    return


@app.cell
def _():
    from pathlib import Path
    import sys, json, gzip, csv
    from fractions import Fraction as Q
    from math import comb, lcm

    import bcov
    PACKAGE = Path(bcov.__file__).resolve().parent
    DATA = PACKAGE / 'data'
    conventions = json.loads((PACKAGE / 'CONVENTIONS.json').read_text())
    print('Package:', PACKAGE)
    print('Convention version:', conventions['schema_version'])
    print('Canonical frame:', conventions['canonical_ring_frame'])
    print('Available data:', sorted(p.name for p in DATA.glob('ff_*.txt'))[:4], '...')
    return DATA, PACKAGE, Path, Q, comb, conventions, csv, gzip, json, sys


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. The three coordinate frames

    The source of the earlier genus-two mismatch was comparing objects with **the same variable names but different underlying coordinates**:

    | Object | Original YY (`yy`) | HKQ (`hkq`) | Gu/Klemm (`small_z`) |
    |:--|:--|:--|:--|
    | Coordinate | $s$ | $\psi=s^5$ | $z=(3125\psi)^{-1}$ |
    | Large radius | $s=\infty$ | $\psi=\infty$ | $z=0$ |
    | Euler derivation | $\theta_s=s\partial_s$ | $\theta_\psi=\theta_s/5$ | $\theta_z=-\theta_\psi$ |
    | Discriminant generator | $X=(1-s^5)^{-1}$ | $X=(1-\psi)^{-1}$ | $X_z=(1-3125z)^{-1}$ |

    Hence $X_{\rm HKQ}=1-X_z$. The connection generators obey

    $$
    (u,v_1,v_2,v_3)_{\rm YY}=(5u,\;5v_1,\;25v_2,\;125v_3)_{\rm HKQ},
    $$

    $$
    (u,v_1,v_2,v_3)_{\rm HKQ}=(-u,\;-v_1,\;v_2,\;-v_3)_{z}.
    $$

    We **always transform the whole normalized amplitude** $P_g=C^{g-1}F_g$, not only its variable names:

    $$
    P_g^{\rm YY}(v_i,X)=5^{3g-3}P_g^{\rm HKQ}(v_i/5^i,X),\qquad
    P_g^{\rm HKQ}(v_i,X)=(-1)^{g-1}P_g^z(-v_1,v_2,-v_3,1-X).
    $$

    The $i$ in $5^i$ means $i=1,2,3$ for $v_1,v_2,v_3$; $u$ scales by 5 too. The prefactor $(-1)^{g-1}$ comes from the cubic Yukawa transformation. The factor $(2\pi i)^3$ occurring in some Yukawa formulas must also be tracked at the period/three-form normalization level; it is not a numerical parameter to tune after the fact.
    """)
    return


@app.cell
def _():
    from bcov.bcov_yy import (Poly, YYRing, BCOV, u, v1, v2, v3, X,
                         YY_GENUS2_AMBIGUITY, reference_genus2)
    from bcov.quintic_frames import (FramedPolynomial, yy_to_hkq, hkq_to_yy,
                                hkq_to_z, z_to_hkq, yy_to_z,
                                theta_hkq, theta_z, GU_GAUGE)

    ring = YYRing()
    assert ring.theta(X) == 5*X*(X-1)         # original YY fifth-root coordinate
    assert theta_hkq(X) == X*(X-1)            # HKQ coordinate
    assert theta_z(X) == X*(X-1)              # this X denotes X_z here!
    print('YY θ(X)  =', ring.theta(X))
    print('HKQ θ(X) =', theta_hkq(X))
    print('z θ(Xz)  =', theta_z(X))
    return (
        BCOV,
        FramedPolynomial,
        GU_GAUGE,
        Poly,
        X,
        YY_GENUS2_AMBIGUITY,
        hkq_to_yy,
        reference_genus2,
        ring,
        v1,
        v2,
        v3,
        yy_to_hkq,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### A rule that prevents silent convention errors

    Use `FramedPolynomial` to **tag** which coordinate frame an expression belongs to. Operations on objects with different frame labels fail rather than returning a plausible-looking wrong polynomial. Unlike a symbolic relabeling, `.to(...)` converts *both* generators and the genus-dependent normalization.
    """)
    return


@app.cell
def _(FramedPolynomial, reference_genus2, yy_to_hkq):
    P2_YY = reference_genus2()
    yy = FramedPolynomial(P2_YY, 2, 'yy')
    hkq = yy.to('hkq')
    zframe = hkq.to('small_z')
    assert zframe.to('yy').poly == P2_YY
    assert yy_to_hkq(P2_YY, 2) == hkq.poly

    try:
        _ = yy + hkq  # deliberately wrong: two frames
    except ValueError as error:
        print('Caught the intended mixed-frame error:', error)
    print('Round trip YY → HKQ → small z → YY: exact.')
    return P2_YY, yy


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Elementary exact polynomial algebra

    A monomial is encoded by five exponents ordered as $(u,v_1,v_2,v_3,X)$. Polynomial arithmetic uses `Fraction`; **floats are deliberately forbidden**. The usual filtration weights are $(1,1,2,3,1)$.

    `Poly.coefficient("X",k)` returns a **polynomial in the remaining variables**, not a number. If we want the coefficient of the *pure* $X^k$ term, all three $v_i$ must be set to zero as a second operation.
    """)
    return


@app.cell
def _(Poly, Q, X, v1, v2):
    example = Q(7,12)*v1*X**2 - Q(5,3)*v2*X + Q(11,8)*X**2
    print('Example:', example)
    print('Weighted degree:', example.degree())
    print('[X^2] example:', example.coefficient('X', 2))
    print('∂/∂v1 example:', example.partial('v1'))
    assert example.coefficient('X', 2) == Q(7,12)*v1+Q(11,8)

    try:
        Poly.constant(0.1)
    except TypeError:
        print('Float rejected; use Q(1,10).')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. BCOV's differential algebra and the genus-one seed

    In the original YY variables, $\theta=s\partial_s$ is a **derivation of the differential ring**, while covariant differentiation also includes a connection term:

    $$
    \mathcal D_{g,n}P=\theta P +\left[-n(v_1-2u)+(2g-2)\left(u-\tfrac52X\right)\right]P.
    $$

    For example, `vertex(1,1)` is the specified **genus-one one-point seed** $P_1^{(1)}$; it is **not** the genus-one free energy $F_1$.
    """)
    return


@app.cell
def _(BCOV, X, ring, v1, v2):
    solver = BCOV()
    print('P_1^(1) =', solver.vertex(1, 1).expand())
    print('θ(v1) =', ring.theta(v1))
    print('D_(2,0)(X) =', ring.covariant(X, genus=2, insertions=0))
    assert ring.theta(v1*v2) == ring.theta(v1)*v2+v1*ring.theta(v2)
    return (solver,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Genus two: anomaly versus holomorphic ambiguity

    The recursion implements the polynomial anomaly equation

    $$
    \left(\partial_{v_1}-u\partial_{v_2}-u(u+5X)\partial_{v_3}\right)P_g=R_g,
    \quad R_g=-\frac12\left(P_{g-1}^{(2)}+\sum_{h=1}^{g-1}P_h^{(1)}P_{g-h}^{(1)}\right).
    $$

    The source $R_g$ is quadratic in $u$. The solver checks mixed partials and integrates all three anomaly gradients, choosing the particular primitive with

    $$P_g^{\mathrm{an}}(v_1=v_2=v_3=0,X)=0.$$

    **That is an integration convention, not a physical boundary condition.** The general answer is $P_g=P_g^{\mathrm{an}}+a_g(X)$, where $a_g(X)$ must be supplied or independently fixed by geometric data. `ambiguity=None` means **unresolved**; zero ambiguity is a separate, explicit choice.
    """)
    return


@app.cell
def _(P2_YY, Poly, YY_GENUS2_AMBIGUITY, reference_genus2, solver):
    g2 = solver.integrate(2)
    assert g2.ambiguity is None
    assert g2.anomaly.substitute({'v1':0,'v2':0,'v3':0}) == Poly()
    assert solver.residual(2, g2.anomaly) == Poly()
    print('Anomaly (non-holomorphic part), terms:', len(g2.anomaly))
    print('Genus-two holomorphic ambiguity: unresolved =', g2.ambiguity is None)
    try:
        g2.polynomial()
    except ValueError as error:
        print('Expected: cannot request the physical polynomial yet:', error)

    # Published YY genus-two boundary data, NOT calculated by the anomaly alone.
    g2 = solver.set_ambiguity(2, YY_GENUS2_AMBIGUITY)
    P2 = g2.polynomial()
    assert P2 == P2_YY == reference_genus2()
    assert solver.residual(2) == Poly()
    print('Physical genus-two P_2 reproduced, monomials:', len(P2))
    print('Pure-X ambiguity vector:', list(g2.ambiguity))
    return (P2,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5. Extracting $X$ coefficients correctly

    There are **three different** questions:

    1. **Full coefficient polynomial** $[X^k]P_g$: may involve $v_1,v_2,v_3$.
    2. **Pure-$X$ coefficient** $a_{g,k}=[X^k]P_g|_{v_1=v_2=v_3=0}$ in our canonical primitive convention.
    3. **Physical $q$-series coefficient** $[q^d]F_g^A$, obtained only after taking the holomorphic limit, multiplying by periods and applying the mirror map. This is *not* generally a GV integer.

    For $P_g$ independent of $u$, the first two are obtained directly from `Poly`. The second equals an entry of `Amplitude.ambiguity` **only when both use the same frame and integration convention**.
    """)
    return


@app.cell
def _(P2, Poly, Q, YY_GENUS2_AMBIGUITY):
    ZERO_V = {'u':0,'v1':0,'v2':0,'v3':0}

    def full_X_coefficient(p: Poly, k: int) -> Poly:
        '''[X^k]p, retaining every other generator.'''
        return p.coefficient('X', k)

    def pure_X_polynomial(p: Poly) -> Poly:
        '''Project to v1=v2=v3=u=0, NOT to the physical holomorphic limit.'''
        return p.substitute(ZERO_V)

    def pure_X_vector(p: Poly, max_degree: int | None = None) -> list[Q]:
        h = pure_X_polynomial(p)
        D = max_degree if max_degree is not None else max(0, h.degree('X'))
        return [h.coefficient('X', k).evaluate({}) for k in range(D+1)]

    print('[X]P2, retaining v:', full_X_coefficient(P2, 1))
    print('P2 at all v=0:', pure_X_polynomial(P2))
    coeff2 = pure_X_vector(P2, 3)
    print('Physical YY pure-X coefficients:', coeff2)
    assert coeff2 == list(YY_GENUS2_AMBIGUITY)
    assert full_X_coefficient(P2, 1) != Poly.constant(coeff2[1])
    return ZERO_V, coeff2, pure_X_vector


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 6. The Klemm/Gu propagator gauge is a separate choice

    The uploaded `ff_g.txt` files are expressed in the Gu/Klemm small-$z$ frame with the **specific holomorphic gauge** recorded in the original notes. For example,

    $$
    S^{zz}=-\frac{z^2}{5X_z}\left(v_{1,z}+\frac35\right),\qquad
    \widetilde S^z=\frac{z}{5X_z}\left(v_{2,z}+\frac{2}{25}\right),
    $$
    $$
    \widetilde S=\frac{1}{10X_z}\left[v_{3,z}+(1-X_z)\left(v_{2,z}+\frac4{125}\right)\right].
    $$

    The $q$-functions here are **holomorphic functions entering the propagator differential equations**, **not** the mirror-map $q=e^{2\pi it}$. The converter assumes this gauge. If you choose a different gauge, you must rederive the propagator identities and rebuild the converted source; **do not simply alter `GU_GAUGE` and reuse the old atomized free energies**.
    """)
    return


@app.cell
def _(GU_GAUGE, P2, conventions, yy):
    print('Gu/Klemm gauge used by the checked conversion:')
    for name, value in GU_GAUGE.items():
        print(f'  {name:9s} = {value}')

    # A quick convention diagnostic, rather than fitting unknown factors.
    assert conventions['mu'] == '1/3125'
    assert conventions['physical_Yukawa_small_z'] == 'Czzz=5*Xz/z^3'
    assert yy.to('hkq').to('small_z').to('yy').poly == P2
    print('Gauge ledger and frame round trip agree.')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 7. Read the high-genus free energies without expanding them

    Fully expanding genus 15 in $v_1,v_2,v_3,X$ would be wasteful. The package instead stores a **factored atomic representation**, one JSON record per genus:

    $$
    P_g^{\rm HKQ}=\sum_{r,n,j,k,\ell}c_{rnjk\ell}
    (1-X)^r X^n
    \underbrace{(v_1-3/5)^j}_{\alpha^j}
    \underbrace{(v_2+2/25)^k}_{\beta^k}
    \underbrace{\big(v_3-X(v_2+4/125)\big)^\ell}_{\gamma^\ell}.
    $$

    The file `data/bcov_hkq_atomic_genus_2_15.jsonl.gz` was obtained from the uploaded `ff_2.txt`–`ff_15.txt` after the coordinate and propagator conversions. The coefficients are stored as numerator and denominator **strings**, preserving exact integers.
    """)
    return


@app.cell
def _(DATA, gzip, json):
    ATOMIC = DATA / 'bcov_hkq_atomic_genus_2_15.jsonl.gz'

    def load_genus_atomic(g: int) -> dict:
        with gzip.open(ATOMIC, 'rt', encoding='utf8') as stream:
            for line in stream:
                obj = json.loads(line)
                if obj['genus'] == g:
                    return obj
        raise KeyError(f'No supplied atomic free energy for genus {g}')

    obj2 = load_genus_atomic(2)
    print('Genus two atomic basis:', obj2['basis'])
    print('Genus two atomic terms:', len(obj2['terms']))
    print('First term [r,n,j,k,ell,num,den]:', obj2['terms'][0])
    return ATOMIC, load_genus_atomic, obj2


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Small genera: reconstruct the entire polynomial as an independent check

    For genus 2 or 3, a full expansion is inexpensive. The following helper serves as a readable reference implementation. **Do not** use this approach blindly at genus 15.
    """)
    return


@app.cell
def _(P2, Poly, Q, X, hkq_to_yy, obj2, v1, v2, v3):
    alpha = v1 - Q(3,5)
    beta  = v2 + Q(2,25)
    gamma = v3 - X*(v2 + Q(4,125))

    def expand_atomic_hkq(obj: dict) -> Poly:
        '''For pedagogical small genera only; uses exact Poly arithmetic.'''
        result = Poly()
        for r,n,j,k,ell,num,den in obj['terms']:
            result += Q(int(num), int(den)) * (1-X)**r * X**n * alpha**j * beta**k * gamma**ell
        return result

    P2_from_klemm_hkq = expand_atomic_hkq(obj2)
    assert hkq_to_yy(P2_from_klemm_hkq, 2) == P2
    print('Genus two: source free energy → HKQ → YY agrees with BCOV exactly.')
    return alpha, beta, expand_atomic_hkq


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 8. Fast extraction of the *pure* $X^k$ terms from atomic data

    After setting $v_1=v_2=v_3=0$, one atomic term becomes

    $$
    c_{rnjk\ell}\left(-\frac35\right)^j
    \left(\frac2{25}\right)^k
    \left(-\frac4{125}\right)^\ell
    X^{n+\ell}(1-X)^r.
    $$

    The power $n+\ell$ is a **hard lower bound** on its possible contribution to $[X^m]$. This makes fixed low-depth extraction very fast even at high genus. We only expand $(1-X)^r$ through the requested degree, never the full multivariate polynomial.
    """)
    return


@app.cell
def _(Q, coeff2, comb, load_genus_atomic, obj2):
    def pure_X_coeffs_atomic_hkq(obj: dict, depth: int) -> list[Q]:
        '''Extract a[0],...,a[depth] at v=0, exact and without full expansion.'''
        if depth < 0: raise ValueError('depth must be nonnegative')
        result = [Q(0)] * (depth+1)
        for r,n,j,k,ell,num,den in obj['terms']:
            start = n + ell
            if start > depth:
                continue  # that term cannot reach low X powers
            prefactor = (Q(int(num),int(den))*Q(-3,5)**j
                         *Q(2,25)**k*Q(-4,125)**ell)
            for t in range(min(r, depth-start)+1):
                result[start+t] += prefactor * (-1)**t * comb(r,t)
        return result

    h2 = pure_X_coeffs_atomic_hkq(obj2, 3)
    y2 = [c*5**(3*2-3) for c in h2]
    print('HKQ genus 2: ', h2)
    print('YY genus 2:  ', y2)
    assert y2 == coeff2  # Independent from Poly's [X^k] extraction

    obj3 = load_genus_atomic(3)
    y3 = [c*5**6 for c in pure_X_coeffs_atomic_hkq(obj3, 6)]
    print('YY genus 3 pure-X vector (all 7 coefficients):')
    for k,c in enumerate(y3):
        print(f'  a_(3,{k}) = {c}')
    return obj3, pure_X_coeffs_atomic_hkq, y3


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### The *full* $[X^m]P_g(v_1,v_2,v_3)$, without expanding all of $P_g$

    The pure-$X$ projection is ideal for holomorphic ambiguity, but sometimes we need the **complete coefficient polynomial** at a particular $X$ power. We can also do that at fixed depth. The only extra step is to binomially expand the $\gamma$ factor:

    $$\gamma^\ell=\sum_{t=0}^{\ell}\binom\ell t(-1)^tX^t v_3^{\ell-t}(v_2+4/125)^t.$$

    Together with $[X^s](1-X)^r=(-1)^s\binom rs$, a term can contribute to $[X^m]P_g$ only when $n+t+s=m$. The algorithm below **never expands higher powers of $X$**; it does expand the required coefficient polynomial in the $v_i$, so it is best used for low $m$ rather than for the full genus-15 polynomial.
    """)
    return


@app.cell
def _(
    Poly,
    Q,
    ZERO_V,
    alpha,
    beta,
    comb,
    expand_atomic_hkq,
    obj3,
    pure_X_coeffs_atomic_hkq,
    v2,
    v3,
):
    def X_coefficient_polynomial_atomic_hkq(obj: dict, m: int) -> Poly:
        # Return [X^m]P_H(v1,v2,v3), retaining all v dependence.
        if m < 0: raise ValueError('m must be nonnegative')
        out = Poly()
        delta = v2 + Q(4,125)
        for r,n,j,k,ell,num,den in obj['terms']:
            if n > m:
                continue
            common = Q(int(num),int(den))*alpha**j*beta**k
            for t in range(min(ell,m-n)+1):
                s = m-n-t
                if s > r:
                    continue
                scalar = (-1)**(s+t)*comb(r,s)*comb(ell,t)
                out += common*scalar*v3**(ell-t)*delta**t
        return out

    P3_from_klemm_hkq = expand_atomic_hkq(obj3)
    for m in (0,1,2):
        direct = X_coefficient_polynomial_atomic_hkq(obj3,m)
        expanded = P3_from_klemm_hkq.coefficient('X',m)
        assert direct == expanded
        assert direct.substitute(ZERO_V) == Poly.constant(pure_X_coeffs_atomic_hkq(obj3,m)[m])
        print(f'g=3: [X^{m}]P_H has {len(direct)} v-monomials; pure coefficient = {pure_X_coeffs_atomic_hkq(obj3,m)[m]}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 9. Genus three: obtain the *physical* polynomial by supplying its ambiguity

    Genus three has degree $3g-3=6$ in $X$, hence **seven** holomorphic coefficients. The BCOV anomaly determines the nonholomorphic part but **not** those seven numbers. We now give the solver the seven values *extracted from the supplied physical genus-three free energy*. This is an independent checkpoint, not a prediction.

    The equality below checks all mixed $v_i$ terms, not merely the pure-$X$ projection.
    """)
    return


@app.cell
def _(Poly, expand_atomic_hkq, hkq_to_yy, obj3, pure_X_vector, solver, y3):
    g3 = solver.integrate(3)
    assert g3.ambiguity is None
    print('Genus-three anomaly terms:', len(g3.anomaly))
    print('Anomaly residual zero:', solver.residual(3, g3.anomaly)==Poly())

    P3_reference_hkq = expand_atomic_hkq(obj3)
    P3_reference_yy = hkq_to_yy(P3_reference_hkq,3)
    g3 = solver.set_ambiguity(3, y3)
    P3 = g3.polynomial()
    assert P3 == P3_reference_yy
    assert solver.residual(3) == Poly()
    assert pure_X_vector(P3, 6) == y3
    print('Physical genus 3: Klemm source = BCOV plus extracted ambiguity.')
    print('YY genus-three monomials:',len(P3),'; weighted degree:',P3.degree())
    return (P3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 10. A table of low-$X$ coefficients through genus 15

    For any genus supplied by Klemm, the same function works in $O(\text{number of atomic terms}\times\text{depth})$ rather than expanding all $(v_1,v_2,v_3,X)$ monomials. Beware that the resulting $a_{g,k}$ are **in our chosen zero-at-$v$ BCOV integration convention**; their numerical values change under coordinate and propagator-gauge changes.

    The first two physical YY coefficients obey the independently derived Bernoulli recursion in `audit_coefficients.py`. We check them against the **newly re-extracted** source values rather than simply plotting a fitted curve.
    """)
    return


@app.cell
def _(load_genus_atomic, pure_X_coeffs_atomic_hkq):
    from bcov.audit_coefficients import bernoulli_abs, constant as bernoulli_constant, linear as bernoulli_linear
    b_abs = bernoulli_abs(15)
    expected_a0 = {g: bernoulli_constant(g, b_abs) for g in range(2, 16)}
    summary = []
    for _g in range(2, 16):
        _aH = pure_X_coeffs_atomic_hkq(load_genus_atomic(_g), depth=4)
        _aY = [c * 5 ** (3 * _g - 3) for c in _aH]
        expected1 = bernoulli_linear(_g, b_abs, expected_a0)
        assert _aY[0] == expected_a0[_g]
        assert _aY[1] == expected1
        summary.append((_g, _aY[1] / _aY[0], _aY[2] / _aY[0]))
    print(' g       a[g,1]/a[g,0]                 a[g,2]/a[g,0]')
    for _g, r1, r2 in summary[:6] + summary[-2:]:
        print(f'{_g:2}  {str(r1):>28}  {str(r2):>28}')
    print('Passed exact a[g,0], a[g,1] Bernoulli checks at every genus 2..15.')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Export the exact coefficient table for further analysis

    Use strings such as `"1953125/36288"` in CSV rather than decimal approximations. The `frame` and `integration_convention` fields are included deliberately: omitting them invites the same normalization confusion that motivated this audit.
    """)
    return


@app.cell
def _(Path, csv, load_genus_atomic, pure_X_coeffs_atomic_hkq):
    EXPORT_DIR = Path.cwd() / 'notebook_output'
    EXPORT_DIR.mkdir(exist_ok=True)
    EXPORT_CSV = EXPORT_DIR / 'low_X_coefficients_YY_genus_2_15.csv'
    with EXPORT_CSV.open('w', newline='', encoding='utf8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['genus', 'frame', 'integration_convention', 'a0', 'a1', 'a2', 'a3', 'a4'])
        for _g in range(2, 16):
            _aH = pure_X_coeffs_atomic_hkq(load_genus_atomic(_g), 4)
            _aY = [x * 5 ** (3 * _g - 3) for x in _aH]
            writer.writerow([_g, 'yy', 'zero_at_v=0'] + [str(x) for x in _aY])
    print('Exported exact coefficients:', EXPORT_CSV.name)
    with EXPORT_CSV.open(newline='', encoding='utf8') as stream:
        assert len(list(csv.DictReader(stream))) == 14
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 11. Why individual denominators are frame dependent

    Between YY and HKQ, the pure-$X$ vector is uniformly rescaled by $5^{3g-3}$. Its normalized ratios are identical, but its powers of $5$ in denominators can differ.

    HKQ $\to$ small-$z$ sends $a_g(X)\mapsto(-1)^{g-1}a_g(1-X)$, mixing the coefficients through an **integer triangular binomial matrix with determinant $\pm1$**. Therefore the least common denominator of the *whole vector* is invariant, although an individual coefficient's denominator may change.

    This is a useful robust denominator pattern. It is **not** a proof of a smallest-integer principle or a theorem about arbitrary gauge transformations.
    """)
    return


@app.cell
def _(obj3):
    from bcov.audit_denominators import coeff_H as all_coeffs_hkq, invert_X, lcd, other_than_five
    _g = 3
    hkq_coeffs = all_coeffs_hkq(obj3)
    small_z_coeffs = invert_X(hkq_coeffs, _g)
    yy_coeffs = [x * 5 ** (3 * _g - 3) for x in hkq_coeffs]
    print('Genus 3 LCD in HKQ:', lcd(hkq_coeffs))
    print('Genus 3 LCD in z:  ', lcd(small_z_coeffs))
    print('Genus 3 LCD in YY: ', lcd(yy_coeffs))
    assert lcd(hkq_coeffs) == lcd(small_z_coeffs)
    assert other_than_five(lcd(hkq_coeffs)) == other_than_five(lcd(yy_coeffs))
    assert hkq_coeffs != small_z_coeffs
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 12. Physical free energies and GW/GV invariants: the independent check

    In HKQ's normalization the *holomorphic limit* is

    $$
    F_g^A(q)=\omega_0^{2g-2}\left(\frac{1-\psi}{5\psi}\right)^{g-1}P_g^{\rm HKQ},
    \qquad q=e^{2\pi it}.
    $$

    The mirror map and period are determined by Picard–Fuchs. The **Gromov–Witten** coefficients $[q^d]F_g^A$ are generally rational; the **Gopakumar–Vafa** integers $n_d^h$ appear after the sine/multicover reorganization. For instance,

    $$[q]F_2^A=\frac{n_1^0}{240}=\frac{2875}{240}=\frac{575}{48},$$

    whereas the *constant-map* contribution is $F_2^{\rm const}=-5/144$. Neither is the YY pure-$X$ constant $a_{2,0}=3125/144$.

    The check below uses the **GV table parsed from the user's Mathematica file**, not a handful of hard-coded sample invariants. It checks coefficients of $q$ through $q^4$ at genera 2 and 3 (the package's comprehensive suite goes through genus 15 and $q^6$).
    """)
    return


@app.cell
def _(ATOMIC, DATA, Q):
    from bcov.verify_physical_qseries import verify as verify_gv_qseries
    GV_FILE = DATA / 'GVTab.csv.gz'
    checks = verify_gv_qseries(ATOMIC, GV_FILE, N=4, gmax=3)
    for result in checks:
        assert all(Q(x)==0 for x in result['GV_residuals'])
        print('g=',result['genus'], 'constant F_g^A=',result['F_constant'],
              ', [q]F_g^A=',result['coefficients'][0])
    assert Q(checks[0]['coefficients'][0]) == Q(575,48)
    assert Q(checks[0]['F_constant']) == Q(-5,144)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 13. Save a boundary-fixed checkpoint and run the regression suite

    A `BCOV` checkpoint saves exactly which ambiguity choices were supplied, along with the convention identifier. `BCOV.load()` rechecks the anomaly recursively. **Do not** copy a checkpoint into a different frame or relabel its coordinates by editing JSON.

    For a full independent audit, set `RUN_FULL_SUITE=True` to execute `run_all_checks.py` (requires SymPy only for the propagator-gauge check). The full suite additionally checks the Bernoulli and GW/GV data through genus 15.
    """)
    return


@app.cell
def _(BCOV, P3, Path, solver):
    OUTPUT = Path.cwd() / 'notebook_output'
    OUTPUT.mkdir(exist_ok=True)
    checkpoint = OUTPUT / 'bcov_genus3_YY.json'
    solver.save(checkpoint)
    reloaded = BCOV.load(checkpoint)
    assert reloaded.amplitudes[3].polynomial() == P3
    print('Saved and reloaded exact YY checkpoint:', checkpoint.name)
    return


@app.cell
def _(PACKAGE, sys):
    RUN_FULL_SUITE = False  # change to True for the full, independent v0.2 audit
    if RUN_FULL_SUITE:
        import subprocess
        completed = subprocess.run([sys.executable, 'run_all_checks.py'],
                                  cwd=PACKAGE.parent, capture_output=True, text=True)
        print(completed.stdout[-3500:])
        if completed.returncode:
            print(completed.stderr)
        assert completed.returncode == 0
    else:
        print('Skipped the optional full-suite rerun; the notebook assertions above passed.')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 14. Diagnostics and choices you should never silently conflate

    **If the $P_2$ comparison fails:** first inspect your frame label and apply `hkq_to_yy(p,2)` (or the corresponding map), then verify the Yukawa-cubic factor and generator scaling. Do **not** refit individual rational coefficients.

    **If your $q$-series differs from GV:** test the Picard–Fuchs coordinate and mirror map, the $\omega_0^{2g-2}$ period factor, the $(1-\psi)/(5\psi)$ multiplier, constant maps, and multi-covers. A GW coefficient need not be an integer.

    **If your pure-$X$ coefficients differ:** check whether you extracted from the *full physical polynomial* or the *zero-at-$v$ particular primitive*. A nonzero polynomial in $v_i$ can contribute to $[X^k]P$ without contributing to the pure-$X$ coefficient. When changing propagator gauge, the separation into anomaly and ambiguity also changes.

    **If you want to go to genus 4:** supply the full **physical** genus-three ambiguity first, as above. Using `ambiguity=()` is a *formal* zero-boundary exercise, not a quintic prediction. The same caution extends to purported all-genus Bernoulli laws.

    **If you want new $q$-series precision:** increase the truncation order in the independent physical verifier and confirm the GV table covers the requested degree; entries marked unknown are **not zeros**. The two high-genus GV tables in the user's source have 29 unresolved disagreements above genus 25; nothing in our genus 2–15 tests adjudicates those cells.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Exercises / suggested experiments

    1. **Another frame:** compute `yy.to('small_z').poly.coefficient('X',1)`. Why is that coefficient *not* simply a sign change of the original YY $[X]P_2$?
    2. **Depth five:** call `pure_X_coeffs_atomic_hkq(load_genus_atomic(15), 5)`. Time the call and compare with the would-be full expansion cost.
    3. **Independent genus four:** reconstruct the genus-four physical pure-$X$ vector, integrate with the genus-three boundary data fixed above, and check exact equality against the HKQ atomic genus-four polynomial.
    4. **Arithmetic:** compare prime-adic denominator valuations of the full ambiguity vectors in the three frames; determine which properties survive unimodular changes of basis.
    5. **Boundary bootstrap:** withhold one physical $a_{g,k}$ as an independent test, rather than fitting it, when exploring Bernoulli/DT/PT/resurgent constraints.

    **Core principle:** the BCOV differential equation, the convention dictionary, and the enumerative boundary data are three separate pieces of information. Keeping them separate is what makes these calculations independently checkable.
    """)
    return


if __name__ == "__main__":
    app.run()

