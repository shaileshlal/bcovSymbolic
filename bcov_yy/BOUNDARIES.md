# Physical quintic boundaries and the genus-53 driver

This extension solves the *physical* ambiguity in the existing engine's
zero-at-v=0 convention. It adds exact local period series, the conifold gap,
orbifold regularity, MUM constant maps, GV multicover conversion, and the
Castelnuovo data required to close the recursion through genus 53.
It does not import a table of higher-genus ambiguities.

## Run and resume

Python >=3.10; standard library only:

```bash
python -m unittest discover -s tests -v
python run_physical.py --plan --genus 53
python run_physical.py --genus 6 --gv-degree 4 --output my_physical_run
python run_physical.py --genus 53 --mu mu --output my_physical_run
```

The last command resumes verified lower-genus data, recomputes their needed GV
series through degree 20, and targets genus 53. Default per-operation limits are
1,000,000 polynomial terms and 50,000,000 pair products. They stop the run with
an explicit error and retain the last committed physical genus. To change them:

```bash
python run_physical.py --genus 53 --max-terms 3000000 --max-products 200000000 --output my_physical_run
```

A zero limit disables that particular limit. These are operation budgets, not
an operating-system memory cap; external process termination or power failure
can still interrupt a run. Amplitude checkpoints use atomic file replacement.
A run never substitutes zero ambiguities or continues past a rank defect.
Changing a lower-genus physical convention requires a fresh output directory.

**Genus 53 is a configured mathematical target, not a demonstrated performance
result.** The expanded sparse-polynomial backend may be impractical there.
No full genus-53 run or high-genus memory/runtime estimate is claimed.
The supplied `validated_physical` directory records exactly what was executed,
including completed genus, boundary ranks, coefficient bit sizes, timings, and
GV integrality checks. Its `plan.json` describes that validation run. The
separate `genus53_plan.json` describes the full target without pretending it ran.

## What keeping mu free means

Set x=z/mu, with mu nonzero, so

```
L = theta_x^4 - x product_(j=1..4)(theta_x+j/5)
X = x/(x-1),   delta=1-x,   psi^(-5)=x.
```

All period and YY ring arithmetic is over QQ in x, with no numerical mu.
The parameter enters only when converting x back to z. The original quintic
coordinate is z_can=x/3125. Its physical instanton variable is
q_can=x exp(h/f)/3125. In an arbitrary z coordinate,
q_can=z exp(h/f)/(3125 mu). Thus changing mu does not change the enumerative
quintic, its Picard-Fuchs exponents, or its intersection numbers.

`--mu` records this external scale, with `mu` as the default symbolic label.
It is deliberately not substituted into rational coefficients. For a jet f(x),
`period_frames.coefficients_in_z(f)` returns triples (n,c_n,-n), meaning
c_n z^n mu^(-n). No general rational-function coefficient field is necessary.

If instead one defines q_z=mu x exp(h/f), its degree-d coefficient is the
canonical coefficient times (3125 mu)^(-d). At mu=1/3125 it coincides with q_can.
Do not impose integer canonical GV invariants in q_z for arbitrary mu.
The boundary equations in this project always use q_can.

Pulling out mu removes artificial powers caused by a coordinate choice.
It **does not** remove Bernoulli denominators, enumerative growth, or intrinsic
rational coefficient growth in P_g. Both sparse-polynomial and Laurent-series
multiplication now collect common denominators and multiply integers internally,
reducing rationals only once per output coefficient.

## Local monodromy and frames

`monodromy.py` provides three exact *local* representations:

* MUM: a maximally unipotent matrix, acting on Frobenius derivatives without
  factorial division, with entries binomial(i,j) below the diagonal.
* Conifold: nu stays invariant and its logarithmic companion changes by +nu.
  Its nilpotent part has rank one and square zero.
* Orbifold: the companion matrix of t^4+t^3+t^2+t+1, of order five.

These matrices use different local bases. They are not global integral
monodromy matrices in a shared symplectic basis. No global product relation or
analytic-continuation connection matrix is asserted. Local monodromy determines
the type of boundary expansion; polarization and topology also enter its
normalization. Monodromy alone does not fix a free energy.

All local holomorphic limits use

```
u  = theta(omega)/omega
v1 = theta(theta(t))/theta(t) + 2u
v2 = theta^2(omega)/omega - u v1
v3 = theta^3(omega)/omega + u v2 - 5uX(v1-2).
```

The evaluated amplitude in a period frame is

```
F_g^frame = [omega^2 (x-1)/625]^(g-1) P_g^hol.
```

This uses the original YY quintic psi gauge, with C=(2pi i)^3(X-1)/125.
It must not be combined with a differently scaled v3 or a z-gauge propagator.

### MUM

f=sum c_n x^n, c_0=1, with
c_n/c_(n-1)=product(n-1+j/5)/n^4.
The logarithmic solution is f log x+h, with
h_n=c_n sum_(k=1..n)[sum_j 1/(k-1+j/5)-4/k].
Here t=log x+h/f and theta=-5x d/dx.
The sign of t is immaterial for the logarithmic derivative defining v1.

The constant-map condition is

```
F_g^A(q=0) = (-1)^g chi |B_(2g) B_(2g-2)| /
             [2 (2g) (2g-2) (2g-2)!],  chi=-200.
```

It gives F_2=-5/144 and F_3=5/36288. Since the frame factor tends to -1/625,
the raw P_2 constant is +3125/144.

### Conifold

In delta=1-x, theta=5(1-delta)d/delta. The analytic period f_c has initial
coefficients (1,0,0). The vanishing period is

```
nu = delta + 7/10 delta^2 + 41/75 delta^3 + ... .
tau = nu/f_c.
```

The coefficient of delta^2 is obtained from the obstruction to an analytic
completion of nu log(delta), not chosen by matching genus two.
`conifold_action` differentiates the indicial action with respect to the
exponent using exact dual arithmetic. The PF recurrence then generates both
periods. The other analytic initial conditions amount to an allowed local
A-period frame choice.

The canonical genus-zero Yukawa is 5/[x^3(1-x)]. In the f_c frame its residue
in the tau coordinate is -5/tau. The unit conifold coordinate is therefore
sqrt(5) tau in the rational frame used here. With the code's coupling convention,

```
F_g^c(tau) = 5^(1-g) (-1)^(g-1) B_(2g) /
             [2g(2g-2) tau^(2g-2)] + O(tau^0).
```

The square root never enters the arithmetic, since 2g-2 is even.
The coefficient is +1/1200 at g=2 and +1/25200 at g=3.
Subtract this full singular series in delta and set every negative delta
coefficient to zero. Keeping only its leading delta pole would lose the gap.

### Orbifold

Let y=psi^5=1/x. Two local periods are
omega_r=psi^r sum c_(r,k)y^k, r=1,2, with c_(r,0)=1 and

```
c_(r,k)/c_(r,k-1) = (r+5k-5)^4 /
                   product_(j=1..4)(r+5k-j).
```

Use omega=omega_1 and t=omega_2/omega_1. Their explicit psi powers are factored
out, so all stored jets have integer y powers. The amplitude is
psi^(-3g+3) times an analytic series S_g(y). Regularity imposes
[y^k]S_g=0 for k < ceil((3g-3)/5).
This is imposed on the full amplitude, not blindly on the pure-X vector.

## Why the boundary data close through genus 53

There are 3g-2 raw pure-X unknowns. At each genus the solver imposes

| Condition | Number |
|---|---:|
| Orbifold negative-power cancellation | ceil((3g-3)/5) |
| Full conifold gap | 2g-2 |
| MUM constant plus low-degree coefficients | 1+floor((2g-2)/5) |

The counts add to 3g-2. Exact elimination also checks actual rank, rather than
assuming that matching counts imply independence.

The positive-degree MUM targets use the multicover formula

```
F(lambda,q) = sum_(h,d,r>=1) n_h^d/r * (2 sin(r lambda/2))^(2h-2) q^(rd),
```

with h starting at zero. Castelnuovo gives n_g^d=0 when
10g>d^2+5d+10. For the required degrees d<=floor((2g-2)/5) and g<=53,
this supplies every needed same-genus value except (g,d)=(51,20).
The extremal value n_51^20=175 is explicitly supplied there.
In particular n_52^20=n_53^20=0; setting n_51^20=0 would be wrong.
These statements are from Liu--Ruan, arXiv:2210.13411v2, Theorems 1.1, 1.3,
and Appendix A. The package implements these constraints in its own YY basis.

The lower-genus GV tables are obtained from computed amplitudes, with genus-zero
Yukawa and genus-one seed as initial data. Lagrange inversion extracts coefficients
of q directly without building z(q). Integral GV values, boundary rank, the full
BCOV residual, and available Castelnuovo values are checked before checkpointing.
At genus 53 the raw ambiguity has 157 slots: 32 orbifold conditions, 104
conifold conditions, and 21 MUM conditions. Only curve degrees through 20 are
needed. This is independent of the cost of constructing the anomaly polynomial.

## Modules and further scaling

* `exact_series.py`: exact Laurent jets with certified absolute precision and
  nested Horner evaluation. Reading an uncertified coefficient raises.
* `period_frames.py`: rational local periods and generators; symbolic mu restore.
* `physical_boundaries.py`: exact rank-aware linear system, gap/constant maps,
  GV multicover kernels, seed and higher-genus GV conversion.
* `monodromy.py`: local representations and nilpotence/order checks.
* `run_physical.py`: physical orchestration, per-genus diagnostics, checkpoints,
  target plan, operation budgets, and resume.
* `BOUNDARIES.tex`: a self-contained mathematical companion.

The main remaining engineering bottlenecks are expanded multivariate products,
repeated frame evaluation, and large rational reductions. A FLINT/Sage backend,
finite-field computation with certified rational reconstruction, or a recursion
that operates directly on boundary jets would improve those costs. None is
silently substituted by this version. Retaining mu symbolically alone does not
solve the high-genus scaling problem.

Sources:
https://arxiv.org/abs/hep-th/0406078
https://arxiv.org/abs/hep-th/0612125
https://arxiv.org/abs/2210.13411
