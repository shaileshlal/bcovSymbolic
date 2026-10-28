# GV lattice walking on the residual quintic YY ambiguity

**27 October 2026 — exact rational computational audit**

## Question

Does the minimum-absolute-GV representative of the remaining x(1-x) ambiguity agree with the physical reference? If it does not, can an integer lattice translation be inferred using the section representation and consistency tests?

## Conventions and independent information

We use x=psi^(-5), X=x/(x-1), K=3g-3, m=ceil(K/5), D=floor((2g-2)/5), and a_g=a_g^part+δa_g. The homogeneous boundary kernel is

δa_g(X)=X(1-X)^m Σ_{r=0}^{D-1} c_{g,r} X^r.

The reduced section is

δf_g(x)=(x-1)^(g-1)δa_g(x/(x-1))
       =Σ_{r=0}^{D-1} (-1)^(g-2-r)c_{g,r} x^(r+1)(1-x)^(D-1-r).

For g=4,...,10, `experiment.py` reconstructs the zero-free-coordinate affine particular directly from the BCOV anomaly polynomial and the conifold gap, orbifold regularity, and MUM **constant map only**; it does not use the target-genus physical pure-X coefficients to build that particular. The reduced jet precision n=3g+8 passed exact reference consistency checks at those genera; full-conservative-precision evaluation was separately checked through genus 8 in the preceding experiment. For g=11,12 the particular is extracted from the published/reference polynomial (retrospective only). The exact MUM response matrices use the physical lower-genus amplitudes, so this is not an independent sequential high-genus bootstrap.

## Exact triangular lattice map

For d=1,...,D, in the standard GV normalization,

    n_{g,d} = n^{part}_{g,d} + Σ_{r=0}^{D-1} H_{g;d,r} c_{g,r},

with H_{g;d,r}=0 if r>=d and

    H_{g;d,d-1}=(-1)^(g-1+d) 5^(5d-4g+4).

These assertions were checked against exact period/mirror-map series and the full GV multi-cover subtraction. Consequently, given an independently determined BPS charge displacement Δn in Z^D, the complete section displacement is

    Δc=H_g^(-1)Δn,
    Δf_g(x)=Σ_r (-1)^(g-2-r)Δc_r x^(r+1)(1-x)^(D-1-r).

This map is fully calculated in `results/unit_charge_steps.json`.

## Minimum-height test

All 22 residual degrees in the corrected g=4,..,12 reference set correspond to geometrically vanishing GV charges. Minimize |n_{g,d}| independently at d=1,...,D, in increasing d. Choosing zero and then solving the exact triangular equations reproduces every reference residual coefficient:

| g | D | reference recovery | independence of particular |
|---:|---:|---|---|
| 4 | 1 | 1/1 | independently boundary solved |
| 5 | 1 | 1/1 | independently boundary solved |
| 6 | 2 | 2/2 | independently boundary solved |
| 7 | 2 | 2/2 | independently boundary solved |
| 8 | 2 | 2/2 | independently boundary solved |
| 9 | 3 | 3/3 | independently boundary solved (reduced jet bound) |
| 10 | 3 | 3/3 | independently boundary solved (reduced jet bound) |
| 11 | 4 | 4/4 | *reference-derived anchor, retrospective* |
| 12 | 4 | 4/4 | *reference-derived anchor, retrospective* |

Independent results: **14/14** coefficients in g=4,..,10. Retrospective: 8/8 in g=11,12. An oracle denominator-cap congruence test also selects all 22/22, but the cap is the least common multiple of the target reference denominators and must NOT be used as independent predictive evidence.

**Interpretation:** Because all these degrees obey known geometric vanishing, this is a strong exact implementation test but is not evidence that minimum height identifies genuinely nonzero charges. In the unconstrained rational coefficient space, any integral value of the first D GV charges is attainable. Minimum height then simply chooses the zero-charge vector whenever allowed.

## Integer transport beyond the D unconstrained coefficients

For each genus, compute H_{g;d,r} out to d=26 and define the response to the kth **unit change of low-degree GV charge**:

    A_{g;d,k} = Σ_r H_{g;d,r}(H_{g;1:D}^{-1})_{r,k}.

All **572 entries** of these matrices across g=4,..,12 and d=1,..,26 are exactly integral. For g=4,6,9,12 the same holds through degree 50 (500 entries). This is an observed finite range, not an all-degree proof.

Thus, in our tests, imposing **integrality of additional GV invariants at higher degrees does not shrink the unit-charge lattice**. One cannot determine the physical displacement just by checking that the transformed x(1-x) polynomial is a valid section or that further GV numbers are integers.

### Constructive counterexample: genus four

The physical residual coefficient is

    c_{4,0} = -101510017565 / 217728.

Change it by Δc_{4,0}=78125. This produces

    Δf_4(x)=78125 x,
    Δn_{4,1}=+1.

The difference in later GV invariants begins

    (Δn_{4,1},Δn_{4,2},Δn_{4,3},Δn_{4,4},Δn_{4,5})
      = (1, -82, -41118, -26613400, -26396534375).

These are integers. Both amplitudes obey the same conifold gap, orbifold conditions and MUM constant map. The physical genus-four choice is selected because n_{4,1}=0 is known geometrically, **not** because the lattice-shifted section fails a consistency condition.

### Constructive counterexample: genus six

Keeping Δn_{6,2}=0 and changing Δn_{6,1}=1 gives

    (Δc_{6,0},Δc_{6,1})=(30517578125,-31318359375),

or equivalently

    Δf_6(x)=30517578125 x(1-x)+31318359375 x^2.

The higher GV response begins (1,0,57962,34162756,25523186625), again entirely integral.

## Comparing the old naive-ML estimates with lattice decoding

Taking the existing g=9–12 model predictions of the residual coefficients (trained only through genus 8) and rounding **sequentially** to the nearest integral GV charges does **not** generally recover the physical polynomial. For the best local-source magnitude + majority-sign baseline (`mixed_plus_genus__majority`), only **6/14** held-out coefficients agree with their references; **2/6** of these are in independently solved genera 9 and 10. Example predicted GV charge vectors are:

- g=9: (0,-1,3400) instead of (0,0,0).
- g=10: (0,-1,3231) instead of (0,0,0).

The triangular GV response is substantially sheared: a small error in a lower-index coefficient can produce a large spurious charge at a higher index. Adding independent low-degree geometric vanishing conditions fixes those errors, but doing so is an analytic constraint, not an ML success.

## What can and cannot fix the lattice displacement

If an independent numerical or analytic computation supplies the physical x(1-x) polynomial, the displacement is uniquely calculable:

    Δn=H(c_target-c_min).

If the independent input supplies the first D midpoint jets, recover c_target by polynomial interpolation first. Either method uses **new physical information**. Algebraic recomposition, coordinate inversion at x=1/2, and GV integrality alone add no new equation for the displacement in the cases checked.

For a genuinely new test, pursue g=51, d=20 (known nonvanishing n_{51,20}=175) as a benchmark and then the unprotected g=54, d=21 sector. This archive contains only full corrected amplitudes to g=12 and low-X jets to g=1000; it does not contain independently validated full affine systems at genus 51 or 54. Hence these two cases have **not** been run, and no blind nonzero lattice step has been predicted.

## Execution

From this package directory:

```
python experiment.py
python charge_integrality_and_ml.py
python tests.py
```

`experiment.py` needs the bundled `Quintic_Intermediate_ML/source` data and exact BCOV engine, and its corrected reference inputs. The output audit records every rational coefficient as an exact numerator/denominator string. `charge_integrality_and_ml.py` reuses the earlier learned-model predictions only for held-out projection diagnostics. No new neural-network training is claimed.

### Explicit correction vectors for the held-out ML model

For the `train_to_8 / mixed_plus_genus__majority` predictor, the known geometric vanishing changes the provisional charges as follows:

| Genus | Predicted GV charges | Exact correction Δn | Induced Δc in kernel basis |
|---:|---|---|---|
| 9 | (0, -1, 3400) | (0, +1, -3400) | (0, 2384185791015625, 499725341796875) |
| 10 | (0, -1, 3231) | (0, +1, -3231) | (0, -1490116119384765625, -346183776855468750) |

These are **retrospective corrections supplied by an independent vanishing theorem**, not offsets learned by the model. They illustrate why tracking integer charge vectors is preferable to fitting huge rational coefficient differences.
