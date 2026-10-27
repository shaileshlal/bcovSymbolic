# bcovSymbolic toolkit

Exact-arithmetic core of the quintic BCOV recursion, extracted from the
`pypackage` branch, pruned, and fully type-annotated and documented.
The public names the Quintic BCOV tutorial imports are unchanged.

| File | Contents |
|---|---|
| `bcov_yy.py` | Sparse rational ring `Poly` in `u, v1, v2, v3, X`; `YYRing` (Euler derivation `theta`, covariant derivative); lazy `Recipe`/`RecipeArena` DAG; `Amplitude`; `integrate_gradient`; `BCOV` recursion solver with exact JSON checkpoints; `reference_genus2`, `YY_GENUS2_AMBIGUITY`; `demo` CLI |
| `quintic_frames.py` | Frame maps `yy_to_hkq`, `hkq_to_yy`, `hkq_to_z`, `z_to_hkq`, `yy_to_z`, `z_to_yy`; Euler derivations `theta_hkq`, `theta_z`; `GU_GAUGE`; frame-tagged `FramedPolynomial` |
| `test_frames.py`, `CONVENTIONS.json` | Frame-change regression test and the conventions ledger it checks against |
| `tools/callgraph.py` | Builds the call graph below |
| `audit_coefficients.py` | Bernoulli constant and linear pure-`X` coefficients, genera 2–15 → `low_x_reaudit_genus_2_15.csv` |
| `audit_denominators.py` | Denominator-lattice invariance of the pure-`X` polynomials → `denominator_invariants.csv` |
| `verify_physical_qseries.py` | Physical A-model `q`-series vs `GVTab` multicover prediction → `qseries_gv_residuals.json` |
| `verify_yy_qseries_legacy.py` | Older SymPy-based `q`-series check from the raw `ff_g` files (genus 2–3 GV assertions) |
| `data/` | Raw free energies `ff_2..ff_15.txt`, `quinticGVDTPT.m`, extracted `GVTab.csv.gz`, HKQ atomic data `bcov_hkq_atomic_genus_2_15.jsonl.gz`, propagators, GV manifest and conflicts |

The core modules and the audit/verification scripts (except the legacy one) need only the
standard library, Python 3.10+. `verify_yy_qseries_legacy.py` needs SymPy.
The audit and verification scripts still use their own polynomial and series
code rather than `Poly`; that is unchanged for now.

```bash
python test_frames.py
python bcov_yy.py --genus 3
python audit_coefficients.py
python audit_denominators.py --genus 8
python verify_physical_qseries.py --genus 15 --degree 6
python verify_yy_qseries_legacy.py
```

The result files the scripts write are not committed; regenerated outputs are
byte-identical to those on `pypackage`.

## Call graph

Built from two sources. **pyan3** (static) guesses edges by name, so it both invents
calls and misses calls on objects other than `self`. A **runtime trace** of
`test_frames.py`, `bcov_yy.py --genus 4` and the tutorial records the calls that
actually happen. The diagrams show:

- `-->` found by pyan3 **and** observed at runtime
- `-.->` observed at runtime only (pyan3 missed it)
- `==>` real call (checked by reading the code) that the traced runs do not exercise

The solver diagram hides edges into `Poly` arithmetic and construction
(`+ - * / **`, `==`, `coerce`, `constant`, `scale`, `terms`, ...) and into the
helpers `rational`, `variable_index`, `_weight`, `_product`, since nearly every
function uses them; those edges appear in the first diagram.

### `Poly` and helpers

```mermaid
flowchart LR
  subgraph Poly
    bcov_yy_Poly__add_["Poly.__add__"]
    bcov_yy_Poly__bool_["Poly.__bool__"]
    bcov_yy_Poly__eq_["Poly.__eq__"]
    bcov_yy_Poly__init_["Poly.__init__"]
    bcov_yy_Poly__len_["Poly.__len__"]
    bcov_yy_Poly__mul_["Poly.__mul__"]
    bcov_yy_Poly__neg_["Poly.__neg__"]
    bcov_yy_Poly__pow_["Poly.__pow__"]
    bcov_yy_Poly__repr_["Poly.__repr__"]
    bcov_yy_Poly__rsub_["Poly.__rsub__"]
    bcov_yy_Poly__str_["Poly.__str__"]
    bcov_yy_Poly__sub_["Poly.__sub__"]
    bcov_yy_Poly__truediv_["Poly.__truediv__"]
    bcov_yy_Poly_coefficient["Poly.coefficient"]
    bcov_yy_Poly_coerce["Poly.coerce"]
    bcov_yy_Poly_constant["Poly.constant"]
    bcov_yy_Poly_degree["Poly.degree"]
    bcov_yy_Poly_evaluate["Poly.evaluate"]
    bcov_yy_Poly_from_data["Poly.from_data"]
    bcov_yy_Poly_integrate["Poly.integrate"]
    bcov_yy_Poly_partial["Poly.partial"]
    bcov_yy_Poly_scale["Poly.scale"]
    bcov_yy_Poly_substitute["Poly.substitute"]
    bcov_yy_Poly_variable["Poly.variable"]
  end
  subgraph bcov_yy
    bcov_yy_product["_product"]
    bcov_yy_weight["_weight"]
    bcov_yy_rational["rational"]
    bcov_yy_variable_index["variable_index"]
  end
  bcov_yy_Poly__add_ --> bcov_yy_Poly__init_
  bcov_yy_Poly__add_ -.-> bcov_yy_Poly_coerce
  bcov_yy_Poly__eq_ -.-> bcov_yy_Poly_constant
  bcov_yy_Poly__init_ --> bcov_yy_rational
  bcov_yy_Poly__mul_ -.-> bcov_yy_Poly__bool_
  bcov_yy_Poly__mul_ --> bcov_yy_Poly__init_
  bcov_yy_Poly__mul_ -.-> bcov_yy_Poly__len_
  bcov_yy_Poly__mul_ --> bcov_yy_Poly_scale
  bcov_yy_Poly__neg_ --> bcov_yy_Poly_scale
  bcov_yy_Poly__pow_ -.-> bcov_yy_Poly__mul_
  bcov_yy_Poly__pow_ -.-> bcov_yy_Poly_constant
  bcov_yy_Poly__repr_ ==> bcov_yy_Poly_degree
  bcov_yy_Poly__rsub_ -.-> bcov_yy_Poly__sub_
  bcov_yy_Poly__rsub_ -.-> bcov_yy_Poly_coerce
  bcov_yy_Poly__str_ --> bcov_yy_weight
  bcov_yy_Poly__sub_ -.-> bcov_yy_Poly__add_
  bcov_yy_Poly__sub_ -.-> bcov_yy_Poly__neg_
  bcov_yy_Poly__sub_ -.-> bcov_yy_Poly_coerce
  bcov_yy_Poly__truediv_ --> bcov_yy_Poly_scale
  bcov_yy_Poly__truediv_ --> bcov_yy_rational
  bcov_yy_Poly_coefficient --> bcov_yy_Poly__init_
  bcov_yy_Poly_coefficient --> bcov_yy_variable_index
  bcov_yy_Poly_coerce -.-> bcov_yy_Poly_constant
  bcov_yy_Poly_constant -.-> bcov_yy_Poly__init_
  bcov_yy_Poly_constant -.-> bcov_yy_rational
  bcov_yy_Poly_degree --> bcov_yy_weight
  bcov_yy_Poly_degree --> bcov_yy_variable_index
  bcov_yy_Poly_evaluate --> bcov_yy_product
  bcov_yy_Poly_evaluate ==> bcov_yy_rational
  bcov_yy_Poly_evaluate ==> bcov_yy_variable_index
  bcov_yy_Poly_from_data -.-> bcov_yy_Poly__init_
  bcov_yy_Poly_integrate --> bcov_yy_Poly__init_
  bcov_yy_Poly_integrate --> bcov_yy_variable_index
  bcov_yy_Poly_partial --> bcov_yy_Poly__init_
  bcov_yy_Poly_partial --> bcov_yy_variable_index
  bcov_yy_Poly_scale -.-> bcov_yy_Poly__bool_
  bcov_yy_Poly_scale --> bcov_yy_Poly__init_
  bcov_yy_Poly_scale --> bcov_yy_rational
  bcov_yy_Poly_substitute -.-> bcov_yy_Poly__add_
  bcov_yy_Poly_substitute --> bcov_yy_Poly__init_
  bcov_yy_Poly_substitute -.-> bcov_yy_Poly__mul_
  bcov_yy_Poly_substitute -.-> bcov_yy_Poly__pow_
  bcov_yy_Poly_substitute -.-> bcov_yy_Poly_coerce
  bcov_yy_Poly_substitute -.-> bcov_yy_Poly_constant
  bcov_yy_Poly_substitute --> bcov_yy_variable_index
  bcov_yy_Poly_variable -.-> bcov_yy_Poly__init_
  bcov_yy_Poly_variable -.-> bcov_yy_variable_index
```

### Solver and frames

```mermaid
flowchart LR
  subgraph Amplitude
    bcov_yy_Amplitude__init_["Amplitude.__init__"]
    bcov_yy_Amplitude_from_data["Amplitude.from_data"]
    bcov_yy_Amplitude_polynomial["Amplitude.polynomial"]
    bcov_yy_Amplitude_to_data["Amplitude.to_data"]
    bcov_yy_Amplitude_with_ambiguity["Amplitude.with_ambiguity"]
  end
  subgraph BCOV
    bcov_yy_BCOV__init_["BCOV.__init__"]
    bcov_yy_BCOV_anomaly_gradient["BCOV.anomaly_gradient"]
    bcov_yy_BCOV_integrate["BCOV.integrate"]
    bcov_yy_BCOV_load["BCOV.load"]
    bcov_yy_BCOV_residual["BCOV.residual"]
    bcov_yy_BCOV_save["BCOV.save"]
    bcov_yy_BCOV_set_ambiguity["BCOV.set_ambiguity"]
    bcov_yy_BCOV_source_recipe["BCOV.source_recipe"]
    bcov_yy_BCOV_vertex["BCOV.vertex"]
  end
  subgraph FramedPolynomial
    quintic_frames_FramedPolynomial__post_init_["FramedPolynomial.__post_init__"]
    quintic_frames_FramedPolynomial_to["FramedPolynomial.to"]
  end
  subgraph Poly
    bcov_yy_Poly__str_["Poly.__str__"]
    bcov_yy_Poly_coefficient["Poly.coefficient"]
    bcov_yy_Poly_degree["Poly.degree"]
    bcov_yy_Poly_from_data["Poly.from_data"]
    bcov_yy_Poly_integrate["Poly.integrate"]
    bcov_yy_Poly_partial["Poly.partial"]
    bcov_yy_Poly_substitute["Poly.substitute"]
    bcov_yy_Poly_to_data["Poly.to_data"]
  end
  subgraph Recipe
    bcov_yy_Recipe__add_["Recipe.__add__"]
    bcov_yy_Recipe__init_["Recipe.__init__"]
    bcov_yy_Recipe__mul_["Recipe.__mul__"]
    bcov_yy_Recipe_expand["Recipe.expand"]
  end
  subgraph RecipeArena
    bcov_yy_RecipeArena__init_["RecipeArena.__init__"]
    bcov_yy_RecipeArena_make["RecipeArena._make"]
    bcov_yy_RecipeArena_binary["RecipeArena.binary"]
    bcov_yy_RecipeArena_clear_values["RecipeArena.clear_values"]
    bcov_yy_RecipeArena_coerce["RecipeArena.coerce"]
    bcov_yy_RecipeArena_covariant["RecipeArena.covariant"]
    bcov_yy_RecipeArena_expand["RecipeArena.expand"]
    bcov_yy_RecipeArena_leaf["RecipeArena.leaf"]
    bcov_yy_RecipeArena_stats["RecipeArena.stats"]
  end
  subgraph YYRing
    bcov_yy_YYRing__init_["YYRing.__init__"]
    bcov_yy_YYRing_monomial_derivative["YYRing._monomial_derivative"]
    bcov_yy_YYRing_covariant["YYRing.covariant"]
    bcov_yy_YYRing_theta["YYRing.theta"]
  end
  subgraph bcov_yy
    bcov_yy_demo["demo"]
    bcov_yy_integrate_gradient["integrate_gradient"]
    bcov_yy_reference_genus2["reference_genus2"]
  end
  subgraph quintic_frames
    quintic_frames_frame["_frame"]
    quintic_frames_genus["_genus"]
    quintic_frames_hkq_to_yy["hkq_to_yy"]
    quintic_frames_hkq_to_z["hkq_to_z"]
    quintic_frames_theta_hkq["theta_hkq"]
    quintic_frames_theta_z["theta_z"]
    quintic_frames_yy_to_hkq["yy_to_hkq"]
    quintic_frames_yy_to_z["yy_to_z"]
    quintic_frames_z_to_hkq["z_to_hkq"]
    quintic_frames_z_to_yy["z_to_yy"]
  end
  bcov_yy_Amplitude__init_ -.-> bcov_yy_Poly_degree
  bcov_yy_Amplitude__init_ -.-> bcov_yy_Poly_substitute
  bcov_yy_Amplitude_from_data -.-> bcov_yy_Amplitude__init_
  bcov_yy_Amplitude_from_data -.-> bcov_yy_Poly_from_data
  bcov_yy_Amplitude_to_data -.-> bcov_yy_Poly_to_data
  bcov_yy_Amplitude_with_ambiguity --> bcov_yy_Amplitude__init_
  bcov_yy_BCOV__init_ --> bcov_yy_RecipeArena__init_
  bcov_yy_BCOV__init_ --> bcov_yy_YYRing__init_
  bcov_yy_BCOV_anomaly_gradient --> bcov_yy_BCOV_source_recipe
  bcov_yy_BCOV_anomaly_gradient -.-> bcov_yy_Poly_coefficient
  bcov_yy_BCOV_anomaly_gradient -.-> bcov_yy_Poly_degree
  bcov_yy_BCOV_anomaly_gradient -.-> bcov_yy_Recipe_expand
  bcov_yy_BCOV_integrate --> bcov_yy_Amplitude__init_
  bcov_yy_BCOV_integrate --> bcov_yy_BCOV_anomaly_gradient
  bcov_yy_BCOV_integrate --> bcov_yy_integrate_gradient
  bcov_yy_BCOV_load -.-> bcov_yy_Amplitude_from_data
  bcov_yy_BCOV_load -.-> bcov_yy_BCOV__init_
  bcov_yy_BCOV_load -.-> bcov_yy_BCOV_integrate
  bcov_yy_BCOV_residual -.-> bcov_yy_Amplitude_polynomial
  bcov_yy_BCOV_residual --> bcov_yy_BCOV_source_recipe
  bcov_yy_BCOV_residual -.-> bcov_yy_Poly_partial
  bcov_yy_BCOV_residual -.-> bcov_yy_Recipe_expand
  bcov_yy_BCOV_save -.-> bcov_yy_Amplitude_to_data
  bcov_yy_BCOV_set_ambiguity -.-> bcov_yy_Amplitude_with_ambiguity
  bcov_yy_BCOV_set_ambiguity --> bcov_yy_BCOV_integrate
  bcov_yy_BCOV_set_ambiguity --> bcov_yy_RecipeArena_clear_values
  bcov_yy_BCOV_source_recipe --> bcov_yy_BCOV_vertex
  bcov_yy_BCOV_source_recipe -.-> bcov_yy_Recipe__add_
  bcov_yy_BCOV_source_recipe -.-> bcov_yy_Recipe__mul_
  bcov_yy_BCOV_vertex -.-> bcov_yy_Amplitude_polynomial
  bcov_yy_BCOV_vertex --> bcov_yy_BCOV_vertex
  bcov_yy_BCOV_vertex --> bcov_yy_RecipeArena_covariant
  bcov_yy_BCOV_vertex --> bcov_yy_RecipeArena_leaf
  bcov_yy_Recipe__add_ --> bcov_yy_RecipeArena_binary
  bcov_yy_Recipe__mul_ --> bcov_yy_RecipeArena_binary
  bcov_yy_Recipe_expand --> bcov_yy_RecipeArena_expand
  bcov_yy_RecipeArena_make --> bcov_yy_Recipe__init_
  bcov_yy_RecipeArena_binary --> bcov_yy_RecipeArena_make
  bcov_yy_RecipeArena_binary --> bcov_yy_RecipeArena_coerce
  bcov_yy_RecipeArena_coerce --> bcov_yy_RecipeArena_leaf
  bcov_yy_RecipeArena_covariant --> bcov_yy_RecipeArena_make
  bcov_yy_RecipeArena_covariant --> bcov_yy_RecipeArena_coerce
  bcov_yy_RecipeArena_expand --> bcov_yy_RecipeArena_coerce
  bcov_yy_RecipeArena_expand --> bcov_yy_YYRing_covariant
  bcov_yy_RecipeArena_leaf --> bcov_yy_RecipeArena_make
  bcov_yy_YYRing_covariant --> bcov_yy_YYRing_theta
  bcov_yy_YYRing_theta -.-> bcov_yy_YYRing_monomial_derivative
  bcov_yy_demo -.-> bcov_yy_Amplitude_polynomial
  bcov_yy_demo --> bcov_yy_BCOV__init_
  bcov_yy_demo --> bcov_yy_BCOV_integrate
  bcov_yy_demo --> bcov_yy_BCOV_residual
  bcov_yy_demo ==> bcov_yy_BCOV_save
  bcov_yy_demo --> bcov_yy_BCOV_set_ambiguity
  bcov_yy_demo --> bcov_yy_BCOV_vertex
  bcov_yy_demo -.-> bcov_yy_Poly__str_
  bcov_yy_demo -.-> bcov_yy_Poly_degree
  bcov_yy_demo -.-> bcov_yy_Poly_partial
  bcov_yy_demo -.-> bcov_yy_Recipe_expand
  bcov_yy_demo --> bcov_yy_RecipeArena_stats
  bcov_yy_demo --> bcov_yy_reference_genus2
  bcov_yy_integrate_gradient --> bcov_yy_Poly_degree
  bcov_yy_integrate_gradient -.-> bcov_yy_Poly_integrate
  bcov_yy_integrate_gradient --> bcov_yy_Poly_partial
  quintic_frames_FramedPolynomial__post_init_ --> quintic_frames_frame
  quintic_frames_FramedPolynomial__post_init_ --> quintic_frames_genus
  quintic_frames_FramedPolynomial_to --> quintic_frames_frame
  quintic_frames_FramedPolynomial_to ==> quintic_frames_hkq_to_yy
  quintic_frames_FramedPolynomial_to -.-> quintic_frames_hkq_to_z
  quintic_frames_FramedPolynomial_to -.-> quintic_frames_yy_to_hkq
  quintic_frames_FramedPolynomial_to -.-> quintic_frames_yy_to_z
  quintic_frames_FramedPolynomial_to ==> quintic_frames_z_to_hkq
  quintic_frames_FramedPolynomial_to -.-> quintic_frames_z_to_yy
  quintic_frames_hkq_to_yy -.-> bcov_yy_Poly_substitute
  quintic_frames_hkq_to_yy --> quintic_frames_genus
  quintic_frames_hkq_to_z -.-> bcov_yy_Poly_substitute
  quintic_frames_hkq_to_z --> quintic_frames_genus
  quintic_frames_theta_hkq -.-> bcov_yy_Poly_substitute
  quintic_frames_theta_hkq --> bcov_yy_YYRing__init_
  quintic_frames_theta_hkq -.-> bcov_yy_YYRing_theta
  quintic_frames_theta_z -.-> bcov_yy_Poly_substitute
  quintic_frames_theta_z --> quintic_frames_theta_hkq
  quintic_frames_yy_to_hkq -.-> bcov_yy_Poly_substitute
  quintic_frames_yy_to_hkq --> quintic_frames_genus
  quintic_frames_yy_to_z --> quintic_frames_hkq_to_z
  quintic_frames_yy_to_z --> quintic_frames_yy_to_hkq
  quintic_frames_z_to_hkq --> quintic_frames_hkq_to_z
  quintic_frames_z_to_yy --> quintic_frames_hkq_to_yy
  quintic_frames_z_to_yy --> quintic_frames_z_to_hkq
```

pyan3 edges rejected as false positives:

- `BCOV.save` → `BCOV.__init__`
- `Recipe.__add__` → `Recipe.__init__`
- `Recipe.__mul__` → `Recipe.__init__`
- `Recipe.expand` → `Recipe.__init__`
- `RecipeArena._make` → `RecipeArena.__init__`
- `RecipeArena.clear_values` → `RecipeArena.__init__`
- `YYRing.__init__` → `YYRing._monomial_derivative`
- `demo` → `integrate_gradient`

To regenerate (the tutorial run needs the `pypackage` data files, so run from a
checkout of that branch with these modules copied in):

```bash
pyan3 bcov_yy.py quintic_frames.py --uses --no-defines --dot > /tmp/cg.dot
python tools/callgraph.py /tmp/cg.dot /tmp/callgraph.md \
    test_frames.py "bcov_yy.py --genus 4" Quintic_BCOV_Tutorial.py
```

## Changes from `pypackage`

Removed (unused by the tutorial, tests, demo, or any other script on `pypackage`):

- `Poly.homogeneous_part`
- `YYRing.theta_z`: unused, and easily confused with the unrelated `quintic_frames.theta_z`
- `YYRing.clear_cache` and `BCOV.compact_recipes`
- `Recipe.__neg__`

Refactored without changing behaviour:

- Every function and method has type hints and a docstring; type aliases `Scalar`, `Variable`, `RecipeSpec`.
- `bcov_yy._weight` replaces three copies of the weighted-degree sum.
- `quintic_frames`: the generator substitutions are named constants
  (`_HKQ_IN_YY`, `_YY_IN_HKQ`, `_HKQ_Z_SWAP`); the frame-label check is `_frame`;
  the conversion table `_CONVERSIONS` is built once at import; `FRAMES` is derived from `Frame`.

Kept as is: all input-validation guards, and every name the tutorial imports.

The audit and verification scripts and `data/` are copied unchanged, except that
`verify_yy_qseries_legacy.py` now reads `data/ff_{g}.txt` next to the script
instead of the hard-coded `/mnt/data/ff_{g}.txt`.

Verified against `pypackage` with these modules dropped in: `test_frames.py` passes,
the `bcov_yy.py --genus 4` checkpoint is byte-identical, the tutorial output is
identical, and `run_all_checks.py` passes.
