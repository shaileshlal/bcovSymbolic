# bcovSymbolic

An installable, exact-arithmetic toolkit for quintic BCOV recursion, coordinate
frames, source-amplitude analysis, and GW/GV verification. Distribution name:
`bcov-symbolic`; Python import: `bcov`. Python 3.10 or newer is required.
Core calculations use `fractions.Fraction` and have no external dependencies.

## Install

```bash
conda activate sage-env                       # optional; ordinary Python also works
python -m pip install -e '.[notebooks,dev]'
```

Sage is available for future symbolic extensions but is not required by the
current algorithms. The `notebooks` extra installs marimo; `dev` installs the
optional validation tools and SymPy needed by the historical gauge check.
For a minimal offline installation with setuptools already installed:

```bash
python -m pip install --no-build-isolation --no-deps -e .
```

## Solve the anomaly with explicit boundary data

```python
from bcov import BCOV, YY_GENUS2_AMBIGUITY, reference_genus2

solver = BCOV()
solver.integrate(2)
solver.set_ambiguity(2, YY_GENUS2_AMBIGUITY)
assert solver.amplitudes[2].polynomial() == reference_genus2()

amplitude3 = solver.integrate(3)
assert solver.residual(3, amplitude3.anomaly) == 0
# amplitude3.polynomial() raises: its holomorphic ambiguity is unresolved.
solver.save("checkpoint.json")
restored = BCOV.load("checkpoint.json")
```

`None` means unresolved ambiguity; `()` explicitly chooses all coefficients to
be zero. A formal zero choice is not a physical boundary condition. Checkpoints
preserve exact rational coefficients and unresolved status.

## Work in explicit frames

```python
from bcov import FramedPolynomial, reference_genus2

p = FramedPolynomial(reference_genus2(), genus=2, frame="yy")
assert p.to("hkq").to("small_z").to("yy") == p
```

Frame labels are `yy`, `hkq`, and `small_z`. The wrapper rejects mixed-frame or
mixed-genus addition. Transformations include generator substitutions **and**
the genus-dependent Yukawa normalization. YY/HKQ share the discriminant
coordinate; small-z uses its complement. See
[the conventions ledger](bcov/CONVENTIONS.json) and
[the derivation notes](docs/notes/bcov_yy_conventions_and_derivations_v2.tex).

## Extract and audit arithmetic patterns

```python
from bcov.datasets import load_amplitude
from bcov.patterns import audit_patterns, pure_x, low_x_predictions

source = load_amplitude(3)
yy = pure_x(source, depth=6, frame="yy")  # inclusive X^0,...,X^6
report = audit_patterns(source, depth=4)
assert report.low_x_matches and report.denominator_matches
assert yy[:2] == low_x_predictions(3)[3]
print(report.ratios)  # YY/HKQ ratios to the constant coefficient
```

Extraction sets all propagators to zero and uses the atomic valuation `n+ell`
to avoid expanding the multivariate polynomial. `depth=None` returns the whole
pure-X vector. Small-z extraction transforms the **full** HKQ vector before
truncation: inversion and truncation do not commute. Use `source.expand("yy")`
only at small genus when a full polynomial is needed.

| Pattern or check | API | Evidence / limitation |
|---|---|---|
| Sparse YY differential ring, anomaly integration and recipe DAG | `BCOV`, `Poly`, `YYRing`, `RecipeArena` | Exact algebra; physical ambiguities remain external inputs |
| YY/HKQ/small-z changes and Gu gauge identities | `FramedPolynomial`, `bcov.quintic_frames`; historical gauge test | Whole-amplitude transformations and convention ledger |
| Bernoulli absolute-value convolution | `bernoulli_abs` | Exact rational recurrence |
| Constant and linear pure-X formulas | `low_x_predictions`, `audit_patterns` | Compared with all 14 supplied genera 2–15; not a formula for all ambiguity coefficients |
| Pure-X jets and normalized ratios | `pure_x`, `PatternReport.ratios` | Source-derived, frame-aware; arbitrary requested depth |
| HKQ/small-z denominator lattice; YY agreement away from 5 | `invert_x`, `denominator_lcm`, `audit_patterns` | Full coefficient vectors; integral inverse binomial transformation |
| Periods, mirror map and physical multicover residuals | `bcov.verify_physical_qseries` | Default 84 comparisons: genera 2–15, degrees 1–6 |

The separate **Develop AI Quintic Framework** chat contains newer shelf-depth,
scale-change-pulse, and residual-recurrence work. Its generated code and formulas
were not accessible through the chat reader during this reorganization. Those
extensions are **pending source import**, not implemented or claimed verified
here. The older notes mention low-X checks through genus 53; this repository's
bundled source amplitudes stop at genus 15.

## Command line

Commands work from any directory after installation. Reports are written only
when an output path is supplied.

```bash
bcov recursion --genus 3 --checkpoint checkpoint.json
bcov patterns --genus 15 --depth 4 --output patterns.json
bcov verify --genus 15 --degree 6 --output residuals.json
# Equivalent: python -m bcov ...
```

Pattern audits check the complete denominator vector even when displaying a
short jet. Verification rejects missing or duplicate requested genera and
missing GV cells; unknown entries are never silently substituted with zero.
A nonzero residual produces a nonzero command exit status. Constant maps are
reported separately, not compared to positive-degree GV invariants.

## Marimo demonstrations

```bash
marimo edit Quintic_BCOV_Tutorial.py
marimo edit notebooks/arithmetic_patterns.py
marimo edit notebooks/physical_verification.py
```

The original guided laboratory covers recursion, frames, source extraction and
checkpoints. The arithmetic notebook adds interactive genus/depth controls and
exact coefficient tables; the physical notebook compares GW/GV coefficients.
Run a notebook headlessly with `python path/to/notebook.py`, or export it with
`marimo export html path/to/notebook.py -o /tmp/demo.html`.
The tutorial writes its generated reports to `notebook_output/` in the working
directory; package data are never modified by the new APIs or demonstrations.

## Validation and development

```bash
python -m unittest discover -s tests -v     # dependency-free package regressions
python run_all_checks.py                   # historical full suite; needs SymPy
python -m bcov patterns --genus 15          # low-X and full denominator audits
python -m bcov verify --genus 15 --degree 6 # complete physical comparison
python -m mypy                            # typed public analysis modules
python -m build                          # wheel and source distribution
```

Tests compare atomic jets with independently expanded polynomials in every
frame, exercise unresolved ambiguities and checkpoint round trips, check the
inverse mirror map, and ensure incomplete data fail explicitly.

```
bcov/
  __init__.py, cli.py       public imports and command line
  bcov_yy.py               exact ring, lazy recipes, BCOV solver
  quintic_frames.py        frame maps and tagged amplitudes
  datasets.py              typed, streaming bundled-data access
  patterns.py              exact jets, recurrences and arithmetic reports
  verify_physical_qseries.py  periods and GW/GV verification
  data/                    raw and derived research data, packaged in wheels
  *.py                     historical conversion/audit/verification scripts
notebooks/                 focused interactive demonstrations
tests/                     package-level regressions
Quintic_BCOV_Tutorial.py   extended guided laboratory
```

Historical scripts remain available for reproducing earlier outputs; direct
execution from `bcov/` is preserved. Their older helpers are not the typed
public analysis API. The existing call graphs in `docs/callgraph/` describe the
pre-packaging core and do not include the new analysis modules.

## Data provenance and rebuilding

The bundled free energies `ff_2.txt`–`ff_15.txt` and `quinticGVDTPT.m` are the
existing supplied research inputs. The source data have not been changed.
`GVTab.csv.gz` is the reference table for physical checks. Extended tables
contain unknown entries and 29 conflicting overlaps; see `gv_manifest.json`
and `GVTab_vs_GVTabExt_disagreements.csv`. These conflicts are not resolved by
this package. Source-amplitude versus GV agreement is an independent algebraic
check, not independent provenance or an all-genus proof.

To rebuild derived data from the existing raw inputs:

```bash
python -m bcov.bcov_sparse_converter --input-dir bcov/data --output /tmp/source_z.jsonl.gz --manifest /tmp/source_manifest.json
python -m bcov.bcov_reconcile_frames /tmp/source_z.jsonl.gz --output /tmp/source_hkq.jsonl.gz
python -m bcov.quintic_invariant_extract bcov/data/quinticGVDTPT.m --out /tmp/gv_tables
python -m bcov.verify_physical_qseries --atomic /tmp/source_hkq.jsonl.gz --gv /tmp/gv_tables/GVTab.csv.gz --genus 15 --degree 6 --output /tmp/rebuilt_residuals.json
```
