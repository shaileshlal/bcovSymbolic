# Quintic BCOV research — unified sequential code package

**Version 1.3, 28 October 2026.** This repo consolidates **18 archived experimental
stages** (v0.1 through the nested-residual v1.2 work) into a reproducible
chronological research workflow. Its primary aims are (i) exact quintic BCOV
holomorphic-ambiguity reconstruction, (ii) prediction of **genuinely free**
GV-lattice integers, and (iii) empirically discovered PT/DT arithmetic and
nested finite-difference recurrences.

It contains the **actual archived Python source, exact rational calibration
labels, independent reference/validation code and marimo notebook scripts**;
not a rewrite of the reported results. The later-stage predictors are kept
separate from their reference validators exactly as in the source archives.

## Quick start

Python 3.10+; install scientific dependencies only for experiments that need
regression or modular-reference checks.

```bash
python -m pip install -e '.[science,test]'
python -m quintic_bcov_lab backend-check
python -m quintic_bcov_lab prepare
python -m quintic_bcov_lab list
python -m quintic_bcov_lab all --action test
python -m pytest -q tests
```

Test one study and then run the original prediction script:

```bash
python -m quintic_bcov_lab stage 05 --action test
python -m quintic_bcov_lab stage 05 --action run
python -m quintic_bcov_lab stage 18 --action test
python -m quintic_bcov_lab stage 18 --action run
```

The runner launches the original stage code from its original working directory
and records logs in `run_logs/`, including errors and elapsed time. All stages
are invoked sequentially in chronological order for `all`. The `--action run`
path reproduces **predictions**, while `--action test` executes each archived
stage's stated verification suite; the original stage code documents when its
validator consumes a withheld reference file. This distinction matters.

If not installing the package, substitute
`PYTHONPATH=src python -m quintic_bcov_lab ...`.

Interactive notebooks: `pip install marimo pandas matplotlib` and
`marimo edit notebooks/overview_marimo.py` or the `*marimo.py` notebook
within any `studies/<stage>/` directory. Notebook availability does not
impact the command-line exact tests.

## Backend and exact conventions

The pinned source of the exact solver is the archived quintic YY/BCOV engine
(`src/quintic_bcov_lab/vendor/engine`). All nine historical copies of
`bcov_yy.py` are byte-identical to this shared canonical copy; the SHA256 is
tracked. The **current external `bcov0.2` structural corpus** was not
available as exportable source during consolidation. Consequently this is
an integration of the accessible **exact YY backend**, not an unsupported
claim of executing an unknown v0.2 API. To opt in to a locally installed BCOV
backend with the same exact convention and `Poly`/`BCOV`/`YYRing` API:

```bash
BCOV_BACKEND_MODULE=your_package.bcov_yy python -m quintic_bcov_lab backend-check
```

This override applies to the **unified analysis modules**; archived stage
scripts continue to use their pinned exact engine so historical results remain
reproducible. An external backend with a different gauge is explicitly rejected.

Canonical gauge: `x=psi**-5=3125*z_PF`, `X=x/(x-1)`,
`P_g=P_g^an + a_g(X)`, and `P_g^an(0,0,0,X)=0`. The genus-four seed is
`P4_printed + (6746/5**12) X**5 (1-X)**2`.

The exact residual geometry is central:

```text
K=3g-3   m=ceil(K/5)   D=floor((2g-2)/5)

delta a_g(X) = X (1-X)^m sum_{r=0}^{D-1} c[g,r] X^r

delta f_g(x) = sum_{r=0}^{D-1} (-1)^(g-2-r) c[g,r]
                    x^(r+1) (1-x)^(D-r-1).
```

The unified `sections.py` verifies these transformations **over the
rationals** and generates exactly 22 intermediate targets at genera 2–12.
Their zero-pivot affine anchor is **reference-derived**, not a blind
prediction. The archived boundary modules separately solve the physical
conifold/orbifold/MUM systems where tested.

## The chronological research pipeline

| Stage | Archived work | Research target |
|---|---|---|
| 01 | v0.1 low-X ML | BCOV jets, rational digit encoding, inverse-genus cumulants, mixed YY features |
| 02 | High-X / GV–RoPE | Reverse cumulants, phase-feature ablations and boundary constraints |
| 03 | Two-stream | Mixed YY and GV/MUM predictor comparisons |
| 04 | Midpoint | Coordinate/section transformation and residual kernel |
| 05 | Intermediate dataset | Exact `x(1-x)` reduced coordinates and ML features |
| 06 | Intermediate ML | Magnitude, sign, local-source and shape baselines |
| 07 | GV lattice walk | Minimum-height displacement and exact triangular transport |
| 08 | GV displacement | Extremal sign periods and nearest-degree transport |
| 09 | Genus 45–53 | Single exceptional charge and exact BCOV sensitivity |
| 10 | DT/PT v0.4 | Degree-residue prediction, blind PT/GV conversion and optional DT solver |
| 11 | PT boundary v0.5 | Curvature/finite differences and DT fingerprints |
| 12 | Defect v0.6 | Sparse third differences and active-front classification |
| 13 | Recursive v0.7 | Newton pulse decomposition, mirrored curvature |
| 14 | Survival v0.8 | Degree 31–53 teacher-forced/autonomous holdout |
| 15 | Deep v0.9 | Shells 6–10, quadratic envelope and gated exact predictions |
| 16 | Pulses v1.0 | Affine inverse-envelope profiles and paired-residue transport |
| 17 | Depth v1.1 | Shells 11–15, anchor rank, recurrence-closure checks |
| 18 | Nested v1.2 | Terminal pulses and the residual of the residual |

Every directory in `studies/` contains the corresponding original README and
research report where available, *and* their supplied input/output snapshots.
The mapping to original ZIP archives and SHA256 hashes is in
`shared/PROVENANCE.json`, `shared/source_files.json`, and
`shared/stage_catalog.json`.

## Math and experiment modules

- `quintic_bcov_lab.backend`: type-checked exact BCOV engine adapter,
  independent genus-2 anomaly residual check.
- `quintic_bcov_lab.sections`: exact homogeneous kernel, rational
  `x(1-x)` reorganization, triangular GV response decoding.
- `quintic_bcov_lab.arithmetic`: lossless Fractions, prime valuations,
  finite differences, period-four integer phases, Newton pulses,
  quadratic envelope and moving-front functions.
- `quintic_bcov_lab.datasets`: reference-gauge checking and provenance-tagged
  intermediate labels.
- `quintic_bcov_lab.runner`: original stage run/test orchestration,
  compatibility aliases, SHA audits, and optional large dataset installation.

The main regression suite uses `pytest`; the original stage-specific tests
are separate and remain available with `all --action test`.

## Optional complete exact genus-1000 training data

The long low-X exact BCOV quotient datasets are **not required for the quick
suite**. To avoid a 20+ MB main download and duplicate enormous fraction
arrays, they are in `Quintic_BCOV_Optional_Exact_LowX_Data.zip` (separate).
The v0.1 stage already includes its 29 KB numerical training cache, exact
Bernoulli/Airy source routines, and saved benchmarks. The optional archive
contains `direct_J4_exact.json.gz` and `all_constant_and_leading.json.gz`.
Install it once extracted:

```bash
python -m quintic_bcov_lab install-large-data \
  --source ../Quintic_BCOV_Optional_Exact_LowX_Data.zip
cd studies/01_v01_low_x
python -m hybrid.model --data results/gmp/direct_J4_exact.json.gz \
  --output hybrid/results
```

The numeric cache is **not** an exact arithmetic certificate; exact results
are in the optional files or may be recomputed with the BCOV quotient code.

## Correctness and research hygiene

1. **Exact versus learned:** BCOV recursions, reference rational polynomials,
   conifold/orbifold/MUM projections, and GV-lattice transport are exact
   computations. Learned signed/logarithmic ratios or recurrent pulse profiles
   are separate empirical hypotheses.
2. **Modularity-free predictors:** PT stages 10–18 are trained on arithmetic
   PT observations only. Their *reference* PT integers were generated under a
   DT/modularity and wall-crossing framework; independence of the empirical
   prediction code is not an unconditional proof of those inputs.
3. **No-oracle tests:** many late archived tests copy the predictor and its
   calibration data to a temporary directory without the reference/validator.
   They do **not** imply that the model family was selected before researchers
   inspected the larger dataset.
4. **Abstentions count:** a model with 59/59 exact and six abstentions is 59
   exact *predictions*, **not** 65 of 65 exact predictions.
5. **Finite generation remains open:** exact affine pulse profiles do not
   show that the mature-scale sequences have a finite recursive description.
6. **Historical result snapshots:** no result recorded under `studies/*/results`
   should be described as newly generated unless its script was rerun and the
   corresponding `run_logs/` entry shows success.

See `BASELINES.md` and `REPRODUCIBILITY.md` for the quantitative results and
how to audit each stage.
