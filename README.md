# bcovSymbolic
My computations of the BCOV equations, Holomorphic Anomaly and Yamaguchi-Yau polynomials
for the quintic Calabi-Yau threefold.

Each folder below is self-contained: run its scripts from inside that folder.
Some small modules (`yy_p2_p3_exact.py`, the `bcov_yy` engine files) are
deliberately duplicated so that every folder runs on its own.

## Contents

| Folder | What is there |
|---|---|
| [`bcov_yy/`](bcov_yy/) | Main exact BCOV engine in Yamaguchi-Yau generators (standard library only). Sparse exact polynomial algebra, the anomaly recursion, the physical boundary solver (conifold gap, orbifold regularity, constant maps, GV/Castelnuovo data), and the checkpointed driver `run_physical.py`. Validated physical amplitudes through genus 8 are in `validated_physical/`. See [README](bcov_yy/README.md) and [BOUNDARIES](bcov_yy/BOUNDARIES.md). |
| [`coefficient_scan/`](coefficient_scan/) | Arithmetic scan of the pure-X coefficients of the physical YY polynomials at genera 2-10: derived edge identities, empirical denominator bounds, and the conifold/binomial triangular check. It bundles a copy of the engine in `engine/` and the genus-10 checkpoint in `physical/`. See [README.txt](coefficient_scan/README.txt) and `REPORT.pdf`. |
| [`yy p2 p3 exact/`](yy%20p2%20p3%20exact/) | `yy_p2_p3_exact.py`: the exact genus-2 and genus-3 YY polynomials, basis transformations and enumerative sanity checks (SymPy). It is the base module the other `yy *` folders use. See [report](yy%20p2%20p3%20exact/yy_p2_p3_report.md). |
| [`yy bootstrap/`](yy%20bootstrap/) | Pullback of published modular-bootstrap / near-Castelnuovo GV and PT data to the YY ring at genera 2-3, with Fourier/GV-dual ambiguity bases. See [report](yy%20bootstrap/yy_modular_bootstrap_report.md). |
| [`yy connections/`](yy%20connections/) | Mirror-quintic connection action on the YY generators: exact monodromy matrices, the numerical LCS-to-conifold Frobenius connection, a ring-closure check, and the induced P2/P3 cocycles. See [report](yy%20connections/yy_connection_report.md). |
| [`yy crossing gap equations/`](yy%20crossing%20gap%20equations/) | Towards Polyakov-type blocks for BCOV: finite-dimensional gap/LCS crossing bases at genera 2-3. See [report](yy%20crossing%20gap%20equations/yy_crossing_gap_report.md). |
| [`metaplectic lift/`](metaplectic%20lift/) | The metaplectic lift of the conifold monodromy (a rank-one symplectic shear) on the BCOV wave function: exact connected Gaussian graph corrections through genus 3. Includes a script, its test output, and the note `bcov_metaplectic_lift_note.pdf`. |
| [`toySymbolic/`](toySymbolic/) | [marimo](https://marimo.io) notebooks prototyping a symbolic differential-algebra engine over Q(z): `toyCSA.py` is a toy model and `bcovCSA.py` the BCOV version. Open them with `marimo edit`. |

## Requirements

Python 3.10+. `bcov_yy/` and `coefficient_scan/` use only the standard library.
The `yy *` folders and `metaplectic lift/` need SymPy, and `yy connections/` also needs mpmath.
