BCOV pure-X coefficient scan
============================

Read REPORT.pdf for the derivations, exact identities, empirical denominator
pattern, unsuccessful recurrence searches, and limitations.

Scope
-----
Physical quintic YY polynomials at genera 2 through 10, with x=z/mu and mu
symbolic. The normalized ring coefficients are independent of mu. The scan
covers 144 pure-X coefficients and 16,384 nonzero mixed/pure monomials.

The following are derived identities in the documented normalization:
  * leading and next-to-leading pure-X coefficients;
  * the mixed coefficient [v1 X^(3g-4)] P_g;
  * the constant-map coefficient;
  * a conifold/binomial recursion for the highest 2g-2 pure-X coefficients.

The large-prime denominator bound for full coefficient arrays is EMPIRICAL.
It agrees with the degree-2g Todd/Hirzebruch exponents for primes >=7 through
genus 10. It is not established for arbitrary genus. The constant-map term
already saturates the bound in the tested range. No compact predictive rule
for all interior numerators was found by the bounded ratio searches.

Run
---
Python 3.10 or newer; standard library only. Run from this directory:

  python scan_coefficients.py physical/amplitudes.json --output results
  python verify_triangular.py

The second command can take several minutes. It predicts the high pure-X
coefficients from the stored anomaly body and local period frame, then compares
to the physical checkpoint. It does not reuse the target pure-X coefficients
to construct its predictions. The bundled engine files are the exact algebra
and boundary modules from the original implementation.

The physical checkpoint was extended from genus 8 to genus 10 with the original
run_physical.py runner, using --gv-degree 4. It is compatible with that runner.
The runner's original boundary checks passed, including the anomaly residual,
conifold gap, MUM data, finite-monodromy regularity, and the imposed GV conditions.
This bundle contains the scan and its dependencies, not the full physical runner.

Files
-----
  REPORT.pdf / REPORT.tex: mathematical report and editable source.
  scan_coefficients.py: arithmetic scan, exact formulas, bounded recurrence fits.
  verify_triangular.py: independent conifold/binomial coefficient check.
  engine/: exact modules used by the triangular check.
  physical/amplitudes.json: full physical checkpoint through genus 10.
  physical/boundary_g09.json and boundary_g10.json: new-genus diagnostics.
  physical/run.json, plan.json, gv.json: computation metadata and GV data.
  results/summary.json: all scan diagnostics and held-out recurrence failures.
  results/pure_x_coefficients.json: exact coefficients and denominator factors;
    numerator factorization is only partial (primes <=97 plus an unfactored
    cofactor, explicitly labelled as such).
  results/triangular_checks.json: all 90 high-coefficient predictions.
  results/exact_edge_predictions_through_53.json: leading, next-to-leading,
    constant, and mixed-v1 edge values through genus 53. These are NOT full
    amplitudes above genus 10.
  results/*.tex: report tables.
  PROVENANCE.json: checkpoint hash, conventions, and verification scope.

Report regeneration requires LaTeX (pdflatex); it is not needed for the scans:
  pdflatex REPORT.tex
  pdflatex REPORT.tex

Verification completed
----------------------
45 edge/arithmetic checks against the nine physical polynomials.
Bernoulli binomial convolution independently checked through genus index 53.
90 conifold/binomial coefficient predictions matched exactly.
Fresh genera 9 and 10 agree with the patterns found at genera 2 through 8.

All arithmetic is exact fractions; no floating-point sequence fitting is used.
