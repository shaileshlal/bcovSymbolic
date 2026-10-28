# Sequential research v1.3 integration

The consolidated archive supplied on 28 October 2026 is integrated into
`bcov-symbolic` 0.4.0 on the existing `toolkit` branch. It includes 18 historical
stages and 528 files. The original README, baseline tables, reproduction notes,
license notice, and historical validation logs remain in
[`quintic_bcov_lab/archive`](../../quintic_bcov_lab/archive/README.md).
The separately supplied Markdown files were byte-identical to those in the ZIP.

## What changed

- `quintic_bcov_lab.backend` now defaults to the live `bcov.bcov_yy` backend,
  checks its convention, and validates the exact genus-two anomaly equation.
  Explicit `BCOV_BACKEND_MODULE` overrides remain supported.
- Shared section, arithmetic, data-access, and runner APIs have type annotations
  and docstrings. Exact section validation compares every reduced polynomial
  coefficient. Reference reads verify the recorded SHA256.
- The original stage code, seeds, calibration tables, validators, result snapshots,
  and pinned engines are preserved byte-for-byte. A new import manifest records
  the supplied ZIP's digest and every extracted file's digest.
- The runner copies study inputs into `.bcov-research/` (or `--workspace`) and
  restores shared references there. Stage commands run from their historical
  working directories, in chronological order, with full logs and JSON status.
  The archive is packaged in the wheel and works outside the checkout.
- Existing and new APIs share one installation. The `bcov` CLI remains intact;
  `quintic-bcov` exposes the research workflow. A new interactive marimo
  overview presents the stage catalog, exact section labels, a calibration-only
  terminal-pulse example, and clearly labeled historical snapshots.

The integrated runner intentionally does **not** replace the original stage
engines with the current engine. Keeping those sources fixed allows exact
reproduction and preserves their original import paths and model assumptions.

## Fresh validation

Machine-readable results are in
[`integration_validation.json`](integration_validation.json). These results
were generated during integration, independently of the archived release's
`RUN_VALIDATION.json`.

- 21 package and integration tests pass.
- Static type checking passes for all 11 maintained public analysis/adapter modules;
  archived historical source is deliberately excluded from this check.
- All 18 original stage test commands pass sequentially in a fresh workspace.
- All 528 imported files match the import manifest. The historical source audit
  reports 480 matching files, ten shared references and two optional files,
  with no changed or missing bundled files.
- Live and pinned solver anomaly polynomials agree through formal genus 3.
- All eleven supplied ambiguity vectors at genera 2–12 agree with coefficients
  extracted independently from the repository's HKQ atomic data.
- The pre-existing full regression suite passes, including all 84 positive-degree
  GW/GV comparisons at genera 2–15 and degrees 1–6.
- The research overview executes headlessly. An installed wheel is checked from
  outside the repository for backend validation, complete archive integrity,
  and the original stage-18 test command.

The research tests ran with Python 3.14 in the existing `jax` environment,
which already provides NumPy, SciPy, scikit-learn, SymPy, pandas, pytest, and
marimo. Wheel construction used Python 3.11. Sage was not required.

## Research interpretation

The imported benchmarks evaluate different targets, fitting windows, supplied
anchors, and abstention policies. They are not a common leaderboard. Preserve
these distinctions when consuming results:

1. The 22 intermediate labels and their zero-pivot affine representatives are
   reference-derived, not blind predictions.
2. PT/DT reference data retain their modularity and wall-crossing assumptions.
   A predictor that does not import those formulas is still being scored
   against conditional reference data.
3. The original no-oracle tests are retained; retrospective model-family
   selection is not erased by a later holdout.
4. Missing predictions and abstentions are not exact successes. The nested
   secondary-pulse benchmark includes two explicit misses.
5. Mature-scale or PT-anchor inputs make the corresponding diagnostics
   conditional; finite recurrence closure and full all-genus ambiguity
   reconstruction remain open.

The optional exact genus-1000 archive was not supplied. Its numerical cache
is present, but the missing exact files are neither fabricated nor inferred.
Use `quintic-bcov install-large-data --source ...` if those files become available.
