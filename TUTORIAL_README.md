# Quintic BCOV — guided notebook

The executed Jupyter notebook `Quintic_BCOV_Tutorial.ipynb` accompanies the
exact-rational BCOV package, release v0.2.0 (canonical original YY/2004 frame).

## Start

1. Extract `Quintic_BCOV_Tutorial_and_Package.zip` (do not flatten the directory).
2. `cd BCOV_conventions_v2`
3. In a Python 3.10+ environment with Jupyter available, run:

   ```bash
   python -m pip install jupyterlab ipykernel
   jupyter lab Quintic_BCOV_Tutorial.ipynb
   ```

4. In the notebook: **Kernel → Restart Kernel and Run All Cells**.

The notebook's core algebra is standard-library only. The optional comprehensive
`run_all_checks.py` needs SymPy for one propagator-gauge regression:

```bash
python -m pip install sympy
python run_all_checks.py
```

## The notebook shows

- Exact YY ring arithmetic, the genus-one one-point seed, BCOV genus 2 and 3.
- YY/HKQ/small-z frame conversion with explicit normalization factors.
- Extraction of `[X^k] P_g(v)` (polynomial) versus its pure-`X` coefficient.
- Streaming fixed-depth extraction from the included factorized Klemm `ff_g` data.
- Independent reconstruction of `P_2` and `P_3`, the Bernoulli a0/a1 checks,
  denominator-lattice checks, and GW/GV series comparisons.
- Exact CSV export to `notebook_output/low_X_coefficients_YY_genus_2_15.csv`.
- Safe checkpoint export and reload.

`examples/` contains the output table and genus-three boundary-fixed checkpoint
from a successfully executed notebook, for comparison without running it first.

## Mathematics and provenance cautions

- Canonical recursion uses *original YY fifth-root* `s`, not HKQ's `psi=s^5`
  or Gu/Klemm's `z=1/(3125 psi)`.
- The atomized `ff_g` files assume the **specific Gu/Klemm propagator gauge**
  documented in `CONVENTIONS.json`. Changing that gauge requires rederivation.
- `Amplitude.ambiguity=None` is unresolved; an empty tuple is an explicit
  zero choice, not quintic boundary data.
- The pure-`X` coefficients are not GW invariants. The physical `q`-expansion
  uses periods, mirror map and multi-covers; constant maps are separate.
- Tests supplied in v0.2.0 reach genus 15 for free energies and degrees 1–6
  for q-expansions, not an all-genus verification.
- User-supplied `GVTab` and `GVTabExt` disagree at 29 high-genus entries.
  The archive preserves both rather than silently selecting one.

All notebook arithmetic is exact `fractions.Fraction`, with zero numerical fitting.
