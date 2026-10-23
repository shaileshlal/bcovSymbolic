# Quintic BCOV exact algebra: convention freeze v0.2.0

**19 October 2026**. Canonical recursion: **original Yamaguchi–Yau (2004) fifth-root `s` frame**. The 2006 Huang–Klemm–Quackenbush (HKQ) frame and Gu/Klemm small-`z` frame are **explicit changes of variables**, not hidden assumptions.

This is a standalone, exact-arithmetic update of `bcov_yy.py` and `bcov_yy_conventions_and_derivations.tex`. It resolves the earlier genus-2 polynomial and GW/GV checks without altering the BCOV anomaly equations. The physical free energies are supplied through genus 15; numerical enumerative tests of those free energies are **limited to genus 15**.

## Run all tests

```bash
python run_all_checks.py
```

Requires Python 3.10+; SymPy is used only for `test_gu_gauge.py`, and the remainder of the tests and reconstruction use the standard library and `fractions.Fraction`. The suite also needs the `data/` directory included here. All assertions fail loudly if a coefficient or frame convention drifts.

Tests executed on the included files:

- Gu/Klemm differential propagator identities: **2 exact symbolic zero identities**.
- YY ↔ HKQ ↔ small-`z` round trips and canonical genus-2 published polynomial: **pass**.
- Retained YY BCOV recursion (`P_2` equality, anomaly source, and `P_3` formal recursion): **pass**.
- Bernoulli constant and linear pure-`X` coefficients from source `ff_2..ff_15`: **28/28 pass**.
- Full pure-`X` denominator-lattice invariance under HKQ ↔ `z`, and non-5 denominators under HKQ ↔ YY: **genera 2–8 pass**; the general basis-invariance statement has a direct algebraic proof.
- Physical A-model limit vs **`GVTab` from the uploaded Mathematica file**, exact `q^1,...,q^6` at genus 2,…,15: **84/84 pass**.

## Conventions ledger

| Item | Original YY | HKQ | Gu/Klemm source |
|---|---|---|---|
| Local coordinate | `s` | `psi=s**5` | `z=1/(3125*psi)` |
| Large-radius point | `s=∞` | `psi=∞` | `z=0` |
| Euler operator | `s*d/ds` | `psi*d/dpsi=theta_s/5` | `z*d/dz=-theta_psi` |
| Discriminant generator | `X=1/(1-s**5)` | `X=1/(1-psi)` | `Xz=1/(1-3125*z)` |
| Metric generators | `uY=5*uH`, `v1Y=5*v1H`, `v2Y=25*v2H`, `v3Y=125*v3H` | `uH=-uz`, `v1H=-v1z`, `v2H=v2z`, `v3H=-v3z` | `XH=1-Xz` |
| Normalized amplitude | `PY=5**(3*g-3)*PH(v1Y/5,v2Y/25,v3Y/125,X)` | `PH=(-1)**(g-1)*Pz(-v1,v2,-v3,1-X)` | `Pz=(5*Xz)**(g-1)*F_g^B` |

Here `C_t=t**3*C_ttt`, and the tensor's *cubic Jacobian* must be applied before comparing Yukawa couplings. The supplied Gu/Klemm holomorphic gauge functions are recorded in `quintic_frames.py` and in the TeX.

The 2004 YY finite-degree polynomial recursion remains in `bcov_yy.py`. Module `quintic_frames.py` provides maps `yy_to_hkq`, `hkq_to_yy`, `hkq_to_z`, `z_to_hkq`, `yy_to_z`, `z_to_yy` and Euler derivations in the other frames. Use `FramedPolynomial(poly,genus,"yy").to("hkq")` to keep frame metadata attached and prevent mixed-frame arithmetic. All transformations use immutable, sparse rational polynomial objects.

## Reproduction and datasets

- `bcov_yy_conventions_and_derivations_v2.tex` / `.pdf`: revised original technical note, with the new convention ledger **inserted** before the existing derivations.
- `bcov_convention_reconciliation_v2.tex`: stand-alone TeX section source used by the main manuscript.
- `bcov_yy.py`: canonical sparse YY ring and BCOV recursion, deliberately unchanged from its previously audited version.
- `quintic_frames.py`: frame-normalization and Euler-operator conversions, plus the immutable `FramedPolynomial` wrapper that refuses to add expressions from different frames or genera.
- `bcov_sparse_converter.py` + `bcov_reconcile_frames.py`: lossless rational conversion from source free energies to factored HKQ atomic polynomials.
- `audit_coefficients.py`: exact low-`X` coefficients (including Bernoulli convolution) of **all** supplied genus files.
- `audit_denominators.py`: tests for denominator-lattice invariance of entire pure-`X` polynomials.
- `verify_physical_qseries.py`: uses `GVTab.csv.gz` from **`quinticGVDTPT.m`**, and evaluates `F_g^A(q)` directly from the atomic factorizations, without full polynomial expansion.
- `quintic_invariant_extract.py`: safe `GridBox` integer parser for the Mathematica `.m` file; no code is executed.
- `data/ff_*.txt`, `data/quinticGVDTPT.m`, `data/GVTab.csv.gz`, `data/bcov_hkq_atomic_genus_2_15.jsonl.gz`: exact input and factored intermediate data for reproducing the results.
- `low_x_reaudit_genus_2_15.csv`, `denominator_invariants.csv`, `qseries_gv_residuals.json`, `regression.log`: independently rerunnable results.
- `data/GVTab_vs_GVTabExt_disagreements.csv`: preserves the **29 conflicting high-genus entries**; our free-energy files cannot adjudicate them. Unknown entries are not zeros.

To rebuild the intermediate HKQ data from the raw `ff_g` expressions:

```bash
python bcov_sparse_converter.py --input-dir data --output /tmp/source_z.jsonl.gz --manifest /tmp/source_manifest.json
python bcov_reconcile_frames.py /tmp/source_z.jsonl.gz --output /tmp/source_hkq.jsonl.gz
python verify_physical_qseries.py --atomic /tmp/source_hkq.jsonl.gz --gv data/GVTab.csv.gz --genus 15 --degree 6
```

The exact physical **genus-g GW coefficients** are not the integer GV numbers. For example, at genus two:

`[q] F_2^A = n^0_1/240 = 2875/240 = 575/48`.

For degree four, both primitive and multiple-cover terms enter:

`[q^4]F_2^A = n^0_4/240+n^2_4 + 2*n^0_2/240 + 4*n^0_1/240 = 1010821250`.

The physical genus-2 constant is `-5/144`, whereas the pure-`X` constant `a_{2,0}` in the *original YY polynomial* is `3125/144`. These must not be compared as the same object.

## Bernoulli and denominator statements: what is and is not invariant

**Exact physical input checked from g=2..15:**

`a[g,0] = 100 * 625**(g-1) * |B[2*g] B[2*g-2]| / (2*g*(2*g-2)*(2*g-2)!)`.

The previously derived exact `a[g,1]` recurrence (with Bernoulli and degree-one curve contributions) passes all the same tests. The normalized ratios `a[g,n]/a[g,0]` agree identically between original YY and HKQ because their pure-`X` polynomials differ by a uniform genus-dependent power of 5. The small-`z` coordinate `X -> 1-X` **mixes the coefficients**; do not compare raw `a[g,n]` across that change. The full vector's least common denominator *is* invariant under this integral unimodular transformation. Between YY/HKQ, its non-5-primary part is invariant.

This does **not** establish the previously attempted smallest-integer uniqueness conjecture, characterize all prime factors, or fix the entire holomorphic ambiguity without geometric boundary information. The first two coefficients survive the audit, and the remaining low-degree numerical patterns are retained as empirical tests (with further fixed-degree Taylor results through genus 53 documented in the earlier Bernoulli notes), not silently promoted to theorems.

## Remaining issue: two incompatible high-genus GV tables

`GVTab` has g=0..60, d=1..44. `GVTabExt` is partial through g=112, d=31, with explicit unknown placeholders. They disagree in 29 cells from g=26 onwards. **Both integer values are preserved.** Their resolution requires an independent PT/DT-to-GV consistency check, high-genus amplitude reconstruction, or authoritative source correction—not a normalizing factor or a low-genus q-series check.
