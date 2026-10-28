# High-X holomorphic ambiguity / GV–RoPE pilot

This archive is self-contained. It includes exact corrected YY amplitudes for g=2..12, a minimal exact BCOV boundary engine, numeric models, held-out results, and a speculative g=13 forecast.

Requirements: Python 3.10+ and NumPy. The exact verifier requires only the standard library.

From this directory:

```bash
python high_degree_phase_study.py
python verify_phases_and_digits.py
python verify_boundary.py
python verify_boundary.py --all
```

See `High_X_GV_RoPE_Study.md` for formulas, exact-versus-predictive provenance, holdout and ablation tables, citations, and limitations.

*Important:* the original genus-4 printed ambiguity discrepancy is corrected in this reference. The high-X boundary is conditional on physically fixed lower-genus amplitudes. The g=13 forecast is numerical, not an exact coefficient reconstruction. Only `a_g,0` is used as an exact oracle in the full-numerical forecast.
