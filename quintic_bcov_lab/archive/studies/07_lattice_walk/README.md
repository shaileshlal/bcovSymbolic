# Quintic reduced-ambiguity GV lattice-walk experiment

Run with Python 3.10+, standard library (NumPy/SciPy are *not* required to replay the exact audit):

```bash
python experiment.py
python charge_integrality_and_ml.py
python tests.py
```

The input `Quintic_Intermediate_ML/source/reference_amplitudes.json` is the corrected genus-2..12 YY archive from the earlier package; exact BCOV and mirror-map Python files are bundled alongside it. The same earlier trained-model prediction CSV is included for the naive-ML lattice-rounding control.

Outputs: `results/exact_lattice_audit.json` (full rational matrix and path), `results/lattice_coordinates.csv` (one row per residual coefficient), `results/unit_charge_steps.json` (exact x(1-x) response to one GV unit), `results/integral_transport_degree26.json` (572 exact integer checks), and `results/ml_lattice_projection.json`.

**Crucial scope:** independent affine BCOV solve for genera 4–10 (with a reduced series precision checked against references), retrospective reference anchor for 11–12; all target lower-degree GV values vanish. This does *not* establish minimum height for nonzero physically admissible BPS charges. See `RESEARCH_REPORT.md`.
