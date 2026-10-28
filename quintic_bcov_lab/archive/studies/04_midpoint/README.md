# Exact midpoint involution test, mirror quintic

Run `python midpoint_test.py --genera 4 5 6 7 8 9` with Python 3.10+.
No NumPy, SciPy, SymPy or network is needed. All physics checks use `Fraction`.

The program checks the factorization of the *homogeneous* ambiguity
kernel and runs exact-reference-oracle examples at genera 7 and 9.
**It does not claim the involution supplies additional physical
boundary conditions**. Detailed derivation: `RESEARCH_NOTE.md`.

`quintic_ai_predictor/` is a limited copy of the corrected exact BCOV
reference and unchanged exact arithmetic/boundary modules from the
preceding experiment, included for standalone reproduction.
