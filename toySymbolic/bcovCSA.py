import marimo

__generated_with = "0.24.2"
app = marimo.App()


@app.cell
def _():
    from __future__ import annotations

    from dataclasses import dataclass
    from fractions import Fraction
    from collections import defaultdict

    import sympy as sp

    return dataclass, sp


@app.cell
def _(dataclass, sp):
    # ----------------------------------------------------------------------
    # Base coefficient field
    #
    # K(z) = rational functions in z.
    # We use SymPy here only for the coefficient field, NOT for the
    # differential-ring expression tree.
    # ----------------------------------------------------------------------

    z = sp.Symbol("z")


    def rat(x):
        """Convert x to an exact SymPy rational function."""
        return sp.cancel(sp.sympify(x))


    def Dz_coeff(f):
        """Derivative in the coefficient field K(z)."""
        return sp.cancel(sp.diff(f, z))


    # ----------------------------------------------------------------------
    # Propagator monomials
    #
    # (a,b,c,d) represents
    #
    #     (Szz)^a (Sz)^b S^c (Kz)^d
    #
    # ----------------------------------------------------------------------

    @dataclass(frozen=True, order=True)
    class Monomial:
        Szz: int = 0
        Sz: int = 0
        S: int = 0
        Kz: int = 0

        def __mul__(self, other):
            return Monomial(
                self.Szz + other.Szz,
                self.Sz + other.Sz,
                self.S + other.S,
                self.Kz + other.Kz,
            )

        def degree(self):
            return self.Szz + self.Sz + self.S + self.Kz

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
