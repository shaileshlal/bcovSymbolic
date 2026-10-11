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

    return dataclass, defaultdict, sp


@app.cell
def _(dataclass, defaultdict, format_term, sp):
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


    ONE = Monomial()


    # ----------------------------------------------------------------------
    # Differential-ring element
    #
    # f = sum_m c_m m
    #
    # where c_m in K(z).
    # ----------------------------------------------------------------------

    class BCOV:
        """
        Sparse element of the BCOV propagator algebra

            K(z)[Szz, Sz, S, Kz].

        Coefficients are exact rational functions of z.
        """

        def __init__(self, terms=None):
            self.terms = defaultdict(lambda: sp.Integer(0))

            if terms:
                for m, c in terms.items():
                    c = rat(c)
                    if c != 0:
                        self.terms[m] += c

            self._cleanup()

        def _cleanup(self):
            zero = []
            for m, c in self.terms.items():
                c = rat(c)
                if c == 0:
                    zero.append(m)
                else:
                    self.terms[m] = c

            for m in zero:
                del self.terms[m]

        # --------------------------------------------------------------
        # Constructors
        # --------------------------------------------------------------

        @classmethod
        def scalar(cls, c):
            return cls({ONE: c})

        @classmethod
        def generator(cls, name):
            names = {
                "Szz": Monomial(Szz=1),
                "Sz":  Monomial(Sz=1),
                "S":   Monomial(S=1),
                "Kz":  Monomial(Kz=1),
            }

            return cls({names[name]: sp.Integer(1)})

        # --------------------------------------------------------------
        # Basic arithmetic
        # --------------------------------------------------------------

        def __add__(self, other):
            other = bcov_ize(other)

            out = defaultdict(lambda: sp.Integer(0))
            out.update(self.terms)

            for m, c in other.terms.items():
                out[m] += c

            return BCOV(out)

        __radd__ = __add__

        def __neg__(self):
            return BCOV({m: -c for m, c in self.terms.items()})

        def __sub__(self, other):
            return self + (-bcov_ize(other))

        def __rsub__(self, other):
            return bcov_ize(other) - self

        def __mul__(self, other):
            other = bcov_ize(other)

            out = defaultdict(lambda: sp.Integer(0))

            for m1, c1 in self.terms.items():
                for m2, c2 in other.terms.items():
                    out[m1 * m2] += c1 * c2

            return BCOV(out)

        __rmul__ = __mul__

        def __pow__(self, n):
            if n < 0:
                raise ValueError("Negative powers are not implemented.")

            out = BCOV.scalar(1)

            for _ in range(n):
                out *= self

            return out

        # --------------------------------------------------------------
        # Coefficient derivative
        # --------------------------------------------------------------

        def Dz(self):
            """
            Total BCOV derivative D_z.

            This preliminary version uses the polynomial rules

                D_z(Szz), D_z(Sz), D_z(S), D_z(Kz)

            supplied by `bcov_derivatives`.

            The coefficient part is always included.
            """

            out = defaultdict(lambda: sp.Integer(0))

            for m, c in self.terms.items():

                # derivative of coefficient
                dc = Dz_coeff(c)

                if dc != 0:
                    out[m] += dc

                # derivative of the monomial
                out_bc = derivative_monomial(m)

                for m2, c2 in out_bc.terms.items():
                    out[m2] += c * c2

            return BCOV(out)

        # --------------------------------------------------------------
        # Utilities
        # --------------------------------------------------------------

        def simplify(self):
            return BCOV(self.terms)

        def is_zero(self):
            return len(self.terms) == 0

        def max_degree(self):
            if not self.terms:
                return 0
            return max(m.degree() for m in self.terms)

        def __repr__(self):
            if not self.terms:
                return "0"

            pieces = []

            for m, c in sorted(self.terms.items()):
                pieces.append(format_term(c, m))

            return " + ".join(pieces)


    def bcov_ize(x):
        if isinstance(x, BCOV):
            return x
        return BCOV.scalar(x)


    # ----------------------------------------------------------------------
    # Generators
    # ----------------------------------------------------------------------

    Szz = BCOV.generator("Szz")
    Sz  = BCOV.generator("Sz")
    S   = BCOV.generator("S")
    Kz  = BCOV.generator("Kz")


    # ----------------------------------------------------------------------
    # Differential relations
    #
    # THESE are deliberately isolated.
    #
    # We should eventually put the precise BCOV/Yamaguchi-Yau relations
    # here rather than hard-wire them throughout the algebra engine.
    # ----------------------------------------------------------------------

    bcov_derivatives = {
        "Szz": None,
        "Sz":  None,
        "S":   None,
        "Kz":  None,
    }


    def set_bcov_derivative(name, expression):
        """
        Set D_z(generator) = expression.
        """
        if name not in bcov_derivatives:
            raise KeyError(name)

        bcov_derivatives[name] = bcov_ize(expression)


    def derivative_generator(name):
        expr = bcov_derivatives[name]

        if expr is None:
            raise ValueError(
                f"D_z({name}) has not yet been specified."
            )

        return expr


    def derivative_monomial(m):
        """
        Leibniz rule for a propagator monomial.
        """

        out = BCOV.scalar(0)

        powers = {
            "Szz": m.Szz,
            "Sz":  m.Sz,
            "S":   m.S,
            "Kz":  m.Kz,
        }

        for name, power in powers.items():

            if power == 0:
                continue

            reduced = {
                "Szz": m.Szz,
                "Sz":  m.Sz,
                "S":   m.S,
                "Kz":  m.Kz,
            }

            reduced[name] -= 1

            reduced_m = Monomial(**reduced)

            derivative = derivative_generator(name)

            out += BCOV({reduced_m: power}) * derivative

        return out


    # ----------------------------------------------------------------------
    # Differential operators
    #
    # We explicitly represent
    #
    #       D_z f = f D_z + D_z(f)
    #
    # rather than pretending the operator commutes with f.
    # ----------------------------------------------------------------------

    @dataclass(frozen=True)
    class DifferentialOperator:
        """
        A sparse differential operator

            sum_k a_k D_z^k.
        """

        terms: dict[int, BCOV]

        def __post_init__(self):
            object.__setattr__(
                self,
                "terms",
                {
                    k: v
                    for k, v in self.terms.items()
                    if not v.is_zero()
                }
            )

        @classmethod
        def Dz(cls):
            return cls({1: BCOV.scalar(1)})

        @classmethod
        def multiplication(cls, f):
            return cls({0: bcov_ize(f)})

        def __add__(self, other):
            out = defaultdict(lambda: BCOV.scalar(0))

            for k, a in self.terms.items():
                out[k] = out[k] + a

            for k, a in other.terms.items():
                out[k] = out[k] + a

            return DifferentialOperator(dict(out))

        def __repr__(self):
            if not self.terms:
                return "0"

            pieces = []

            for k, a in sorted(self.terms.items()):
                if k == 0:
                    pieces.append(f"({a})")
                elif k == 1:
                    pieces.append(f"({a}) Dz")
                else:
                    pieces.append(f"({a}) Dz^{k}")

            return " + ".join(pieces)


    Dz = DifferentialOperator.Dz()


    # ----------------------------------------------------------------------
    # Operator action
    # ----------------------------------------------------------------------

    def apply_operator(op, f):
        """
        Apply a differential operator to f.
        """

        out = BCOV.scalar(0)

        for k, coeff in op.terms.items():

            g = f

            for _ in range(k):
                g = g.Dz()

            out += coeff * g

        return out

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
