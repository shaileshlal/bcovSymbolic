"""Exact, convention-aware arithmetic patterns supported by the source data.

The Bernoulli recurrences below are verified through genus 15, not certified
all-genus formulas for the entire ambiguity. Inverting X mixes coefficients:
small-z jets must be obtained from the full HKQ pure-X polynomial.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction as Q
from math import comb, factorial, isqrt, lcm

from .datasets import AtomicAmplitude
from .quintic_frames import Frame


def bernoulli_abs(max_genus: int) -> dict[int, Q]:
    """Return |B_(2g)| for 1 <= g <= max_genus by rational convolution."""
    if type(max_genus) is not int or max_genus < 1:
        raise ValueError("max_genus must be a positive integer")
    b = {1: Q(1, 6)}
    for g in range(2, max_genus + 1):
        b[g] = sum((comb(2*g, 2*h)*b[h]*b[g-h] for h in range(1, g)), Q(0))/(2*g+1)
    return b


def low_x_predictions(max_genus: int) -> dict[int, tuple[Q, Q]]:
    """Return predicted YY (constant, linear) coefficients through max_genus.

    Implements equations (2026bernoulli) and (2026a1) in the reconciliation
    notes, including the separate genus-two linear seed.
    """
    if type(max_genus) is not int or max_genus < 2:
        raise ValueError("max_genus must be an integer >= 2")
    b = bernoulli_abs(max_genus)
    c = {g: Q(100*625**(g-1))*b[g]*b[g-1]/(2*g*(2*g-2)*factorial(2*g-2))
         for g in range(2, max_genus+1)}
    result = {2: (c[2], Q(-28795, 144))}
    for g in range(3, max_genus+1):
        convolution = sum(((2*h-2)*(2*g-2*h-2)*c[h]*c[g-h] for h in range(2, g-1)), Q(0))
        inhomogeneous = Q(23, 25)*(-625)**(g-1)*b[g]/(2*g*factorial(2*g-2))
        linear = -Q(577, 625)*(g-1)*c[g] + (2*g-4)*(62-Q(24*g, 5))*c[g-1] - Q(12, 5)*convolution - inhomogeneous
        result[g] = c[g], linear
    return result


def invert_x(coefficients: Iterable[Q], genus: int) -> tuple[Q, ...]:
    """Transform a complete HKQ vector into (-1)^(g-1) a(1-X)."""
    if type(genus) is not int or genus < 2:
        raise ValueError("genus must be an integer >= 2")
    f = tuple(coefficients)
    return tuple((-1)**(genus-1+s)*sum((f[k]*comb(k, s) for k in range(s, len(f))), Q(0))
                 for s in range(len(f)))


def pure_x(amplitude: AtomicAmplitude, depth: int | None = None, frame: Frame = "yy") -> tuple[Q, ...]:
    """Extract the pure-X jet through depth inclusively, without expanding Poly.

    YY/HKQ extraction uses valuation n+ell to prune terms. For small_z the
    complete vector is transformed first: truncation and X inversion do not
    commute. A depth beyond 3g-3 is padded with zeros.
    """
    degree = 3*amplitude.genus-3
    depth = degree if depth is None else depth
    if type(depth) is not int or depth < 0:
        raise ValueError("depth must be a nonnegative integer")
    if frame not in ("yy", "hkq", "small_z"):
        raise ValueError(f"Unknown frame: {frame}")
    limit = degree if frame == "small_z" else min(depth, degree)
    result = [Q(0)]*(limit+1)
    for t in amplitude.terms:
        start = t.n+t.ell
        if start > limit:
            continue
        factor = t.coefficient*Q(-3, 5)**t.j*Q(2, 25)**t.k*Q(-4, 125)**t.ell
        for s in range(min(t.r, limit-start)+1):
            result[start+s] += factor*(-1)**s*comb(t.r, s)
    if frame == "yy":
        result = [x*5**degree for x in result]
    elif frame == "small_z":
        result = list(invert_x(result, amplitude.genus))
    return tuple((result + [Q(0)]*max(0, depth+1-len(result)))[:depth+1])


def denominator_lcm(coefficients: Iterable[Q]) -> int:
    """Return the least positive common denominator, including 1 for zeros."""
    return lcm(*(x.denominator for x in coefficients))


def prime_free_part(value: int, prime: int = 5) -> int:
    """Remove every factor of a supplied prime (5 for the YY/HKQ rescaling)."""
    if value < 1 or prime < 2 or any(prime % d == 0 for d in range(2, isqrt(prime)+1)):
        raise ValueError("Require a positive value and a prime >= 2")
    while value % prime == 0:
        value //= prime
    return value


@dataclass(frozen=True)
class PatternReport:
    """Exact observations and predictions for one source genus."""

    genus: int
    yy_jet: tuple[Q, ...]
    predicted_low_x: tuple[Q, Q]
    lcd_hkq: int
    lcd_small_z: int
    lcd_yy: int

    @property
    def low_x_matches(self) -> bool:
        """Whether both source coefficients equal their Bernoulli predictions."""
        return self.yy_jet[:2] == self.predicted_low_x

    @property
    def denominator_matches(self) -> bool:
        """Test the full denominator lattice in all three frames."""
        return (self.lcd_hkq == self.lcd_small_z and
                prime_free_part(self.lcd_hkq) == prime_free_part(self.lcd_yy))

    @property
    def ratios(self) -> tuple[Q, ...]:
        """Frame-independent YY/HKQ ratios to the constant coefficient."""
        return tuple(x/self.yy_jet[0] for x in self.yy_jet)


def audit_patterns(amplitude: AtomicAmplitude, depth: int = 4) -> PatternReport:
    """Check low-X formulas and full denominator lattices against source data."""
    if type(depth) is not int or depth < 1:
        raise ValueError("Auditing constant and linear coefficients requires depth >= 1")
    h = pure_x(amplitude, frame="hkq")
    yy = tuple(x*5**(3*amplitude.genus-3) for x in h)
    return PatternReport(
        amplitude.genus, pure_x(amplitude, depth), low_x_predictions(amplitude.genus)[amplitude.genus],
        denominator_lcm(h), denominator_lcm(invert_x(h, amplitude.genus)), denominator_lcm(yy),
    )
