"""Reusable exact arithmetic for the empirical stages. No floating-point labels."""
from __future__ import annotations
from collections.abc import Iterable
from typing import Literal, TypeAlias
from fractions import Fraction
from math import comb, isqrt

Exact: TypeAlias = int | str | Fraction


def fraction(value: Exact) -> Fraction:
    """Convert an exact label, rejecting floating-point inputs."""
    if isinstance(value, float):
        raise TypeError('A float cannot be used as an exact arithmetic label')
    return Fraction(value)


def difference(seq: Iterable[Exact], order: int = 1) -> list[Fraction]:
    """Return forward finite differences of the requested nonnegative order."""
    values = list(map(fraction, seq))
    if type(order) is not int or order < 0:
        raise ValueError('order must be a nonnegative integer')
    for _ in range(order):
        values = [b-a for a,b in zip(values, values[1:])]
    return values


def newton_forecast(seq: Iterable[Exact], order: int = 2, steps: int = 1) -> list[Fraction]:
    """Polynomial extrapolation in Newton/binomial basis, exact fractions."""
    seq=list(map(fraction,seq))
    if steps<1 or len(seq)<order+1 or order<0:
        raise ValueError('insufficient data or invalid extrapolation parameters')
    # Use the latest order+1 observations for a LOCAL, not global, polynomial.
    a=seq[-order-1:]
    d=[difference(a,i)[0] for i in range(order+1)]
    return [sum((d[i]*comb(order+delta,i) for i in range(order+1)),Fraction())
            for delta in range(1,steps+1)]


def valuation(n: int, prime: int) -> int | None:
    """Return the prime valuation of an integer; zero has no finite valuation."""
    if type(n) is not int or type(prime) is not int or prime < 2:
        raise ValueError("Require an integer and a prime >= 2")
    if any(prime % d == 0 for d in range(2, isqrt(prime)+1)):
        raise ValueError("The base must be prime")
    n=abs(n)
    if n==0:return None
    count=0
    while n%prime==0:n//=prime;count+=1
    return count


def sign(n: int | Fraction) -> int:
    """Return the exact sign, including zero."""
    return (n>0)-(n<0)


def period_four_phase(k: int, variant: Literal['plus', 'minus'] = 'plus') -> int:
    """Integer-valued sin(pi*k/2) +/- cos(pi*k/2) without numerical trig."""
    if variant not in ("plus", "minus"):
        raise ValueError("variant must be plus or minus")
    c=(1,0,-1,0)[k%4]
    s=(0,1,0,-1)[k%4]
    return c+s if variant=='plus' else c-s


def second_difference_pulse(seq: Iterable[Exact]) -> list[Fraction]:
    """Third finite differences: departures from constant curvature."""
    return difference(seq,3)


def sparse_pulse_reconstruct(first_three: Iterable[Exact], pulses: Iterable[Exact]) -> list[Fraction]:
    """Reconstruct a Newton pulse sequence from initial values and third differences.

    pulses[j] is delta^3 at the first index j (i.e. new value at j+3).
    """
    seq=list(map(fraction, first_three))
    if len(seq)!=3:raise ValueError('Exactly three seeds required')
    for p in pulses:
        seq.append(fraction(p)+3*seq[-1]-3*seq[-2]+seq[-3])
    return seq


def envelope(s: int, t: int, k: int) -> int:
    """Evaluate the empirically fitted quadratic shell envelope."""
    c=(10,6,4,4,6)
    if s not in range(5):raise ValueError('s must be residue modulo 5')
    return 5*k*k+(2*s-5)*k+c[s]-2*t


def front(s: int, t: int) -> int:
    """Return the moving front for one residue class modulo five."""
    if s not in range(5):raise ValueError('s must be 0..4')
    return t-1 if s==0 else t+3-s


def terminal_pulse_depth(s: int, r: int) -> int:
    """Return the empirical terminal pulse depth for residue s and shell r."""
    if s not in range(5):raise ValueError("s must be 0..4")
    return 2*r-1 if s==1 else 2*r+3
