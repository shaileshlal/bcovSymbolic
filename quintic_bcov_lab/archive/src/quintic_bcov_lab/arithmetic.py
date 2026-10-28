"""Reusable exact arithmetic for the empirical stages. No floating-point labels."""
from __future__ import annotations
from fractions import Fraction
from math import comb


def fraction(value) -> Fraction:
    if isinstance(value, float):
        raise TypeError('A float cannot be used as an exact arithmetic label')
    return Fraction(value)


def difference(seq, order=1):
    seq=list(map(fraction,seq))
    if order < 0: raise ValueError('order must be >= 0')
    for _ in range(order):
        seq=[b-a for a,b in zip(seq,seq[1:])]
    return seq


def newton_forecast(seq, order=2, steps=1):
    """Polynomial extrapolation in Newton/binomial basis, exact fractions."""
    seq=list(map(fraction,seq))
    if steps<1 or len(seq)<order+1 or order<0:
        raise ValueError('insufficient data or invalid extrapolation parameters')
    # Use the latest order+1 observations for a LOCAL, not global, polynomial.
    a=seq[-order-1:]
    d=[difference(a,i)[0] for i in range(order+1)]
    return [sum((d[i]*comb(order+delta,i) for i in range(order+1)),Fraction())
            for delta in range(1,steps+1)]


def valuation(n,prime):
    if prime<=1:raise ValueError('prime>1 required')
    n=abs(int(n))
    if n==0:return None
    count=0
    while n%prime==0:n//=prime;count+=1
    return count


def sign(n):return (n>0)-(n<0)


def period_four_phase(k,variant='plus'):
    """Integer-valued sin(pi*k/2) +/- cos(pi*k/2) without numerical trig."""
    c=(1,0,-1,0)[k%4]
    s=(0,1,0,-1)[k%4]
    return c+s if variant=='plus' else c-s


def second_difference_pulse(seq):
    """Third finite differences: departures from constant curvature."""
    return difference(seq,3)


def sparse_pulse_reconstruct(first_three, pulses):
    """Reconstruct a Newton pulse sequence from initial values and third differences.

    pulses[j] is delta^3 at the first index j (i.e. new value at j+3).
    """
    seq=list(map(fraction, first_three))
    if len(seq)!=3:raise ValueError('Exactly three seeds required')
    for p in pulses:
        seq.append(fraction(p)+3*seq[-1]-3*seq[-2]+seq[-3])
    return seq


def envelope(s,t,k):
    c=(10,6,4,4,6)
    if s not in range(5):raise ValueError('s must be residue modulo 5')
    return 5*k*k+(2*s-5)*k+c[s]-2*t


def front(s,t):
    if s not in range(5):raise ValueError('s must be 0..4')
    return t-1 if s==0 else t+3-s


def terminal_pulse_depth(s,r):
    return 2*r-1 if s==1 else 2*r+3
