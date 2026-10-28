"""Exact Laurent jets. Every operation propagates an absolute O(t**N) bound.

A coefficient beyond the certified bound raises, rather than silently becoming
zero. Exact finite polynomials are supplied with a finite working bound too.
"""
from fractions import Fraction as Q
from math import lcm


class Jet:
    def __init__(self, coefficients=None, precision=32):
        self.precision = precision
        self.c = {k: Q(v) for k, v in (coefficients or {}).items()
                  if v and k < precision}

    @property
    def valuation(self):
        return min(self.c, default=self.precision)

    def coeff(self, k):
        if k >= self.precision:
            raise ValueError(f'Uncertified coefficient {k}: precision {self.precision}')
        return self.c.get(k, Q(0))

    def truncate(self, precision):
        if precision > self.precision:
            raise ValueError('Cannot increase certified precision')
        return Jet(self.c, precision)

    def _coerce(self, value):
        return value if isinstance(value, Jet) else Jet({0: Q(value)}, self.precision + abs(self.valuation) + 8)

    def __add__(self, other):
        other = self._coerce(other)
        n = min(self.precision, other.precision)
        c = dict(self.c)
        for k, v in other.c.items():
            c[k] = c.get(k, 0) + v
        return Jet(c, n)

    __radd__ = __add__

    def __neg__(self):
        return Jet({k: -v for k, v in self.c.items()}, self.precision)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return Jet({k: v*Q(other) for k, v in self.c.items()}, self.precision)
        n = min(self.precision + other.valuation, other.precision + self.valuation)
        # Reduce rationals once per output coefficient, not per pair product.
        da=lcm(*(a.denominator for a in self.c.values()))
        db=lcm(*(b.denominator for b in other.c.values()))
        aa={i:a.numerator*(da//a.denominator) for i,a in self.c.items()}
        bb={j:b.numerator*(db//b.denominator) for j,b in other.c.items()}
        c = {}
        for i, a in aa.items():
            for j, b in bb.items():
                if i+j < n:
                    c[i+j] = c.get(i+j, 0) + a*b
        return Jet({k:Q(v,da*db) for k,v in c.items()}, n)

    __rmul__ = __mul__

    def inverse(self):
        if not self.c:
            raise ZeroDivisionError('Cannot invert a zero jet')
        v = self.valuation
        count = self.precision-v
        a0 = self.c[v]
        b = [1/a0]
        for k in range(1, count):
            b.append(-sum((self.c.get(v+i, 0)*b[k-i] for i in range(1,k+1)), Q(0))/a0)
        return Jet({i-v: x for i, x in enumerate(b)}, self.precision-2*v)

    def __truediv__(self, other):
        return self * (other.inverse() if isinstance(other, Jet) else 1/Q(other))

    def __pow__(self, n):
        if n < 0:
            return self.inverse() ** (-n)
        result = Jet({0: 1}, self.precision + abs(self.valuation*n) + 8)
        base = self
        while n:
            if n & 1:
                result = result*base
            n >>= 1
            if n:
                base = base*base
        return result

    def derivative(self):
        return Jet({k-1: k*v for k,v in self.c.items() if k}, self.precision-1)

    def euler(self):
        return Jet({k:k*v for k,v in self.c.items()}, self.precision)

    def exp(self):
        if self.valuation < 1:
            raise ValueError('exp requires zero constant and no negative powers')
        a = [Q(1)]
        for k in range(1, self.precision):
            a.append(sum((i*self.c.get(i,0)*a[k-i] for i in range(1,k+1)), Q(0))/k)
        return Jet(dict(enumerate(a)), self.precision)

    def compose(self, other):
        if self.valuation < 0 or other.valuation < 1:
            raise ValueError('Composition requires power series and positive inner valuation')
        result = Jet({}, min(self.precision*other.valuation, other.precision))
        # Horner avoids computing every power independently.
        for k in range(self.precision-1, -1, -1):
            result = result*other + self.c.get(k, 0)
        return result


def evaluate_poly(poly, values):
    """Nested sparse Horner evaluation, with certified Laurent precision."""
    n = min(v.precision for v in values)
    powers = {}
    def power(i,e):
        key=(i,e)
        if key not in powers: powers[key]=values[i]**e
        return powers[key]
    def horner(terms,i):
        if not terms: return Jet({},n)
        if i==5: return Jet({0:sum((c for m,c in terms),Q(0))},n)
        groups={}
        for m,c in terms: groups.setdefault(m[i],[]).append((m,c))
        ordered=sorted(groups,reverse=True)
        last=ordered[0]
        result=horner(groups[last],i+1)
        for exponent in ordered[1:]:
            result=result*power(i,last-exponent)+horner(groups[exponent],i+1)
            last=exponent
        if last: result=result*power(i,last)
        return result
    return horner(list(poly.terms.items()),0)
