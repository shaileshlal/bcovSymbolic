"""Lossless arithmetic tokens for very large signed rational labels.

Do not feed decimal floats of thousand-digit integers to a transformer.
Digits (base 1000) + small-prime valuations separate magnitude and arithmetic.
"""
from fractions import Fraction as Q

def digits_base(n,base=1000):
    if base<2:raise ValueError('base>=2 required')
    n=abs(int(n));out=[]
    while n:out.append(n%base);n//=base
    return out or [0]

def unpack_digits(digits,base=1000):
    return sum(int(d)*base**i for i,d in enumerate(digits))

def valuations(n,primes=(2,3,5,7,11,13)):
    n=abs(int(n))
    if not n:return {str(p):0 for p in primes}
    result={}
    for p in primes:
        k=0
        while n%p==0:n//=p;k+=1
        result[str(p)]=k
    return result

def encode_fraction(f,base=1000):
    f=Q(f)
    return {'sign':(1 if f>0 else -1 if f<0 else 0),
            'num_be_digits':digits_base(f.numerator,base)[::-1],
            'den_be_digits':digits_base(f.denominator,base)[::-1],
            'num_valuations':valuations(f.numerator),
            'den_valuations':valuations(f.denominator),
            'base':base}

def decode_fraction(tok):
    base=int(tok['base'])
    p=unpack_digits(tok['num_be_digits'][::-1],base)
    d=unpack_digits(tok['den_be_digits'][::-1],base)
    return Q(tok['sign']*p,d)

def coefficient_record(g,j,coefficient,alpha=None,a0=None):
    """A complete, lossless record and independent physics labels."""
    c=Q(coefficient)
    return {'genus':g,'X_degree':j,'rational':encode_fraction(c),
            'alpha':None if alpha is None else encode_fraction(alpha),
            'constant':None if a0 is None else encode_fraction(a0),
            'relative_to_a0':None if a0 is None else encode_fraction(c/Q(a0))}
