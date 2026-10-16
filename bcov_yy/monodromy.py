"""Local monodromy representations; these matrices use DIFFERENT local bases.

They are not global integral monodromies in a common symplectic basis. Their
product must not be tested without a connection matrix. Boundary constraints
only need the local nilpotent/finite-order structures and polarization residue.
"""
from math import comb

MUM = tuple(tuple(comb(i,j) if j<=i else 0 for j in range(4)) for i in range(4))
CONIFOLD = ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,1,0,1))
ORBIFOLD = ((0,0,0,-1),(1,0,0,-1),(0,1,0,-1),(0,0,1,-1))
IDENTITY = tuple(tuple(int(i==j) for j in range(4)) for i in range(4))


def multiply(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)) for i in range(4))


def power(a,n):
    r=IDENTITY
    for _ in range(n): r=multiply(r,a)
    return r


def subtract(a,b):
    return tuple(tuple(a[i][j]-b[i][j] for j in range(4)) for i in range(4))


def validate():
    zero=((0,)*4,)*4
    n0=subtract(MUM,IDENTITY)
    nc=subtract(CONIFOLD,IDENTITY)
    assert power(n0,4)==zero and power(n0,3)!=zero
    assert power(nc,2)==zero and nc!=zero
    assert power(ORBIFOLD,5)==IDENTITY and ORBIFOLD!=IDENTITY
    return True
