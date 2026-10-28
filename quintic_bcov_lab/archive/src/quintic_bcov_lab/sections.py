"""Exact residual x(1-x) section and its X-basis lift in canonical YY gauge.

The reference-derived split is an *analysis of physical labels*, not an
independent prediction. Gap/orbifold/MUM kernel is checked separately.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as Q
from math import comb
from .arithmetic import fraction


def dimensions(g):
    if type(g) is not int or g<2:raise ValueError('genus >= 2')
    K=3*g-3;m=(K+4)//5;D=(2*g-2)//5
    assert m+D==g-1
    return K,m,D


def kernel(g,r):
    K,m,D=dimensions(g)
    if not (0<=r<D):raise ValueError('invalid residual index')
    out=[Q(0)]*(K+1)
    for j in range(m+1):out[1+r+j]=Q((-1)**j*comb(m,j))
    return tuple(out)


@dataclass(frozen=True)
class Residual:
    genus:int
    reference:tuple[Q,...]
    anchor:tuple[Q,...]
    coefficients:tuple[Q,...]
    x_bernstein:tuple[Q,...]
    x_power:tuple[Q,...]

    def validate(self):
        K,m,D=dimensions(self.genus)
        assert len(self.reference)==K+1 and len(self.coefficients)==D
        assert len(self.x_power)==D+1 and not self.x_power[0]
        assert all(self.anchor[j]==0 for j in range(m+1,self.genus))
        a=[Q(0)]*(K+1)
        for r,c in enumerate(self.coefficients):
            for j,b in enumerate(kernel(self.genus,r)):a[j]+=c*b
        assert tuple(z+b for z,b in zip(a,self.anchor))==self.reference
        for x in (Q(1,2),Q(1,3),Q(3,7)):
            X=x/(x-1)
            q1=(x-1)**(self.genus-1)*sum((a[j]*X**j for j in range(K+1)),Q(0))
            q2=sum((v*x**j for j,v in enumerate(self.x_power)),Q(0))
            assert q1==q2
        return True


def from_reference(g,coefficients):
    """Zero-pivot affine gauge from a physical reference; no oracle-free claim."""
    a=tuple(map(fraction,coefficients));K,m,D=dimensions(g)
    if len(a)!=K+1:raise ValueError('inconsistent reference length')
    c=[Q(0)]*D
    for r in range(D-1,-1,-1):
        p=m+1+r
        higher=sum((c[t]*kernel(g,t)[p] for t in range(r+1,D)),Q(0))
        c[r]=(-1)**m*(a[p]-higher)
    delta=[Q(0)]*(K+1)
    for r,cr in enumerate(c):
        for j,v in enumerate(kernel(g,r)):delta[j]+=cr*v
    anchor=tuple(av-v for av,v in zip(a,delta))
    b=tuple((-1)**(g-2-r)*cr for r,cr in enumerate(c))
    poly=[Q(0)]*(D+1)
    for r,br in enumerate(b):
        for j in range(D-r):poly[1+r+j]+=br*(-1)**j*comb(D-1-r,j)
    out=Residual(g,a,anchor,tuple(c),b,tuple(poly))
    assert out.validate()
    return out


def difference_from_charges(g,charge_delta,triangular_response):
    """Solve H c = delta_n with exact forward substitution.

    Responds only to a supplied, exactly computed triangular GV matrix.
    Never treats unknown charges as zero.
    """
    D=dimensions(g)[-1]
    n=list(map(fraction,charge_delta))
    if len(n)!=D or len(triangular_response)!=D:raise ValueError('invalid response dimensions')
    c=[]
    for i in range(D):
        row=list(map(fraction,triangular_response[i]))
        if len(row)!=D or any(row[j]!=0 for j in range(i+1,D)):
            raise ValueError('nontriangular matrix')
        v=row[i]
        if not v:raise ValueError('singular GV response')
        c.append((n[i]-sum(row[j]*c[j] for j in range(i)))/v)
    assert all(sum(fraction(triangular_response[i][j])*c[j] for j in range(D))==n[i] for i in range(D))
    return tuple(c)
