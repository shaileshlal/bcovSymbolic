"""Certified arithmetic decoding: rational congruence + integer coset search.

Given n=b+R*a and a=m/D, find ALL m such that n is an integer, then
choose one nearest a proposed BPS index or minimum |n|. D is a genuine
independent assumption; NEVER derive D from a held-out answer.
"""
from fractions import Fraction as Q
from math import gcd
from dataclasses import dataclass

@dataclass(frozen=True)
class IntegralCoset:
    m0: int
    m_step: int
    n0: int
    n_step: int
    denominator: int

    def member(self,k):
        return Q(self.m0+self.m_step*k,self.denominator),self.n0+self.n_step*k

    def closest_index(self,predicted_n=0):
        if self.n_step==0:return self.member(0)
        goal=(Q(predicted_n)-self.n0)/self.n_step
        k=goal.numerator//goal.denominator
        options=[self.member(k),self.member(k+1)]
        return min(options,key=lambda z:(abs(z[1]-predicted_n),abs(z[1]),abs(z[0])))


def integral_coset(intercept,response,denominator):
    b=Q(intercept);R=Q(response);D=int(denominator)
    if D<=0 or R==0:raise ValueError('D>0 and R nonzero required')
    # b+R*m/D in Z. Exact congruence as in the genus-2/3 audit.
    A=R.numerator*b.denominator
    B=-b.numerator*D*R.denominator
    mod=b.denominator*D*R.denominator
    common=gcd(A,mod)
    if B%common:raise ValueError('No integer GV in this proposed coefficient lattice')
    period=mod//common
    m0=((B//common)*pow(A//common,-1,period))%period if period>1 else 0
    n0=b+R*Q(m0,D);step=R*Q(period,D)
    assert n0.denominator==1 and step.denominator==1
    return IntegralCoset(m0,period,int(n0),int(step),D)


def quintic_response(g,d):
    return (-1)**(g+d-1)*Q(5)**(5*d-4*g+4)
