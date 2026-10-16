#!/usr/bin/env python3
"""Independent differential-ring closure test at an ordinary point.

Compute the period jet up to D^5 using the quintic PF equation. Construct
W=x_2 D x_1-x_1 D x_2 and its D^0,...,D^4 derivatives, both before
and after exact conifold monodromy. Test the nontrivial fourth Bell-jet
relation implied by Yamaguchi--Yau's Dv3 identity. This verifies that
holomorphic monodromy acts within the YY differential ring (localized by W).
"""
import math
import mpmath as mp
import sympy as sp
from yy_connection_numeric import period_matrix
from yy_connection_action import Mc


def check(dps=55,nterms=180):
    mp.mp.dps=dps
    z=sp.symbols('s')
    coeff=[sp.Rational(24,625)*z/(1-z), sp.Rational(2,5)*z/(1-z),
           sp.Rational(7,5)*z/(1-z),2*z/(1-z)]
    operators=[]
    for n in range(6):
        if n<4:
            operators.append([sp.Integer(n==j) for j in range(4)])
        else:
            prev=operators[-1]
            operators.append([sp.cancel(z*sp.diff(prev[j],z)
                +(prev[j-1] if j else 0)+prev[3]*coeff[j]) for j in range(4)])
    V=period_matrix(nterms,mp.mpf('0.5'))
    jets=mp.matrix(4,6)
    for n,op in enumerate(operators):
        for i in range(4):
            jets[i,n]=(-5)**n*sum(mp.mpf(str(p.subs(z,sp.Rational(1,2))))*V[i,j] for j,p in enumerate(op))
    def bell_jets(P):
        a=list(P[0,:]);b=list(P[1,:])
        w=[sum(math.comb(n,j)*(b[j]*a[n-j+1]-a[j]*b[n-j+1]) for j in range(n+1)) for n in range(5)]
        j1,j2,j3,j4=[w[n]/w[0] for n in range(1,5)]
        x=mp.mpf(-1) # X=1-(1-s)^-1, s=1/2
        j4_formula=-(25*j1**2*x**2-70*j1**2*x+10*j1*j2*x-2*j1*j3
            -100*j1*x**2-350*j1*x+j2**2+190*j2*x-30*j3*x
            +100*x**2+404*x)/2
        return [j1,j2,j3,j4],abs(j4-j4_formula)
    J,err0=bell_jets(jets)
    M=mp.matrix([[mp.mpc(str(Mc[i,j])) for j in range(4)] for i in range(4)])
    Jc,err1=bell_jets(M*jets)
    assert err0 < mp.mpf('1e-45') and err1 < mp.mpf('1e-45')
    print('Bell J4 closure at LCS branch:',mp.nstr(err0,8))
    print('Bell J4 closure at conifold-continued branch:',mp.nstr(err1,8))
    print('J1 before and after:',mp.nstr(J[0],15),mp.nstr(Jc[0],15))
    return err0,err1

if __name__=='__main__':check()
