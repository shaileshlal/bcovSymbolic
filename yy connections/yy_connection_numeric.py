#!/usr/bin/env python3
"""High-precision LCS->conifold Frobenius connection for mirror quintic.

Local coordinate t=1-s, s=3125*z.  Frobenius local basis:
 f0=1+0*t+0*t**2+...,
 f1=t+7/10*t**2+... (vanishing conifold period),
 f2=t**2+...,
 fl=(f1*log(t)+g(t))/(2*pi*i), g(0)=g'(0)=g''(0)=0,
 with L g = - L(f1*log t).

The analytic basis is fixed by the first three Taylor coefficients, while
fl undergoes fl -> fl + f1 around the conifold.

Connection x=C*(f0,f1,f2,fl), where x is Movasati's integral symplectic
basis rescaled globally to x2 = fundamental hypergeometric period.

Numerics use high-precision Taylor series at overlapping regular point s=1/2;
all local recurrences are rational. M_c = C L_c C^{-1} tests the result.
"""
from fractions import Fraction
from pathlib import Path
import mpmath as m
import sympy as s
from yy_connection_action import Mc

alpha,t=s.symbols('alpha t'); s_var=1-t
# theta t^alpha = -alpha t^(alpha-1)+alpha t^alpha.
def theta_shift(d):
    result={}
    for k,v in d.items():
        result[k-1]=result.get(k-1,0)-(alpha+k)*v
        result[k]=result.get(k,0)+(alpha+k)*v
    return {k:s.expand(v) for k,v in result.items()}
th=[{0:s.Integer(1)}]
for j in range(4):th.append(theta_shift(th[-1]))
L={}
def add_shift(d,sh,fac):
    for k,v in d.items():L[k+sh]=L.get(k+sh,0)+fac*v
add_shift(th[4],1,1) # (1-s)*theta^4 = t*theta^4
for k,fac in [(3,s.Integer(2)),(2,s.Rational(7,5)),(1,s.Rational(2,5)),(0,s.Rational(24,625))]:
    add_shift(th[k],0,-fac);add_shift(th[k],1,fac)
L={k:s.factor(v) for k,v in L.items()}
Ld={k:s.diff(v,alpha) for k,v in L.items()}
assert s.expand(L[-3]-alpha*(alpha-2)*(alpha-1)**2)==0

def fr(q):return Fraction(int(s.numer(q)),int(s.denom(q)))

def ratcoeff(nmax):
    """Exact power series coefficients for three analytic and one log solution."""
    shift={k:[fr(v.subs(alpha,n)) for n in range(nmax+1)] for k,v in L.items()}
    dshift={k:[fr(v.subs(alpha,n)) for n in range(nmax+1)] for k,v in Ld.items()}
    def analytic(init):
        a=[Fraction(0)]*nmax
        a[:3]=[Fraction(x) for x in init]
        for n in range(3,nmax):
            num=sum((a[j]*shift[n-j-3][j] for j in range(max(0,n-4),n)))
            a[n]=-num/shift[-3][n]
        return a
    a0=analytic([1,0,0]); a1=analytic([0,1,Fraction(7,10)]);a2=analytic([0,0,1])
    g=[Fraction(0)]*nmax
    for n in range(3,nmax):
        num=sum(g[j]*shift[n-j-3][j] for j in range(max(0,n-4),n))
        forcing=sum(a1[j]*dshift[n-j-3][j] for j in range(max(0,n-4),n+1))
        g[n]=-(num+forcing)/shift[-3][n]
    for n in range(3):
        assert sum(g[j]*shift[n-j-3][j] for j in range(max(0,n-4),n)) + sum(a1[j]*dshift[n-j-3][j] for j in range(max(0,n-4),n+1))==0
    return a0,a1,a2,g

def mpq(q):return m.mpf(q.numerator)/q.denominator

def jet_t(series,t_value,logpart=None):
    """theta^k jets at t, optionally add logpart(t)*log(t)."""
    t=m.mpf(t_value); lnt=m.log(t)
    out=[m.mpc(0)]*4
    for n in range(len(series)):
        a=mpq(series[n]);b=mpq(logpart[n]) if logpart else m.mpf(0)
        cs={0:m.mpf(1)};ds={0:m.mpf(0)}
        for k in range(4):
            out[k]+=sum((a*cs[sh]+b*(ds[sh]+cs[sh]*lnt))*t**(n+sh) for sh in cs)
            csnew={};dsnew={}
            for sh,c in cs.items():
                dc=ds[sh];deg=n+sh
                csnew[sh-1]=csnew.get(sh-1,0)-deg*c
                dsnew[sh-1]=dsnew.get(sh-1,0)-deg*dc-c
                csnew[sh]=csnew.get(sh,0)+deg*c
                dsnew[sh]=dsnew.get(sh,0)+deg*dc+c
            cs,ds=csnew,dsnew
    return out

def local_matrix(nmax,tval):
    f0,f1,f2,g=ratcoeff(nmax)
    v=[jet_t(f0,tval),jet_t(f1,tval),jet_t(f2,tval),
       [e/(2*m.pi*1j) for e in jet_t(g,tval,logpart=f1)]]
    return m.matrix(v)

def pmul(a,b):
    o=[m.mpc(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b[:4-i]):o[i+j]+=x*y
    return o

def period_matrix(nmax,sval):
    """Movasati x_i and theta^k x_i via epsilon-Frobenius at LCS."""
    sv=m.mpf(sval);logz=m.log(sv/3125); a=[m.mpf(1),m.mpf(0),m.mpf(0),m.mpf(0)]
    expz=[logz**k/m.factorial(k) for k in range(4)]
    psijets=m.matrix(4,4)
    for n in range(nmax):
        for k in range(4):
            poly=pmul(a,expz)
            nk=[m.mpf(0)]*4
            for j in range(k+1):nk[j]=m.binomial(k,j)*n**(k-j)
            val=pmul(poly,nk)
            for i in range(4):psijets[i,k]+=sv**n*val[i]
        j=n+1
        ak=[m.mpf(j-1)+m.mpf(k)/5 for k in (1,2,3,4)]
        r0=m.fprod(ak)/j**4
        l1=sum(1/x for x in ak)-4/m.mpf(j)
        l2=-m.mpf('0.5')*(sum(1/x**2 for x in ak)-4/m.mpf(j)**2)
        l3=(sum(1/x**3 for x in ak)-4/m.mpf(j)**3)/3
        rat=[r0,r0*l1,r0*(l2+l1**2/2),r0*(l3+l1*l2+l1**3/6)]
        a=pmul(a,rat)
    # Transform Frobenius psi_j (i=0...3) to x, globally x2=psi0.
    p2i=2*m.pi*1j
    C=m.matrix([[0,1/p2i,0,0],
                [1,0,0,0],
                [-m.mpf(25)/12,m.mpf(5)/(2*p2i),5/p2i**2,0],
                [200*m.zeta(3)/p2i**3,-m.mpf(25)/(12*p2i),0,-5/p2i**3]])
    return C*psijets

def compute(dps=65,nterms=250):
    m.mp.dps=dps
    Vlocal=local_matrix(nterms,m.mpf('0.5'))
    Vmov=period_matrix(nterms,m.mpf('0.5'))
    C=Vmov*Vlocal**-1
    T=m.eye(4);T[3,1]=1
    Ms=C*T*C**-1
    err=max(abs(Ms[i,j]-Mc[i,j]) for i in range(4) for j in range(4))
    # M1 only x2 has log part, equals x4/f1. Thus C[1,3]=C[3,1].
    err2=abs(C[1,3]-C[3,1])
    analytic_err=max(abs(C[i,3]) for i in (0,2,3))
    digits=max(0,int(-m.log10(err)))
    return C,err,err2,analytic_err,digits

if __name__=='__main__':
    C,err,err2,ae,digits=compute()
    lines=['Numerical LCS -> conifold connection (n=250, 65 dps, s=1/2)',
           'Periods x=C*(f0,f1,f2,f_log) with f_log=(f1 log t+g)/(2pi i)',
           'C =']
    for i in range(4):
        lines.append('  ['+',  '.join(m.nstr(C[i,j],18) for j in range(4))+']')
    lines+=['max |C*T_c*C^-1 - Mc| = '+m.nstr(err,9),
            '|C[x2,flog]-C[x4,f1]| = '+m.nstr(err2,9),
            'max unexpected log component = '+m.nstr(ae,9),
            'estimated good digits = '+str(digits)]
    Path(__file__).with_name('yy_connection_numeric_results.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
