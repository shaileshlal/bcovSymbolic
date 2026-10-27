"""Exact-fraction quintic HKQ holomorphic-limit q-series verification.
Reads uploaded ff_g, normalizes via propagator dictionary + coordinate inversion.
No floating-point arithmetic. SymPy only used for symbolic g<=4 construction.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import sympy as sp
N=6
zz=sp.symbols('z'); v1,v2,v3,X=sp.symbols('v1 v2 v3 X')

def trunc(a):return (list(a)+[Q(0)]*(N+1))[:N+1]
def c(n):return trunc([Q(n)])
def add(a,b):return [x+y for x,y in zip(a,b)]
def neg(a):return [-x for x in a]
def sub(a,b):return add(a,neg(b))
def scale(a,x):return [x*y for y in a]
def mul(a,b):return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N+1)]
def powp(a,n):
 r=c(1)
 for _ in range(n):r=mul(r,a)
 return r
def inv(a):
 assert a[0]
 res=[Q(0)]*(N+1);res[0]=1/a[0]
 for n in range(1,N+1):res[n]=-sum(a[k]*res[n-k] for k in range(1,n+1))/a[0]
 return res
def div(a,b):return mul(a,inv(b))
def theta(a):return [i*x for i,x in enumerate(a)]
def exp(a):
 assert not a[0]
 out=[Q(0)]*(N+1);out[0]=Q(1)
 for n in range(1,N+1):out[n]=sum(k*a[k]*out[n-k] for k in range(1,n+1))/n
 return out
def comp(a,b):
 out=c(0);bk=c(1)
 for ak in a:out=add(out,scale(bk,ak));bk=mul(bk,b)
 return out

def convert_g(g):
 z=zz;k=sp.Integer(5);mu=sp.Rational(1,3125)
 Cf,S2v,T1v,T0v=sp.symbols('Cf S2 tS1 tS0')
 t=(Path(__file__).resolve().parent/'data'/f'ff_{g}.txt').read_text()
 for old,new in [('Cf[z]','Cf'),('S2[z]','S2'),('tS1[z]','tS1'),('tS0[z]','tS0')]:t=t.replace(old,new)
 F=sp.sympify(t.replace('^','**'))
 S2=-z*z/(k*X)*(v1+sp.Rational(3,5))
 T1=z/(k*X)*(v2+sp.Rational(2,25))
 T0=1/(2*k*X)*(v3+(1-X)*(v2+sp.Rational(4,125)))
 C=k*X/z**3
 Pz=sp.Poly(sp.cancel((k*X)**(g-1)*F.subs({Cf:C,S2v:S2,T1v:T1,T0v:T0}).subs(z,mu*(X-1)/X)),v1,v2,v3,X)
 Ph=sp.Poly(sp.expand((-1)**(g-1)*Pz.as_expr().subs({v1:-v1,v2:v2,v3:-v3,X:1-X},simultaneous=True)),v1,v2,v3,X)
 return Ph

def ps(a):return ' + '.join(f'{x} q^{i}' for i,x in enumerate(a) if x)

w0=[Q(factorial(5*n),factorial(n)**5) for n in range(N+1)]
def H(n):return sum((Q(1,i) for i in range(1,n+1)),Q(0))
sigma=trunc([0]+[w0[n]*5*(H(5*n)-H(n)) for n in range(1,N+1)])
f=div(sigma,w0)
qoverz=exp(f)
# theta(t)=1+theta(sigma/w0)
theta_t=add(c(1),theta(f))
# HKQ holomorphic limit: B1 = -theta(w0)/w0, A = -1-theta(theta_t)/theta_t
u=neg(div(theta(w0),w0))
A=sub(neg(c(1)),div(theta(theta_t),theta_t))
w2=div(theta(theta(w0)),w0)
w3=neg(div(theta(theta(theta(w0))),w0))
V1=add(add(c(1),A),scale(u,Q(2)))
V2=sub(w2,mul(u,V1))
# x_HKQ = -3125*z/(1-3125*z)
xH=trunc([0]+[-Q(3125)**n for n in range(1,N+1)])
V3=add(sub(add(w3,mul(u,V2)),mul(mul(u,V1),xH)),scale(mul(u,xH),Q(2,5)))
Pvals=(V1,V2,V3,xH)
assert all(p[0]==0 for p in Pvals)
# ratio (1-psi)/(5psi) = -(1-3125z)/5
factor=scale(mul(powp(w0,2),trunc([1,-3125])),-Q(1,5))
# inverse q(z)=z exp(sigma/w0) by reversion
qofz=trunc([0]+qoverz[:N])
zin=[Q(0),Q(1)]+[Q(0)]*(N-1)
for n in range(2,N+1):
 err=comp(qofz,zin)[n]; zin[n]-=err
assert comp(qofz,zin)==trunc([0,1])
print('q(z) =',ps(qofz));print('z(q)=',ps(zin))
print('v1h:',ps(V1))
for g in (2,3,4):
 P=convert_g(g)
 powers=[[powp(vals,n) for n in range(g*3)] for vals in Pvals]
 result=c(0)
 for mon,co in P.terms():
  term=c(Q(int(sp.numer(co)),int(sp.denom(co))))
  for j,ex in enumerate(mon):term=mul(term,powers[j][ex])
  result=add(result,term)
 Fz=mul(powp(factor,g-1),result)
 Fq=comp(Fz,zin)
 print('\ng',g,'Pmonomials',len(P.terms()))
 print('F_g(q) through q^6:',ps(Fq))
 if g in (2,3):
  n0={1:2875,2:609250,3:317206375,4:242467530000,5:229305888887625,6:248249742118022000}
  n2={1:0,2:0,3:0,4:534750,5:75478987900,6:871708139638250}
  n3={1:0,2:0,3:0,4:8625,5:-15663750,6:3156446162875}
  expect=[Q(0)]*(N+1)
  for d in range(1,N+1):
   for k in range(1,N//d+1):
    if g==2:expect[d*k]+=Q(k)*(Q(n0[d],240)+n2[d])
    if g==3:expect[d*k]+=Q(k**3)*(Q(n0[d],6048)-Q(n2[d],12)+n3[d])
  print('GV-predicted coefficients:',ps(expect))
  print('GV residual:',ps(sub(Fq,add(expect,c(Fq[0])))))
  assert all(Fq[i]==expect[i] for i in range(1,N+1))
