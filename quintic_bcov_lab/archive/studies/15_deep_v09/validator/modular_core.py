"""Reconstruct quintic rank-zero DT modular series from its SEVEN polar terms.

Independent exact implementation of theta/Eisenstein/eta modular ansatz (paper eq 5.5).
No nonpolar h_mu coefficients or Q>20 GV/PT data enter the fit.
"""
from fractions import Fraction as F
from math import comb
import sympy as sp

N=14

def add(a,b):return [a[i]+b[i] for i in range(N+1)]
def mul(a,b):
 c=[F(0)]*(N+1)
 for i,ai in enumerate(a):
  if not ai:continue
  for j,bj in enumerate(b[:N+1-i]):
   if bj:c[i+j]+=ai*bj
 return c

def powser(a,k):
 c=[F(1)]+[F(0)]*N
 while k:
  if k&1: c=mul(c,a)
  a=mul(a,a);k//=2
 return c

def E(n):
 a=[F(1)]+[F(0)]*N
 coeff={2:-24,4:240,6:-504}[n]
 for t in range(1,N+1):
  a[t]=coeff*sum(d**(n-1) for d in range(1,t+1) if t%d==0)
 return a
E2,E4,E6=E(2),E(4),E(6)

def eta_inverse(p):
 a=[F(1)]+[F(0)]*N
 for m in range(1,N+1):
  b=[F(0)]*(N+1)
  for j in range(0,N//m+1):b[j*m]=F(comb(p+j-1,j))
  a=mul(a,b)
 return a
ETA=eta_inverse(70)

def theta(mu):
 # theta_mu^{(5)} = q^e * theta[], with e=(5,9/5,1/5)/8
 e=[F(5,8),F(9,40),F(1,40)][mu]
 a=[F(0)]*(N+1)
 for n in range(-20,21):
  # k = n+mu/5+5/2, 5k half-integer
  k=F(n,1)+F(mu,5)+F(5,2)
  power=F(5,2)*k*k
  relative=power-e
  assert relative.denominator==1,(mu,n,relative)
  if 0<=relative<=N:
   # (-i*5)*exp(i*pi*5*k)*k = (-1)^(5n+mu+12)*(5k)
   coefficient=(1 if (n+mu)%2==0 else -1)*5*k
   a[int(relative)]+=coefficient
 return e,a

def serre(e,a,weight):
 out=[F(0)]*(N+1)
 for i in range(N+1):out[i]=(e+i)*a[i]
 conv=mul(E2,a)
 return [out[i]-F(weight,12)*conv[i] for i in range(N+1)]

def modular_features(mu):
 e,t0=theta(mu)
 t1=serre(e,t0,F(3,2))
 t2=serre(e,t1,F(7,2))
 polys=[mul(powser(E4,8),[F(1)]+[F(0)]*N),
        mul(powser(E4,5),powser(E6,2)),
        mul(powser(E4,2),powser(E6,4)),
        mul(powser(E4,6),E6),
        mul(powser(E4,3),powser(E6,3)),
        powser(E6,5),
        powser(E4,7),
        mul(powser(E4,4),powser(E6,2)),
        mul(E4,powser(E6,4))]
 return e-F(70,24), [mul(ETA,mul(polys[i],[t0,t0,t0,t1,t1,t1,t2,t2,t2][i])) for i in range(9)]

POLAR={0:{0:5,1:-800,2:58500},
       1:{0:0,1:0,2:8625},
       2:{0:0,1:0,2:0}}
SERIES={mu:modular_features(mu) for mu in [0,1,2]}
rows=[];targets=[]
for mu,coeffs in POLAR.items():
 for n,v in coeffs.items():
  rows.append([sp.Rational(b[n].numerator,b[n].denominator) for b in SERIES[mu][1]])
  targets.append(v)
A=sp.Matrix(rows);b=sp.Matrix(targets)
# The seven polar values, plus the forced leading zeros in the odd-theta series,
# provide nine equations. Their exact coefficient matrix has full rank.
assert A.rank()==9 and A.row_join(b).rank()==9
vec=[sp.simplify(x) for x in A.inv()*b]
assert all(not expr.free_symbols for expr in vec)

def h_series(mu):
    e,feat=SERIES[mu]
    coef=[]
    for n in range(N+1):
        v=sum(vec[i]*sp.Rational(feat[i][n].numerator,feat[i][n].denominator) for i in range(9))
        assert v.q==1,(mu,n,v)
        coef.append(int(v))
    return e,coef

# Purely evaluation-only oracle, not used when solving for the modular form.
EXPE={0:[5,-800,58500,5817125,75474060100,28096675153255],
      1:[0,0,8625,-1138500,3777474000,3102750380125,577727215123000],
      2:[0,0,0,-1218500,441969250,953712511250]}

def verify_postfit():
    for mu,known in EXPE.items():
        e,got=h_series(mu)
        assert got[:len(known)]==known,(mu,got[:len(known)],known)
    return True

if __name__=="__main__":
    print('Rank of nine rational polar equations:',A.rank())
    print('Exact modular basis coefficients:',list(map(str,vec)))
    for mu in [0,1,2]:
        e,got=h_series(mu)
        print(f'h_{mu}(q) = q^({e})*',got[:8])
    print('Independent nonpolar comparison:',verify_postfit())
