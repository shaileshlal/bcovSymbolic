"""Gu/Klemm gauge identities as formal algebraic identities; optional SymPy dependency."""
import sympy as s
z,X,v1,v2,v3=s.symbols('z X v1 v2 v3')
K=s.Integer(5)
DX=X*(X-1)
Dv1=-2*v2-v1**2-v1+X*v1+s.Rational(2,5)*(X-1)
Dv2=v3-v1*v2
D=lambda f:z*s.diff(f,z)+DX*s.diff(f,X)+Dv1*s.diff(f,v1)+Dv2*s.diff(f,v2)
C=K*X/z**3
Szz=-z**2/(K*X)*(v1+s.Rational(3,5))
Sz=z/(K*X)*(v2+s.Rational(2,25))
St=(v3+(1-X)*(v2+s.Rational(4,125)))/(2*K*X)
qzz=-s.Rational(8,5)/z
qzss=z/25
qzzlower=s.Rational(2,25)/z**2
qzzup=-s.Rational(6,625)
r1=s.cancel((D(Szz)/z-C*Szz**2+2*qzz*Szz-qzss)/2-Sz)
r2=s.cancel((D(Sz)/z-C*Sz*Szz+qzz*Sz+qzzlower*Szz-qzzup)/2-St)
assert r1==0 and r2==0, (r1,r2)
print('PASS: Gu/Klemm q-gauge propagator differential identities')
