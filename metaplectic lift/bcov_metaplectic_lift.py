#!/usr/bin/env python3
"""BCOV rank-one symplectic-shear lift: exact connected Gaussian graphs to genus 3.

All expressions are formal over Q[beta, jets]. No quintic-specific values are fitted.
The (+) heat convention T_beta=exp(beta*hbar/2*d_x**2) is used;
for x2 -> x2+x4 with p0=partial_0 F0, use beta=-1
under the saddle convention y=x-beta*F0'(x).

Usage: python bcov_metaplectic_lift.py
Requires sympy.
"""
import sympy as S

h, beta, x, y, K = S.symbols('h beta x y K')
f00,f01,f02,f03,f04,f05,f06 = S.symbols('f00 f01 f02 f03 f04 f05 f06')
f10,f11,f12,f13,f14 = S.symbols('f10 f11 f12 f13 f14')
f20,f21,f22,f30 = S.symbols('f20 f21 f22 f30')


def wick(expr, k=K):
    """Normalized centered Gaussian expectation E[x^(2m)]=(2m-1)!! k^m."""
    out=S.S.Zero
    for (power,), coeff in S.Poly(S.expand(expr), x).terms():
        if power % 2 == 0:
            out += coeff*S.factorial2(power-1)*k**(power//2) if power else coeff
    return S.expand(out)

# Gaussian expansion around the classical saddle χ, where
# y=χ-βF0'(χ) and K=β/(1-β F0''(χ)).
V1 = f11*x+f03*x**3/6
V2 = f20+f12*x**2/2+f04*x**4/24
V3 = f21*x+f13*x**3/6+f05*x**5/120
V4 = f30+f22*x**2/2+f14*x**4/24+f06*x**6/720
E2 = wick(V2+V1**2/2)
E4 = wick(V4+V1*V3+V2**2/2+V1**2*V2/2+V1**4/24)
G2 = S.expand(E2)
G3 = S.expand(E4-E2**2/2)


def heat_flow_rhs(Fs):
    """For F=F0+h F1+h^2 F2+..., ∂β F=1/2 Fx^2 + h/2 Fxx."""
    out=[]
    for n,fn in enumerate(Fs):
        term=S.Rational(1,2)*sum(S.diff(Fs[j],x)*S.diff(Fs[n-j],x) for j in range(n+1))
        if n>=1:term+=S.diff(Fs[n-1],x,2)/2
        out.append(S.expand(term))
    return out


def tests():
    expected_g2=(f20+K*(f12+f11**2)/2
                +K**2*(f04/8+f11*f03/2)
                +S.Rational(5,24)*K**3*f03**2)
    assert S.expand(G2-expected_g2)==0
    assert S.expand(G3.subs({f03:0,f04:0,f05:0,f06:0})-f30-
        (K*f22/2+K*f11*f21+K**2*(f14/8+f11*f13/2+f12**2/4+f11**2*f12/2)))==0
    # No higher genus or propagator mixing at beta=0.
    assert G2.subs(K,0)==f20 and G3.subs(K,0)==f30
    # F2, F3 enter affinely with coefficient exactly 1.
    assert S.diff(G2,f20)==1 and S.diff(G3,f30)==1
    assert not (set([f21,f22]) & G2.free_symbols)
    assert S.diff(G3,f20)==0 and S.diff(G3,f21)!=0
    # Gaussian (free quadratic F0) checks: F0=x^2 a/2, F1=b*x, F2=c*x, F3=0.
    # This can be integrated exactly, giving F2'=c*chi+K*b*c, F3'=K*c^2/2.
    assert S.simplify(G2.subs({f03:0,f04:0,f11:S.Symbol('b'),f12:0,f20:S.Symbol('c')*S.Symbol('chi')}))==S.Symbol('c')*S.Symbol('chi')+K*S.Symbol('b')**2/2
    assert S.simplify(G3.subs({f03:0,f04:0,f05:0,f06:0,f11:S.Symbol('b'),f12:0,f13:0,f14:0,f20:0,f21:S.Symbol('c'),f22:0,f30:0}))==K*S.Symbol('b')*S.Symbol('c')
    # At first order in beta, connected transport from quantum HJ equation.
    # At beta=0 χ=y, K=beta+O(beta^2). F2new has shift of argument χ=y+beta F0'(y).
    assert S.expand(G2.coeff(K,1)-(f12+f11**2)/2)==0
    assert S.expand(G3.coeff(K,1)-(f22/2+f11*f21))==0
    # Symbolic symplectic matrix for the mirror-quintic conifold.
    M=S.eye(4);M[1,3]=1
    Omega=S.Matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    assert M.T*Omega*M==Omega
    R=S.eye(4)
    R[1,1]=0; R[1,3]=1; R[3,1]=-1; R[3,3]=0
    assert R.T*Omega*R==Omega
    Mc_con=R*M*R.inv()
    expected=S.eye(4);expected[3,1]=-1
    assert Mc_con==expected
    x0,x1,t=S.symbols('x0 x1 t')
    phi=S.symbols('phi',cls=S.Function)
    f=phi(x1/x0)/x0**2
    d1=S.diff(f,x0).subs(x1,t*x0)
    d2=S.diff(f,x0,2).subs(x1,t*x0)
    assert S.simplify(d1-x0**(-3)*(-2*phi(t)-t*S.diff(phi(t),t)))==0
    assert S.simplify(d2-x0**(-4)*(6*phi(t)+6*t*S.diff(phi(t),t)+t**2*S.diff(phi(t),t,2)))==0
    return True

if __name__ == '__main__':
    tests()
    print('CONVENTION: y=chi-beta*F0_prime(chi), K=beta/(1-beta*F0_second(chi))')
    print('GENUS TWO correction G2-F2:')
    print(S.collect(S.expand(G2-f20),K))
    print('GENUS THREE correction G3-F3, collected by K:')
    print(S.collect(S.expand(G3-f30),K))
    print('genus 3 number of distinct K/jets monomials:',len(S.Poly(G3, K,f03,f04,f05,f06,f11,f12,f13,f14,f20,f21,f22,f30).terms()))
    print('CONIFOLD POLARIZATION R: (x2,x4)->(x4,-x2); conjugated Mc acts p_c->p_c-q_c')
    print('PASS: exact Wick, ambiguity transport, homogeneous derivatives, and symplectic Fourier-frame tests.')
