#!/usr/bin/env python3
"""Exact mirror-quintic monodromy/YY-generator connection action.

Conventions: s=3125*z, theta=s d/ds, D=psi d/dpsi=-5 theta,
X=1-1/(1-s). Period vector (x1,x2,x3,x4)^T uses Movasati's
integral symplectic basis, with x2 normalized to omega_0.

Demonstrates:
* symplectic MUM/conifold monodromy, Gepner fifth-order product;
* exact exterior-square action and invariant quotient;
* explicit semilinear conifold action on J1,J2,J3,H1,H2,H3;
* covariant action on genus-2 and genus-3 YY polynomials;
* finite jet / composition laws for S = W' / W.

Requires sympy and sibling yy_p2_p3_exact.py. Run as python yy_connection_action.py.
"""
from itertools import combinations
from pathlib import Path
import sympy as s
import yy_p2_p3_exact as yy

M0 = s.Matrix([[1,1,0,0],[0,1,0,0],[5,5,1,0],[0,-5,-1,1]])
Mc = s.eye(4); Mc[1,3] = 1
Omega=s.Matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
pairs=list(combinations(range(4),2))

def wedge2(M):
    return s.Matrix([[M[i,a]*M[j,b]-M[i,b]*M[j,a]
        for a,b in pairs] for i,j in pairs])

svar=s.symbols('s')
companion=s.Matrix([[0,0,0,svar*s.Rational(24,625)/(1-svar)],
                    [1,0,0,svar*s.Rational(2,5)/(1-svar)],
                    [0,1,0,svar*s.Rational(7,5)/(1-svar)],
                    [0,0,1,svar*2/(1-svar)]])
# theta Xmat = Xmat * companion, Xmat[i,j] = theta**j xi

J1,J2,J3,X=yy.j1,yy.j2,yy.j3,yy.X
A,B,C=s.symbols('sigma1 sigma2 sigma3') # sigma_k=D^k S / S
Jchange={J1:J1+A, J2:J2+2*A*J1+B,
         J3:J3+3*A*J2+3*B*J1+C}

def connection(P):
    """Conifold-continued formal polynomial via the Wronskian cocycle."""
    return s.expand(P.subs(Jchange, simultaneous=True))

# v's in original YY variables, expressed by J's. Source: YY ring Dmap.
v1=J1
v2=(5*X*J1-10*X-J2)/2
v3=(-J3+5*X*J2+(25*X**2-35*X)*J1+50*X-50*X**2)/2
V=(v1,v2,v3)
Vnew=tuple(connection(e) for e in V)
P2=s.Poly(yy.P2J,J1,J2,J3,X).as_expr()
P3=s.Poly(yy.P3J,J1,J2,J3,X).as_expr()


def assert_exact():
    assert M0.T*Omega*M0==Omega and Mc.T*Omega*Mc==Omega
    assert (M0*Mc)**5==s.eye(4)
    assert (Mc-s.eye(4))**2==s.zeros(4)
    assert wedge2(M0*Mc)==wedge2(M0)*wedge2(Mc)
    basisfixed=(wedge2(M0)-s.eye(6)).col_join(wedge2(Mc)-s.eye(6)).nullspace()
    assert basisfixed==[s.Matrix([0,1,0,0,1,0])]
    assert (V[0],s.simplify(V[1].subs({J1:yy.v1,J2:yy.J2})-yy.v2),
        s.simplify(V[2].subs({J1:yy.v1,J2:yy.J2,J3:yy.J3})-yy.v3))==(J1,0,0)
    # The genus polynomials commute with changing representation.
    for Pg,original in [(P2,yy.P2),(P3,yy.P3)]:
        original_J=s.expand(original.subs({yy.v1:Vnew[0],yy.v2:Vnew[1],yy.v3:Vnew[2]},simultaneous=True))
        assert s.expand(connection(Pg)-original_J)==0
    # M0: W12 is fixed and X is unchanged, hence all J_i fixed.
    assert wedge2(M0)[0,:]==s.Matrix([[1,0,0,0,0,0]])
    # Mc: W12 -> W12 + W14, W14 -> W14.
    assert wedge2(Mc)[0,:]==s.Matrix([[1,0,1,0,0,0]])
    assert wedge2(Mc)[2,:]==s.Matrix([[0,0,1,0,0,0]])
    # All pure-X ambiguity polynomials are left alone by raw holomorphic ring continuation.
    assert connection(s.Symbol('b0')+s.Symbol('b1')*X+s.Symbol('b2')*X**2)==s.Symbol('b0')+s.Symbol('b1')*X+s.Symbol('b2')*X**2
    # Triple derivative connection follows product rule and is triangular.
    assert s.expand(connection(J3)-(J3+3*A*J2+3*B*J1+C))==0


def summary():
    assert_exact()
    d2=s.Poly(connection(P2)-P2, A,B,C)
    d3=s.Poly(connection(P3)-P3, A,B,C)
    print('SYMLECTIC: M0,Mc preserve Omega; (M0 Mc)^5=1; (Mc-I)^2=0: PASS')
    print('GAUSS-MANIN theta companion =\n', companion)
    print('wedge basis indices (12,13,14,23,24,34)')
    print('wedge M0=\n',wedge2(M0))
    print('wedge Mc=\n',wedge2(Mc))
    print('simultaneous fixed space = span(W13+W24), identically zero on Griffiths-transversal periods')
    print('Holomorphic J-action: J1->J1+A, J2->J2+2*A*J1+B, J3->J3+3*A*J2+3*B*J1+C')
    print('YY v-action differences:')
    for i, diff in enumerate(Vnew):print('v%s correction:'%(i+1),s.factor(diff-V[i]))
    print('Genus-2 connection cocycle P2_new-P2 (factored over sigma):')
    print(s.collect(d2.as_expr(),[A,B,C]))
    print('Genus-3 cocycle degrees in sigma:', d3.total_degree(), 'sigma monomials',len(d3.terms()))
    print('genus 3 sigma3 derivative =', s.factor(s.diff(P3,J3)))
    print('genus 3 coefficient sigma3^2 =',d3.coeff_monomial(C**2))
    print('Tests passed: exact period monodromy, exterior square, transformation of P2 and P3.')
    return d2,d3

if __name__=='__main__':
    summary()
