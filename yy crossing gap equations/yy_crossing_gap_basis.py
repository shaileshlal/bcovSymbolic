#!/usr/bin/env python3
"""Exact two-boundary ('Polyakov-like') basis for the quintic YY ambiguity.

Requires sympy plus yy_p2_p3_exact.py (the previously verified source module).
All calculations are over Q.  This is NOT the full nonlinear BCOV crossing map:
- The conifold pole filtration uses a local frame A_g=Delta^(g-1)*U_g,
  and Delta=1-3125*z; the kernel is exact for any analytic unit U_g.
- A pole-DUAL basis with exact conifold t_c and U_g needs their local jets;
  a symbolic genus-two extension is supplied.
- The MUM q-jets and GV-dual basis are fully exact in the LCS frame.
- D4-D2-D0 modular forms and the PT/DT transfer are not implemented.
"""
from fractions import Fraction
from math import comb
from pathlib import Path
import sympy as sp
import yy_p2_p3_exact as yy

X = yy.X
Delta = sp.symbols('Delta')
xi = sp.symbols('xi')


def rational(f):
    return sp.Rational(f.numerator, f.denominator)


def rational_poly_coeff(g, r, n):
    """q^0..q^n coefficients of E_g*(1-X)^r in the LCS frame.

    E_g = omega0(z)^(2g-2)*(1-3125*z)^(g-1),
    X=-3125*xi, xi=z/(1-3125*z), q=z+O(z^2).
    """
    m = n + 3
    q, Ebase = yy.mirror_data(m)
    E = yy.pw(Ebase, g-1)
    coeffs = [Fraction(0)] * (m+1)
    for k in range(r+1):
        coeffs[k] = Fraction(comb(r,k)*3125**k)
    result = yy.compose(yy.mul(E, coeffs), yy.revert(q))
    return [rational(c) for c in result[:n+1]]


def basis(g):
    assert g >= 2
    m = 3*g-3
    # The gap's homogeneous kernel is span{(1-X)^r, r = 0,...,g-1}.
    # This remains true for every analytic nonzero unit U_g and every
    # conifold uniformizer t_c=const*Delta + O(Delta^2).
    T = sp.Matrix(g, g, lambda d,r: rational_poly_coeff(g, r, g-1)[d])
    Ti = T.inv()
    K = [sp.expand(sum(Ti[r,d]*(1-X)**r for r in range(g))) for d in range(g)]
    # In the first g Fourier coefficients these are unit vectors.
    for d in range(g):
        v = sp.Matrix([sum(T[k,r]*Ti[r,d] for r in range(g)) for k in range(g)])
        assert v == sp.eye(g)[:,d]
    C = [K[0]]
    for d in range(1,g):
        C.append(sp.expand(sum(r**(2*g-3)*K[r*d] for r in range(1,(g-1)//d+1))))
    # 'C' is GV dual: after the standard genus-g multiple-cover sieve
    # on the first g-1 positive Fourier coefficients, delta n_{g,d}=1.
    for d in range(1,g):
        # Explicit evaluation through the closed basis-change T.
        power = sp.Poly(C[d].subs(X,1-Delta),Delta)
        coef = [power.nth(r) for r in range(g)]
        qjet = T*sp.Matrix(coef)
        for e in range(1,g):
            actual = sum(sp.mobius(r)*r**(2*g-3)*qjet[e//r] for r in sp.divisors(e))
            assert actual == int(e==d)
        assert qjet[0] == 0
    return T, Ti, K, C


def gap_pole_basis(g):
    """Boundary pole modes in the leading discriminant model U=1, t_c=Delta.

    Delta^(g-1) * (1-X)^(g-1+j) = Delta^(-j), j=1,...,2g-2.
    """
    return [sp.expand((1-X)**(g-1+j)) for j in range(1,2*g-1)]


def exact_conifold_pole_basis_genus_two():
    """Retain unknown exact conifold connection jets as symbols.

    U(Delta)=1+u1 Delta+u2 Delta^2+...;
    t_c(Delta)=Delta (1+c1 Delta+c2 Delta^2+...).
    B_{2,j}=PP_Delta[Delta^(-1) U^-1 t_c^(-j)].
    """
    u1,u2,c1,c2=sp.symbols('u_1 u_2 c_1 c_2')
    U = 1+u1*Delta+u2*Delta**2
    t = Delta*(1+c1*Delta+c2*Delta**2)
    out=[]
    for j in (1,2):
        h=j+1
        unit=sp.series(1/U*(t/Delta)**(-j),Delta,0,h).removeO().expand()
        B=sum(sp.expand(unit).coeff(Delta,k)*Delta**(k-h) for k in range(h))
        out.append(sp.expand(B))
        # Exact symbolic gap-dual condition to given Laurent order.
        expr=sp.series(Delta*U*B-t**(-j),Delta,0,0).removeO().expand()
        assert all(expr.coeff(Delta,-k)==0 for k in (1,2))
    assert sp.simplify(out[0] -(Delta**(-2) -(u1+c1)/Delta)) == 0
    assert sp.simplify(out[1] - (Delta**(-3) -(u1+2*c1)*Delta**(-2)
                      +(u1**2-u2+2*u1*c1+3*c1**2-2*c2)*Delta**(-1))) == 0
    return out


def projection_data(g, T, K, C):
    f=sp.Poly((yy.P2 if g==2 else yy.P3).subs({yy.v1:0, yy.v2:0, yy.v3:0}),X).as_expr()
    p=sp.Poly(sp.expand(f.subs(X,1-Delta)),Delta)
    polar=sp.expand(sum(p.nth(r)*(1-X)**r for r in range(g,3*g-2)))
    ker=sp.expand(sum(p.nth(r)*(1-X)**r for r in range(g)))
    assert sp.expand(ker+polar-f)==0
    qvec=T*sp.Matrix([p.nth(r) for r in range(g)])
    gvvec=sp.Matrix([qvec[0]]+[sum(sp.mobius(r)*r**(2*g-3)*qvec[d//r] for r in sp.divisors(d)) for d in range(1,g)])
    assert sp.expand(ker-sum(gvvec[d]*C[d] for d in range(g)))==0
    return f,polar,ker,[p.nth(r) for r in range(3*g-2)],list(gvvec)


def run():
    lines=[]
    def emit(t=''):
        lines.append(str(t));print(t)
    emit('TWO-BOUNDARY CONIFOLD/LCS BASIS, EXACT RATIONAL ARITHMETIC')
    emit('Convention: X=1/(1-psi^5)=1-Delta^-1, Delta=1-3125*z.')
    emit('Normalized LCS amplitude is N_g * E_g * f_g, N_g=(-1/625)^(g-1).')
    emit('Conifold exact gap KERNEL is span[(1-X)^r, 0<=r<=g-1].')
    emit('The simple conifold pole DUAL basis below assumes t_c=Delta,U_g=1.')
    for g in (2,3):
        m=3*g-3
        T,Ti,K,C=basis(g)
        emit('\n=== GENUS %d ===' % g)
        emit('ambiguity dim=%d; conifold gap rank=%d, kernel dim=%d'%(m+1,2*g-2,g))
        emit('single-MUM filtration quotient dims (degrees 0..m): %s'%([1]*(m+1)))
        emit('joint gap-kernel / MUM filtration dims: %s'%([g-d for d in range(g+1)]))
        emit('MUM q-jet matrix (rows degree 0..g-1, columns (1-X)^r):')
        emit(T)
        emit('det T=%s'%T.det())
        emit('MUM-crossing dual K[g,d] (gap regular, q^d dual):')
        for d,k in enumerate(K): emit('K[%d,%d] = %s'%(g,d,k))
        emit('GV-crossing dual C[g,d] (gap regular, primitive n_g(d) dual):')
        for d,c in enumerate(C): emit('C[%d,%d] = %s'%(g,d,c))
        emit('leading conifold pole modes:')
        for j,pole in enumerate(gap_pole_basis(g),1):
            emit('j=%d: (1-X)^%d   =>   Delta^(-%d)'%(j,g-1+j,j))
        f,pol,ker,coeff,coords=projection_data(g,T,K,C)
        emit('pure-X projection f_g(X): %s'%f)
        emit('decomposition f_g = gap-polar part + gap-regular part:')
        emit('pole-mode coefficients r=0..m: %s'%coeff)
        emit('gap-polar part = %s'%sp.factor(pol))
        emit('gap-regular part = %s'%sp.factor(ker))
        emit('gap-regular coordinates in GV-crossing basis C[g,0..g-1]: %s'%coords)
        # Direct check of joint filtration K: K_d has q^0..q^(g-1)=delta_{jd}
        for d,k in enumerate(K):
            jet=sp.Matrix([sum(sp.Poly(k.subs(X,1-Delta),Delta).nth(r)*T[e,r] for r in range(g)) for e in range(g)])
            assert jet==sp.eye(g)[:,d]
    pole=exact_conifold_pole_basis_genus_two()
    emit('\nExact symbolic genus-two pole DUAL correction with arbitrary conifold jets:')
    for j,b in enumerate(pole,1):emit('B[2,%d](Delta) = %s'%(j,b))
    emit('\nAll identities and exact jet / Möbius checks PASS.')
    Path('/mnt/data/yy_crossing_gap_results.txt').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':
    run()
