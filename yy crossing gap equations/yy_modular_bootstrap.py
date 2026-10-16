#!/usr/bin/env python3
"""Exact quintic YY -> GV bootstrap-adapted bases, genus 2 and 3.

Run alongside yy_p2_p3_exact.py:
    python yy_modular_bootstrap.py

WARNING: This constructs the *pullback* of published modular-bootstrap GV data.
It does not implement the D4-D2-D0 modular-form / PT wall-crossing computation.
Original quintic YY convention X=1/(1-psi^5), xi=-X/3125.
All arithmetic is rational; no numerical approximation.
"""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from sympy import symbols, Poly, expand, mobius, divisors, factor
import yy_p2_p3_exact as yy

xi = symbols('xi')
N0 = lambda g: Q((-1)**(g-1), 625**(g-1))
bootstrap_table = {
 2: {1:0,2:0,3:0,4:534750,5:75478987900,6:871708139638250},
 3: {1:0,2:0,3:0,4:8625,5:-15663750,6:3156446162875},
}
N0_gv = {1:2875, 2:609250,3:317206375,4:242467530000,
         5:229305888887625,6:248249742118022000}

def poly_from_series(series, var=xi):
    return sum(Q(c) * var**i for i,c in enumerate(series) if c)

def truncate(series, n):
    return series[:n]

@lru_cache(None)
def source_series(g, cutoff):
    """Build sufficiently high order to differentiate without loss of last term."""
    q, Ebase = yy.mirror_data(cutoff + 2)
    q=truncate(q,cutoff+1)
    Ebase=truncate(Ebase,cutoff+1)
    E=yy.pw(Ebase,g-1)
    qi=yy.revert(q)
    return q,E,qi

def gv_basis(g, cutoff):
    """Unscaled GV-dual basis: delta n_{g,d} = N0(g) * b_d.

    B[g,d](xi) = [ E_g(xi)^(-1) sum_{r>=1} r^(2g-3)q(xi)^(rd) ]_<=cutoff.
    For d=0 use Fourier constant-map dual B[g,0] = E_g^-1.
    """
    q,E,qi=source_series(g,cutoff)
    Psi=[yy.mul(yy.pw(q,k),yy.inv(E)) for k in range(cutoff+1)]
    B=[Psi[0]]
    for d in range(1, cutoff+1):
        coef = yy.zero(cutoff+1)
        for r in range(1,cutoff//d+1):
            coef=yy.add(coef,yy.scale(Psi[d*r],Q(r**(2*g-3))))
        B.append(coef)
    assert yy.mul(E,B[0]) == [Q(1)]+yy.zero(cutoff)
    for d in range(1,cutoff+1):
        lhs=yy.compose(yy.mul(E,B[d]),qi)
        expect=yy.zero(cutoff+1)
        for r in range(1,cutoff//d+1):
            expect[r*d]=Q(r**(2*g-3))
        assert lhs == expect, (g,d,lhs,expect)
    return B

def ambiguity_projection(g):
    P=yy.P2 if g==2 else yy.P3
    return P.subs({yy.v1:0, yy.v2:0, yy.v3:0})

def dual_coordinates(g):
    cutoff=3*g-3
    B=gv_basis(g,cutoff)
    f=expand(ambiguity_projection(g).subs(yy.X,-3125*xi))
    original=[Q(Poly(f,xi).nth(k)) for k in range(cutoff+1)]
    # B[0] has constant one; B[k] starts at xi^k.
    rem=list(original)
    a=[]
    for k in range(cutoff+1):
        a.append(rem[k]); rem=yy.add(rem,yy.scale(B[k],-rem[k]))
    assert rem==yy.zero(cutoff+1)
    return a

def yy_hol_data(g, cutoff):
    # yy.hol_yy includes derivatives; take higher order to avoid truncation artefacts.
    _, gens, N = yy.hol_yy(g,cutoff+3)
    n=cutoff+1
    return [val[:n] for val in gens], N[:n]

def jet_yy(P, g, cutoff):
    q,E,qi=source_series(g,cutoff)
    gens,N=yy_hol_data(g,cutoff)
    val=yy.eval_poly(P,gens)
    return yy.compose(yy.mul(N,val),qi)

def gv_difference(g, F, d):
    """Variation in n_g(d) for fixed lower genera and free energy delta F."""
    return sum(Q(int(mobius(r))*r**(2*g-3))*F[d//r] for r in divisors(d))

def extract_gv(maxdeg=6):
    F={}
    for g, poly in ((2,yy.P2),(3,yy.P3)):
        F[g] = jet_yy(poly,g,maxdeg)
    g2={}
    for d in range(1,maxdeg+1):
        lower = [Q(0)]+[
            F[2][n]-sum(Q(int(k)*N0_gv[n//int(k)],240) for k in divisors(n))
            for n in range(1,maxdeg+1)]
        g2[d]=gv_difference(2,lower,d)
    g3={}
    for d in range(1,maxdeg+1):
        lower=[Q(0)]+[
            F[3][n]-sum(Q(int(k)**3)*(Q(N0_gv[n//int(k)],6048)-Q(g2[n//int(k)],12))
                        for k in divisors(n))
            for n in range(1,maxdeg+1)]
        g3[d]=gv_difference(3,lower,d)
    for g,dat in ((2,g2),(3,g3)):
        assert dat=={k:Q(v) for k,v in bootstrap_table[g].items()}, (g,dat)
    return F,{2:g2,3:g3}

def pullback(g, cutoff, terms):
    """For YY monomials M, return exact rows d=1..cutoff of delta n_g[d].
    Lower-genus data are held fixed. We compute the *homogeneous* response,
    so the output is just Möbius sieve of the F_g shift.
    """
    columns = [jet_yy(M,g,cutoff) for M in terms]
    return [[gv_difference(g, F, d) for F in columns] for d in range(1,cutoff+1)]

def h_monomials():
    # centerings expressed in original v symbols
    return [yy.v1, yy.h2, yy.h3,
            yy.X*yy.v1, yy.X*yy.h2, yy.X*yy.h3,
            yy.X**2*yy.v1, yy.X**2*yy.h2, yy.X**2*yy.h3]

def reduction_polynomial(M, cutoff):
    """Unique X polynomial matching holomorphic YY M through xi^cutoff."""
    gens,_ = yy_hol_data(2,cutoff)
    jet=yy.eval_poly(M,gens)
    return poly_from_series(jet)

def make_report():
    report=[]
    def emit(msg=''):
        print(msg)
        report.append(str(msg))
    emit('EXACT BOOTSTRAP-ADAPTED BASIS FOR QUINTIC YY P2, P3')
    emit('D4-D2-D0 modular inputs: Table 2 of arXiv:2301.08066 (not computed here).')
    emit('E_g=omega0^(2g-2)(1-3125z)^(g-1), xi=-X/3125')
    emit('F_g^A = (-1/625)^(g-1) E_g P_g^hol')
    F, GV=extract_gv()
    for g in (2,3):
        cutoff=3*g-3
        emit('\n--- GENUS %d ---'%g)
        B=gv_basis(g,cutoff)
        for k,b in enumerate(B):
            emit('B[%d,%d] = %s' %(g,k,expand(poly_from_series(b))))
        a=dual_coordinates(g)
        emit('X-only P_%d projection = %s' %(g,ambiguity_projection(g)))
        for d,c in enumerate(a):
            emit('X-only coefficient in B[%d,%d] = %s; physical delta GV coefficient = %s' %(g,d,c,N0(g)*c))
    emit('\nFull genus 2,3 Gopakumar-Vafa data, through degree six:')
    for d in range(1,7):
        emit('d=%d, n2=%s, n3=%s' %(d,GV[2][d],GV[3][d]))
    emit('\nMIXED X^k H_i MONOMIALS (bootstrap pullback)')
    mixed=h_monomials()
    emit('Columns: '+ ', '.join(map(str,mixed)))
    for g in (2,3):
        M=pullback(g,6,mixed)
        emit('genus %d: entries are n_g shift from coefficient 1 in monomial:'%g)
        for d,row in enumerate(M,1):
            emit('d=%d: %s'%(d,', '.join(map(str,row))))
    emit('\nHOLOMORPHIC-JET-NULL MIXED COMBINATIONS (through degree 3g-3)')
    for g in (2,3):
        m=3*g-3
        gens,_=yy_hol_data(g,m)
        for label,mon in [('H1', yy.v1), ('H2', yy.h2),('H3',yy.h3),('XH1', yy.X*yy.v1),('XH2',yy.X*yy.h2)]:
            jet=yy.eval_poly(mon,gens)
            correction=poly_from_series(jet)
            # polynomial in xi; use substitution xi=-X/3125 in final comparison
            pure=expand(correction.subs(xi,-yy.X/3125))
            assert yy.eval_poly(expand(mon-pure),gens)==yy.zero(m+1)
            emit('g=%d, %s: YY combination %s minus {%s} has F_g q^1..q^%d zero' %
                 (g,label,mon,correction,m))
    Path('/mnt/data/yy_modular_bootstrap_results.txt').write_text('\n'.join(report)+'\n')
    emit('\nGIVEN PURE X AMBIGUITY COEFFICIENTS, PT CONSTRAINTS AT DEGREE FOUR:')
    assert GV[3][4] == 8625
    assert GV[2][4] + 4*GV[3][4] == 569250
    emit('PT(4,-2) = n4^3 = 8625')
    emit('PT(4,-1) = n4^2 + 4*n4^3 = 569250')
    for g in (2,3):
        P=yy.P2 if g==2 else yy.P3
        f=ambiguity_projection(g)
        proj=jet_yy(f,g,6)
        d=4
        amb=gv_difference(g,proj,d)
        part=GV[g][d]-amb
        emit('genus %d, d4: fixed BCOV part %s, pure-X ambiguity %s, total %s'
             %(g,part,amb,GV[g][d]))
    emit('\nGENUS 2 DEGREE-FOUR MODULAR RESPONSE IN GV-DUAL COORDINATES:')
    B2=gv_basis(2,3)
    basis2=[expand(poly_from_series(b).subs(xi,-yy.X/3125)) for b in B2]
    delta_row=pullback(2,4,basis2)[3]
    emit('delta n4^2 = (-1/625)*('+ ' + '.join('%s*a_%d'%(v/N0(2),i) for i,v in enumerate(delta_row))+')')
    assert [v/N0(2) for v in delta_row] == [627948300000,-2203986464,-3708877,4180]
    emit('\nA NONTRIVIAL MIXED YY COMBINATION INVISIBLE THROUGH DEGREE 3:')
    theta=yy.h3+5*yy.h2-Q(14636,3125)*yy.X**2-Q(160657632,1953125)*yy.X**3
    emit('Theta = H3 + 5 H2 - (14636/3125) X^2 - (160657632/1953125) X^3')
    for g in (2,3):
        Fth=jet_yy(theta,g,6)
        assert Fth[:4]==yy.zero(4)
        emit('genus %d: [q^0..q^3] F(Theta) = 0; [q^4] F(Theta) = %s' %(g,Fth[4]))
    emit('\nAll identities, mirror-map inversions, PT/GV checks, and jet null checks PASSED.')

if __name__=='__main__':
    make_report()
