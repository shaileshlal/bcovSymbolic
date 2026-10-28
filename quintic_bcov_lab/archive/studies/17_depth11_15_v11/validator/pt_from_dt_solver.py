"""Exact single-/multi-wall DT->PT solver in the quintic.

The DT input is reconstructed from 7 polar coefficients in modular_core.py.
The wall inequalities and signs follow Alexandrov et al. arXiv:2301.08066v2
Eq. (4.12)--(4.15). The BMT/wall-crossing hypotheses are inherited.

Only degree-one PT(1,m), m=0,1,2, is required for Q=21, delta<=5.
No degree-21 PT or GV reference coefficients are read by this module.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import floor, ceil
from modular_core import h_series

KAPPA=5
CHI_DIV=5
CHI_DIV_HOL=F(55,1)

def parity_sign(a:int)->int:return 1 if a%2==0 else -1

def pt_mC(Q:int)->int:
    return -floor(F(Q*Q+5*Q,10))

def h_component(mu):
    return min(mu%5, (-mu)%5)

def omega_rankzero(mu:int,qhat:F)->int:
    mu=h_component(mu)
    base,coeffs=h_series(mu)
    # h_mu = sum c_n q^(base+n) = sum Omega(qhat) q^(-qhat).
    ix=-qhat-base
    assert ix.denominator==1,(mu,qhat,base,ix)
    idx=int(ix)
    if idx<0 or idx>=len(coeffs):
        raise ValueError(f'Rank-zero DT power q^{idx} lies outside exact q-series range')
    return coeffs[idx]

def qhat(Q,m,Qp=0,mp=0):
    return F(mp-m,1)-F((Qp-Q)**2,2*KAPPA)-F(Q+Qp,2)+F(55,24)

def simple_vacuum_pt(Q:int,m:int)->int:
    chi=m+Q-CHI_DIV
    return parity_sign(chi)*chi*omega_rankzero(Q,qhat(Q,m))

def wall_candidates(Q,m):
    """All potentially contributing lower-degree PT sites, exact rational bounds.

    The enumeration is useful only where alpha>f(Q/5) so the theorem applies.
    Candidates with PT=0 still appear and are not dropped before evaluation.
    """
    alpha=F(-3*m,2*Q)
    maxQp=floor(Q+KAPPA*(F(1,2)-alpha))
    rows=[]
    for Qp in range(1,maxQp+1):
        low=-F(Qp*Qp+5*Qp,10)
        high=F(m)+F((Q-Qp)**2,10)+F(Q+Qp,2)
        for mp in range(ceil(low),floor(high)+1):
            rows.append((Qp,mp))
    return rows

def applicable(Q,m):
    """Exact condition alpha>piecewise f from the paper's Eq (4.10)."""
    x=F(Q,5)
    alpha=F(-3*m,2*Q)
    if x<1: return alpha>x+F(1,2)
    if x<F(15,8):return alpha>0 and alpha*alpha>2*x+F(1,4)
    if x<F(9,4):return alpha>F(2,3)*x+F(3,4)
    if x<3:return alpha>x/F(3)+F(3,2)
    return alpha>x/F(2)+1

PT_SMALL={
    (1,0):0, (1,1):2875,(1,2):-5750,
}

def pt_from_rankzero_dt(Q,m):
    if not applicable(Q,m):
        raise ValueError(f'Wall-crossing relation not certified for Q={Q}, m={m}')
    base=simple_vacuum_pt(Q,m)
    contributions=[]
    for Qp,mp in wall_candidates(Q,m):
        if (Qp,mp) not in PT_SMALL:
            raise KeyError('Need independent lower-degree PT input',Qp,mp)
        chip=m-mp+Q+Qp-CHI_DIV
        Om=omega_rankzero(Q-Qp,qhat(Q,m,Qp,mp))
        val=parity_sign(chip)*chip*PT_SMALL[Qp,mp]*Om
        contributions.append(dict(degree=Qp,chi=mp,pt_lower=PT_SMALL[Qp,mp],
                                  index=chip,rankzero_dt=Om,contribution=val))
    return base+sum(r['contribution'] for r in contributions),base,contributions

if __name__=='__main__':
    EXPECT={-53:319125,-52:40986000,-51:132211590000,
            -50:-105493513413000,-49:19064998024136500}
    for m,t in EXPECT.items():
        val,base,terms=pt_from_rankzero_dt(21,m)
        print(f'PT(21,{m}) predicted {val}; vacuum {base}; lower walls', terms)
        assert val==t,(m,val,t)
    print('Five exact DT->PT held-out matches, without PT(21,m) targets in fit')
