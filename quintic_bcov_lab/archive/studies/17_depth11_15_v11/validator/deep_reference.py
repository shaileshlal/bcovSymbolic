#!/usr/bin/env python3
"""Separate conditional exact PT reference for deep shells t=6..10.

Uses the exact rank-zero DT ansatz (9 polar equations; N=14) and the applicable
wall-crossing formula. For Q'=1, use independently known degree-one PT(q,Q)
=2875*q/(1+q)^2; Q'>1 candidates are NOT silently discarded, and those
sites are marked incomplete and excluded from scoring.
"""
from __future__ import annotations
import csv,sys,hashlib
from pathlib import Path
from fractions import Fraction as F
from pt_from_dt_solver import (applicable,wall_candidates,pt_mC,h_component,qhat,
   omega_rankzero,parity_sign, simple_vacuum_pt)
ROOT=Path(__file__).resolve().parent.parent
DELTA={0:0,1:1,2:2,3:2,4:1}

def pt_degree_one(m):
    # 2875 q/(1+q)^2, m>=1 (m=0: 0)
    return 0 if m<=0 else (-1)**(m-1)*2875*m

def evaluate(Q,t):
    m=pt_mC(Q)+DELTA[Q%5]+t
    if not applicable(Q,m):
        return dict(degree=Q,shell=t,m=m,status='outside_applicability',pt='',max_wall_degree='',wall_candidates='')
    cs=wall_candidates(Q,m)
    max_qp=max((x for x,_ in cs),default=0)
    if max_qp>1:
        return dict(degree=Q,shell=t,m=m,status='missing_lower_PT',pt='',max_wall_degree=max_qp,wall_candidates=len(cs))
    try:
        val=simple_vacuum_pt(Q,m)
        for qp,mp in cs:
            assert qp==1
            chi=m-mp+Q+qp-5
            Om=omega_rankzero(Q-qp,qhat(Q,m,qp,mp))
            val+=parity_sign(chi)*chi*pt_degree_one(mp)*Om
    except (ValueError,AssertionError) as e:
        return dict(degree=Q,shell=t,m=m,status='DT_order_missing:'+str(e),pt='',max_wall_degree=max_qp,wall_candidates=len(cs))
    return dict(degree=Q,shell=t,m=m,status='certified_within_conditional_framework',pt=str(val),max_wall_degree=max_qp,wall_candidates=len(cs))

def main():
    rows=[evaluate(d,t) for d in range(20,54) for t in range(6,11)]
    out=ROOT/'validator/reference_Q20_Q53.csv'
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    from collections import Counter
    print('Reference statuses overall:',Counter(r['status'] for r in rows))
    for lo,hi in [(20,30),(31,35),(36,40),(41,45),(46,50),(51,53)]:
        x=[r for r in rows if lo<=r['degree']<=hi]
        print(lo,hi,Counter(r['status'] for r in x))
    assert all(r['status']=='certified_within_conditional_framework' for r in rows if r['degree']>=41 and r['shell']<=7), 'needed safe cells missing'
if __name__=='__main__':main()
