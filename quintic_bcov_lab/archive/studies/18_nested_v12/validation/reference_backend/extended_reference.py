"""Independent, conditional DT/wall-crossing reference for PT shells t=11..15.

Reference code must not be imported by or exposed to the empirical learner.
Uses rank-zero DT q-series (separate modular-core validator) and small-degree
PT(1,m), PT(2,m), PT(3,m) derived from genus 0/1 GV data.
Refuses to fill a candidate involving degree >3 or missing DT order.
"""
from __future__ import annotations
import csv
from math import comb
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from pt_from_dt_solver import (applicable,wall_candidates,pt_mC,qhat,
    omega_rankzero,parity_sign, simple_vacuum_pt)

ROOT=Path(__file__).resolve().parent.parent
DELTA={0:0,1:1,2:2,3:2,4:1}
N01,N02,N13=2875,609250,609250

def pt_small(d,m):
    if d==1:return 2875*m*(-1)**(m-1) if m>=1 else 0
    if d==2:
        if m<1:return 0
        # [Q^2]Z_PT: n0,2*q/(1+q)^2+ (n0,1^2/2)*q^2/(1+q)^4
        #                            - (n0,1/2)*q^2/(1-q^2)^2.
        ans=F(N02*m*(-1)**(m-1))
        if m>=2:
            ans+=F(N01*N01,2)*(-1)**(m-2)*comb(m+1,3)
            if m%2==0:ans-=F(N01*m,4)
        assert ans.denominator==1,(d,m,ans)
        return ans.numerator
    if d==3:
        if m<0:return 0
        if m==0:return N13
        raise ValueError('PT(3,m>0) needs higher-degree GV data')
    raise ValueError('Unsupported lower-degree PT input')

def evaluate(d,t):
    m=pt_mC(d)+DELTA[d%5]+t
    if not applicable(d,m):
        return (d,t,m,'outside_applicability','','',0)
    cs=wall_candidates(d,m)
    max_qp=max([0]+[qp for qp,mp in cs])
    if max_qp>3:
        return (d,t,m,'missing_lower_PT','',max_qp,len(cs))
    try:
        val=simple_vacuum_pt(d,m)
        for qp,mp in cs:
            v=pt_small(qp,mp)
            if not v:continue
            chi=m-mp+d+qp-5
            DT=omega_rankzero(d-qp,qhat(d,m,qp,mp))
            val+=parity_sign(chi)*chi*v*DT
    except (ValueError,AssertionError) as e:
        return (d,t,m,'DT_order_missing:'+str(e),'',max_qp,len(cs))
    return (d,t,m,'certified_conditional',str(val),max_qp,len(cs))

def run():
    rows=[]
    for d in range(20,54):
        for t in range(6,16):
            if t<=10:continue
            rows.append(evaluate(d,t))
    out=ROOT/'validator/reference_deep_11_15.csv'
    with out.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(('degree','shell','m','status','pt','max_wall_degree','nwall'))
        writer.writerows(rows)
    for a,b in [(20,30),(31,40),(41,53)]:
        print(a,b,Counter(r[3] for r in rows if a<=r[0]<=b))
    print('Examples',[(d,t,pt) for d,t,m,s,pt,qp,n in rows if (d,t) in [(35,11),(41,11),(43,13),(53,15)]])
    return rows
if __name__=='__main__':run()
