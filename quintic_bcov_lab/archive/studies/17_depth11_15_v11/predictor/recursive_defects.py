#!/usr/bin/env python3
"""Purely empirical, exact-integer PT finite-difference and shell-residual models.

Train only: data/calibration_Q1_Q20.csv and data/update_Q21_Q25.csv.
Forecast sites: data/sites_Q26_Q30.csv or user-supplied degree sites.
This module NEVER reads a heldout reference. Python standard library only.
"""
from __future__ import annotations
import argparse,csv,math,json
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PEER={1:4,4:1,2:3,3:2}

def read(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def write(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:raise ValueError('empty table')
    with path.open('w',newline='') as f:
        fieldnames=list(dict.fromkeys(k for row in rows for k in row))
        w=csv.DictWriter(f,fieldnames);w.writeheader();w.writerows(rows)

def phase_fit(rows):
    # Uses ONLY first-shell signs from calibration and searches periods in increasing order.
    lead={(int(z['degree']),int(z['shell'])):int(z['pt']) for z in rows if z.get('pt') not in ('',None)}
    d0={0:0,1:1,2:2,3:2,4:1}
    for p in range(1,13):
        pat={};valid=True
        for d in sorted({k[0] for k in lead}):
            v=lead.get((d,0))
            if v is None or v==0:continue
            key=(d%5,(d//5)%p);sg=1 if v>0 else -1
            if key in pat and pat[key]!=sg:valid=False;break
            pat[key]=sg
        if valid and all((s,r) in pat for s in range(5) for r in range(p)):
            return p,pat
    raise RuntimeError('sign period not identifiable')

def groups(rows,phase,period):
    d=defaultdict(list)
    for z in rows:
        if z.get('pt') in ('',None):continue
        deg=int(z['degree']); s=deg%5;t=int(z['shell']);k=deg//5
        d[s,t].append((k,phase[s,k%period]*int(z['pt'])))
    return {key:sorted(x) for key,x in d.items()}

def curvature(seq):
    a,b,c=[v for _,v in seq[-3:]];return c-2*b+a

def previous_curvature(seq):
    return curvature(seq[-4:-1]) if len(seq)>=4 else None

def defect(seq):
    return curvature(seq)-previous_curvature(seq) if len(seq)>=4 else None

def gauss(A,b):
    n=len(A)
    a=[[Q(x) for x in row]+[Q(bb)] for row,bb in zip(A,b)]
    for j in range(n):
        i=next((i for i in range(j,n) if a[i][j]),None)
        if i is None:return None
        a[j],a[i]=a[i],a[j]
        fac=a[j][j]
        for z in range(j,n+1):a[j][z]/=fac
        for i in range(n):
            if i==j:continue
            r=a[i][j]
            for z in range(j,n+1):a[i][z]-=r*a[j][z]
    return [a[i][n] for i in range(n)]

def evalpoly(points,k):
    ans=Q(0)
    for i,(xi,yi) in enumerate(points):
        f=Q(yi)
        for j,(xj,_) in enumerate(points):
            if j!=i:f*=Q(k-xj,xi-xj)
        ans+=f
    return ans

def nearest_integer(x):
    x=Q(x)
    if x>=0:return (2*x.numerator+x.denominator)//(2*x.denominator)
    return -nearest_integer(-x)

def next_pred(g,s,t,k,phase,period):
    seq=g[s,t]
    if not seq or seq[-1][0]!=k-1 or len(seq)<3:raise ValueError((s,t,k,seq))
    cur=curvature(seq);prev=previous_curvature(seq)
    prev_delta=(cur-prev if prev is not None else None)
    base_val=seq[-1][1]*2-seq[-2][1]
    pred_quad=base_val+cur
    pred_cubic=base_val+cur+prev_delta if prev_delta is not None else None
    peer=PEER.get(s);pseq=g.get((peer,t))
    peer_c=curvature(pseq) if pseq and len(pseq)>=3 else None
    bracket=(peer_c is not None and prev is not None and prev_delta!=0 and min(prev,cur)<=peer_c<=max(prev,cur))
    mirror_curv=peer_c if bracket else cur
    pred_mirror=base_val+mirror_curv
    product=None
    if t>0:
        lead=g.get((s,0));
        if lead and len(lead)>=3:
            # Model T_t(k) = A k * T_0(k) + B T_0(k)+C,
            # fixed exactly from three previously observed points.
            recent=seq[-3:];lv=dict(lead)
            if all(ki in lv for ki,_ in recent):
                coeff=gauss([[ki*lv[ki],lv[ki],1] for ki,yi in recent],[yi for _,yi in recent])
                if coeff is not None:
                    L=evalpoly(lead[-3:],k)
                    product=coeff[0]*k*L+coeff[1]*L+coeff[2]
    # model intervals are empirical, NOT probabilistically certified:
    if prev_delta is None:lo=hi=None
    elif prev_delta==0:lo=hi=0
    else:lo=-abs(prev_delta);hi=abs(prev_delta)
    sign=phase[s,k%period]
    return dict(degree=5*k+s,residue=s,k=k,shell=t,sign=sign,
     quadratic_PT=int(sign*pred_quad), cubic_PT=int(sign*pred_cubic) if pred_cubic is not None else '',
     linear_times_quadratic_PT=str(nearest_integer(sign*product)) if product is not None else '',
     linear_times_quadratic_rational=str(sign*product) if product is not None else '',
     mirror_PT=int(sign*pred_mirror),mirror_used=int(bracket),
     observed_curvature=str(cur),previous_curvature=str(prev) if prev is not None else '',
     paired_curvature=str(peer_c) if peer_c is not None else '',
     previous_defect=str(prev_delta) if prev_delta is not None else '',
     defect_interval_low=str(lo) if lo is not None else '',
     defect_interval_high=str(hi) if hi is not None else '',
     pulse_state='stable_observed' if prev_delta==0 else 'active_observed' if prev_delta is not None else 'undetermined')

def predict(train,sites,period,pat):
    g=groups(train,pat,period);out=[]
    for z in sites:
        d=int(z['degree']);s=d%5;t=int(z['shell']);k=d//5
        out.append(next_pred(g,s,t,k,pat,period))
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--calibration',default=str(ROOT/'data/calibration_Q1_Q20.csv'))
    p.add_argument('--update',default=str(ROOT/'data/update_Q21_Q25.csv'))
    p.add_argument('--sites',default=str(ROOT/'data/sites_Q26_Q30.csv'))
    p.add_argument('--out',default=str(ROOT/'results/predictions_Q26_Q30.csv'))
    args=p.parse_args()
    cal=read(args.calibration);upd=read(args.update)
    period,pat=phase_fit(cal)
    out=predict(cal+upd,read(args.sites),period,pat)
    write(args.out,out)
    Path(args.out).with_suffix('.meta.json').write_text(json.dumps({'period':period,'signed_phases':{f'{s},{r}':v for (s,r),v in pat.items()},'inputs':[args.calibration,args.update,args.sites],'reference_not_read':True},indent=2))
    print('Generated',len(out),'predictions, sign period',period,'(no target reference opened)')
if __name__=='__main__':main()
