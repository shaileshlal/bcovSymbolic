#!/usr/bin/env python3
"""Empirical deep-shell predictor: signed quadratic envelope and paired-residue scale.

DO NOT import validator code or future reference. The calibration table is a
separately supplied observed PT dataset; no wall-crossing/modularity formulas
are used. Extrapolation is in curve degree; confidence is NOT certified.
"""
import csv,sys
from fractions import Fraction as F
from collections import defaultdict,Counter
from pathlib import Path
from recursive_defects import phase_fit,read,write
BASE=Path(__file__).resolve().parent
C={0:10,1:6,2:4,3:4,4:6} # learned from three-point shallow quadratic fits
PEER={1:4,4:1,2:3,3:2}

def W(s,t,k):return 5*k*k+(2*s-5)*k+C[s]-2*t

def build_amplitudes(rows,phase,period):
    A=defaultdict(list)
    for r in rows:
        if r.get('pt','') in ('',None):continue
        d=int(r['degree']);t=int(r['shell']);s=d%5;k=d//5
        w=W(s,t,k)
        if w==0:continue
        a=F(int(r['pt']),phase[s,k%period]*w)
        A[(s,t)].append((k,a))
    return {st:sorted(vals) for st,vals in A.items()}

def select_amplitude(A,s,t,mode):
    own=A.get((s,t),[])
    if mode=='own_last':
        return (own[-1][1],'own_last',len(own)) if own else (None,'no_anchor',0)
    if mode=='paired_last':
        pool=[(k,a,ss) for ss in [s,PEER.get(s)] if ss is not None for k,a in A.get((ss,t),[])]
        if not pool:return (None,'no_anchor',0)
        pool.sort(key=lambda q:(q[0],q[2]==s))
        return pool[-1][1],'paired_last',len(pool)
    if mode=='paired_repeated':
        ss=[s,PEER.get(s)] if s in PEER else [s]
        pool=[(k,a,sub) for sub in ss for k,a in A.get((sub,t),[])]
        if not pool:return (None,'no_anchor',0)
        freq=Counter(a for k,a,sub in pool)
        # Repetition of EXACT scales across different degrees is evidence for
        # a stabilized polynomial, unlike mean/float-based smoothing.
        multi=[a for a,n in freq.items() if n>=2]
        if multi:
            def quality(a):return (freq[a],max(k for k,z,sub in pool if z==a),-abs(a.numerator))
            best=max(multi,key=quality)
            return best,'paired_repeated',len(pool)
        pool.sort(key=lambda q:(q[0],q[2]==s))
        return pool[-1][1],'paired_fallback_last',len(pool)
    raise ValueError(mode)

def forecast(history,sites,mode,phase,period):
    A=build_amplitudes(history,phase,period)
    rows=[]
    for z in sites:
        d=int(z['degree']);s=d%5;t=int(z['shell']);k=d//5
        scale,how,n=select_amplitude(A,s,t,mode)
        if scale is None:pred='';integer=False
        else:
            v=F(phase[s,k%period]*W(s,t,k))*scale
            # Rounding is reported explicitly; exact rational predicted value
            # is also provided, and 'exact_from_scale' records integrality.
            if v.denominator==1: pred=str(v.numerator);integer=True
            elif v>=0: pred=str((2*v.numerator+v.denominator)//(2*v.denominator));integer=False
            else:
                p=-v;pred=str(-((2*p.numerator+p.denominator)//(2*p.denominator)));integer=False
        rows.append(dict(degree=d,shell=t,residue=s,k=k,prediction=pred,
         scale=str(scale) if scale is not None else '',method=how,number_of_anchors=n,
         raw_prediction_integral=int(integer) if scale is not None else '',
         denominator=str(v.denominator) if scale is not None else ''))
    return rows

def main():
    root=BASE.parent
    cal=read(root/'data/shallow_calibration_Q1_Q20.csv')
    period,phase=phase_fit(cal)
    train=read(root/'data/calibration_deep_Q20_Q40.csv')
    assert all(int(r['degree'])<=40 for r in train)
    sites=read(root/'data/forecast_sites_Q41_Q53.csv')
    for mode in ('own_last','paired_last','paired_repeated'):
        out=forecast(train,sites,mode,phase,period)
        write(root/'results'/f'predictions_{mode}.csv',out)
        print(mode,'predicted',sum(bool(r['prediction']) for r in out),'of',len(out))
if __name__=='__main__':main()
