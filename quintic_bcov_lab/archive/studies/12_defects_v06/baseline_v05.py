#!/usr/bin/env python3
"""Empirical, non-modular quintic near-extremal PT/GV experiments.

Training labels: Q <= 20, whole-degree split. Holdout: Q 21..30.
NO modular forms, wall-crossing equations, rank-zero DT series or BCOV oracle
are loaded in these predictors. Only integers in calibration.csv are used.

Usage:
    python run_experiments.py --train calibration.csv --sites sites_unlabeled.csv --out predictions.csv
"""
from __future__ import annotations
import csv, json, hashlib, argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path

DELTA0={0:0,1:1,2:2,3:2,4:1}

def read_csv(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f))

def exact_phase_fit(training):
    """Search shortest noncontradictory sign period 1..12, all five residues."""
    pt={(int(r['degree']),int(r['delta'])):int(r['pt']) for r in training if r['pt'] not in ('','None')}
    for period in range(1,13):
        patterns={}
        conflict=False
        for d in sorted({int(r['degree']) for r in training}):
            s=d%5; k=d//5; val=pt.get((d,DELTA0[s]))
            if val is None or val==0: continue
            key=(s,k%period)
            sg=1 if val>0 else -1
            if key in patterns and patterns[key]!=sg: conflict=True;break
            patterns[key]=sg
        if not conflict:
            if all((s,j) in patterns for s in range(5) for j in range(period)):
                return period,patterns
    raise ValueError('Sign period not uniquely identifiable from available training data')

def lagrange_quad(xs,ys):
    """Lagrange interpolation in fractions (any 3 training points)."""
    assert len(xs)==len(ys)==3
    def f(z):
        total=Fraction(0)
        for i in range(3):
            f_i=Fraction(ys[i])
            for j in range(3):
                if i!=j:f_i*=Fraction(z-xs[j],xs[i]-xs[j])
            total+=f_i
        return total
    return f

def nearest_int(q):
    q=Fraction(q)
    if q>=0:return (2*q.numerator+q.denominator)//(2*q.denominator)
    return -nearest_int(-q)

def train_predict(rows,site_rows,criterion=Fraction(1,10000)):
    # Period is selected using Q1..20, labels PT at FIRST NONZERO shell only.
    period,signs=exact_phase_fit(rows)
    tab={(int(r['degree']),int(r['delta'])):int(r['pt']) for r in rows if r['pt'] not in ('','None')}
    available=sorted({int(r['degree']) for r in rows if 6<=int(r['degree'])<=20})
    train_degrees={s:[q for q in available if q%5==s] for s in range(5)}
    assert all(len(v)==3 for v in train_degrees.values())
    def phase(d):return signs[(d%5,(d//5)%period)]
    def y(d,t):return phase(d)*tab[(d,DELTA0[d%5]+t)]
    models={}
    lead_models={}
    for s,qs in train_degrees.items():
        lead_models[s]=lagrange_quad([q//5 for q in qs],[y(q,0) for q in qs])
    # Make preliminary training-only affine anchor dictionary, no ref labels.
    affine={}
    for s,qs in train_degrees.items():
        base=[y(q,0) for q in qs]
        for t in range(1,6-DELTA0[s]):
            v=[y(q,t) for q in qs]
            if base[1]==base[0]:continue
            A=Fraction(v[1]-v[0],base[1]-base[0]);C=Fraction(v[0])-A*base[0]
            if all(A*u+C==v_i for u,v_i in zip(base,v)):
                affine[(s,t)]=(A,C)
    # Highest-residue anchor on reflection pair, only when exact relation
    # exists on anchor train and latest own training mismatch <1e-4.
    for s,qs in train_degrees.items():
        for t in range(6-DELTA0[s]):
            vals=[y(q,t) for q in qs]
            fallback=lagrange_quad([q//5 for q in qs], vals)
            model=dict(type='quadratic',coeffs=[str(fallback(i)) for i in (0,1,2)])
            if s in (1,2) and t>0:
                mirrored=5-s
                anchor=affine.get((mirrored,t))
                if anchor:
                    A,C=anchor
                    q_last=qs[-1]
                    y0=y(q_last,0)
                    err=abs(A*y0+C-y(q_last,t))
                    if err<=criterion*abs(y(q_last,t)):
                        model=dict(type='mirror_affine',anchor=mirrored,A=str(A),C=str(C),latest_training_relative_error=str(Fraction(err,abs(y(q_last,t)))))
            models[(s,t)]=(model,fallback)
    prediction=[]
    for rr in site_rows:
        d=int(rr['degree']); t=int(rr['shell']);s=d%5
        k=d//5
        assert 21<=d<=30
        model,fallback=models[(s,t)]
        if model['type']=='quadratic':
            unphased=fallback(k)
        else:
            unphased=Fraction(model['A'])*lead_models[s](k)+Fraction(model['C'])
        val=unphased*phase(d)
        prediction.append(dict(degree=d,residue=s,k=k,shell=t,delta=DELTA0[s]+t,
                               predicted_PT=str(nearest_int(val)), exact_rational_prediction=str(val),
                               model=model['type'],phase=phase(d)))
    return prediction,dict(period=period,signs={f'{k[0]},{k[1]}':v for k,v in signs.items()},affine={f'{k[0]},{k[1]}':(str(a),str(c)) for k,(a,c) in affine.items()},models={f'{k[0]},{k[1]}':a for k,(a,b) in models.items()})

def main():
    p=argparse.ArgumentParser();p.add_argument('--train',default='data/calibration.csv');p.add_argument('--sites',default='data/sites_unlabeled.csv');p.add_argument('--out',default='results/pt_predictions.csv');args=p.parse_args()
    rows=read_csv(args.train);sites=read_csv(args.sites)
    results,info=train_predict(rows,sites)
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,results[0]);w.writeheader();w.writerows(results)
    (out.parent/'model_audit.json').write_text(json.dumps(info,indent=2))
    (out.parent/'input_hashes.json').write_text(json.dumps({p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (args.train,args.sites)},indent=2))
    print(f'Predicted {len(results)} PT integers; selected sign period {info["period"]}; no holdout reference file opened.')
    print('Transfer models:',[(k,v['type']) for k,v in info['models'].items() if v['type']!='quadratic'])

if __name__=='__main__':main()
