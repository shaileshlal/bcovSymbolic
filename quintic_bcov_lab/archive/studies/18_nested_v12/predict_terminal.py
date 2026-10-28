#!/usr/bin/env python3
"""Pure empirical integer terminal-pulse predictor.

No modular forms, DT, wall crossing or future PT reference is used.
The calibration labels are supplied exact integers of normalized PT-pulse signals.
Reported predictive gains are conditional on that normalization (B_t scale).
"""
from __future__ import annotations
import csv,sys
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
BASE=Path(__file__).resolve().parent
DATA=BASE/'data'

def read(path):
    return list(csv.DictReader(open(path,newline='')))

def terminal_t(s,r):return 2*r-1 if s==1 else 2*r+3

def second_difference(vals):
    return vals[2]-2*vals[1]+vals[0]

def build_model(training):
    d=defaultdict(dict)
    for row in training:
        s=int(row['residue']);r=int(row['r']);t=int(row['t']);tau=int(row['terminal_pulse'])
        assert t==terminal_t(s,r)
        assert r not in d[s]
        d[s][r]=tau
    for s,n in {0:3,1:3,2:3,3:2,4:2}.items():
        assert len(d[s])==n,(s,len(d[s]),n)
    curv={s:second_difference([d[s][r] for r in sorted(d[s])]) for s in (0,1,2)}
    # Cross-family curvature relations, empirically selected during model discovery.
    curv[3]=-curv[1]
    curv[4]=curv[0]
    return d,curv

def predict(training,sites):
    d,curv=build_model(training)
    out=[]
    for site in sites:
        s=int(site['residue']);r=int(site['r']);t=int(site['t']);assert t==terminal_t(s,r)
        rs=sorted(d[s]);a,b=rs[:2]
        assert b==a+1 and r>=a
        n=r-a
        v=d[s][a]+n*(d[s][b]-d[s][a])+curv[s]*n*(n-1)//2
        assert (curv[s]*n*(n-1))%2==0
        out.append(dict(residue=s,r=r,t=t,terminal_pulse=str(v),source='quadratic_r_with_paired_curvature'))
    return out,curv

def write(path,rows):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def run():
    training=read(DATA/'terminal_pulse_calibration.csv')
    sites=read(DATA/'terminal_pulse_holdout_sites.csv')
    preds,curv=predict(training,sites)
    BASE.joinpath('results').mkdir(exist_ok=True)
    write(BASE/'results/terminal_pulse_predictions.csv',preds)
    print('Training labels:',len(training),'Predictions:',len(preds),'Curvatures:',curv)
    return preds
if __name__=='__main__':run()
