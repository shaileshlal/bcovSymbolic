#!/usr/bin/env python3
"""No-oracle finite-difference-in-front-distance predictor.
Fits a local quadratic pulse to three *previous* r observations.
The evaluation grid is frozen at r=10,11; missing exact reference remains unscored.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
from predict_terminal import read,write
BASE=Path(__file__).resolve().parent

def lagrange(points,x):
    ans=F(0)
    for i,(u,v) in enumerate(points):
        a=F(v)
        for j,(w,_) in enumerate(points):
            if i!=j:a*=F(x-w,u-w)
        ans+=a
    return ans

def predict(train,sites):
    history=defaultdict(list)
    for a in train:
        key=(int(a['residue']),int(a['offset']))
        history[key].append((int(a['r']),F(a['epsilon'])))
    out=[]
    for a in sites:
        s,ell,r=int(a['residue']),int(a['offset']),int(a['r'])
        key=s,ell;old=sorted(history[key])
        assert len(old)>=3 and all(x<r for x,y in old)
        p=lagrange(old[-3:],r)
        out.append(dict(residue=s,offset=ell,r=r,depth=int(a['depth']),pulse_prediction=str(p),source='latest_three_prior_r'))
    return out

def run():
    out=predict(read(BASE/'data/diagonal_pulse_calibration.csv'),read(BASE/'data/diagonal_holdout_sites.csv'))
    write(BASE/'results/diagonal_pulse_predictions.csv',out)
    print('Diagonal pulse predictions',len(out))
    return out
if __name__=='__main__':run()
