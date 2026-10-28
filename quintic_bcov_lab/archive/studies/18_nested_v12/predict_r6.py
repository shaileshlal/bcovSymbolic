#!/usr/bin/env python3
"""Conditional first-failure r=6 shape learner.
Input is Z(t) normalized using an externally *observed mature PT scale* B_group,t.
Predicts Z, not the unknown scale itself.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from predict_terminal import read,write
from predict_diagonal import lagrange
BASE=Path(__file__).resolve().parent

def forecast(training,sites):
    hist=defaultdict(list)
    for a in training:hist[int(a['residue'])].append(a)
    out=[]
    for a in sites:
        s=int(a['residue']);t=int(a['t']);w=int(a['W']);old=sorted(hist[s],key=lambda q:int(q['t']))
        assert len(old)>=3 and all(int(q['degree'])<=53 for q in old)
        early=[(int(q['W']),F(q['Z'])) for q in old[:2]]
        late=[(int(q['t']),F(q['Z'])) for q in old[-3:]]
        late2=[(int(q['W']),F(q['Z'])) for q in old[-2:]]
        p0=lagrange(early,w)
        p1=lagrange(late,t)
        p2=lagrange(late2,w)
        out.append(dict(residue=s,t=t,degree=int(a['degree']),W=w,
            earliest_affine_Z=str(p0),latest_quadratic_Z=str(p1),latest_affine_Z=str(p2),integer_correction=str(p1-p0)))
    return out

def run():
    out=forecast(read(BASE/'data/r6_Z_calibration.csv'),read(BASE/'data/r6_Z_holdout_sites.csv'))
    write(BASE/'results/r6_Z_predictions.csv',out)
    print('r6 corrected normalized predictions',len(out))
    return out
if __name__=='__main__':run()
