#!/usr/bin/env python3
"""Autonomous, entirely empirical roll-forward, degrees 31..53.

Inputs: the FROZEN v0.7 model and archived degree<=30 PT observations.
Neither reference labels nor DT/modular/wall-crossing functions are imported.
Uses blocks of five degrees, appending its own predicted values only.
"""
import argparse
from pathlib import Path
from recursive_defects import read,write,phase_fit,predict

ROOT=Path(__file__).resolve().parent
DEPTH={0:6,1:5,2:4,3:4,4:5}

def sites(lo,hi):
    return [{'degree':str(d),'shell':str(t)} for d in range(lo,hi+1) for t in range(DEPTH[d%5])]

def run(stop=53):
    if not 31<=stop<=100:raise ValueError('stop must be 31..100')
    cal=read(ROOT/'data/calibration_Q1_Q20.csv')
    train=read(ROOT/'data/observed_Q6_Q30.csv')
    period,phase=phase_fit(cal)
    if period!=4:raise AssertionError(f'Archive phase was {period}, not four')
    outputs=[]
    for start in range(31,stop+1,5):
        end=min(start+4,stop)
        before=len(train)
        preds=predict(train,sites(start,end),period,phase)
        outputs.extend(preds)
        for z in preds:
            # NEVER attach true labels: new training rows contain predictions only.
            train.append({'degree':str(z['degree']), 'shell':str(z['shell']), 'pt':str(z['mirror_PT'])})
        assert len(train)-before==len(preds)
    return outputs

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--stop',type=int,default=53)
    p.add_argument('--out',default=str(ROOT.parent/'results/autonomous_predictions_Q31_Q53.csv'))
    opt=p.parse_args()
    pred=run(opt.stop)
    write(opt.out,pred)
    print('Predicted',len(pred),'PT integers to degree',opt.stop,'without reference labels')
