#!/usr/bin/env python3
"""Transport terminal pulse forecasts to the next normalized Z coefficient.
Requires three previous observed normalized Z values, but no target Z.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
from predict_terminal import read,write
BASE=Path(__file__).resolve().parent

def predict(pulses,anchors):
    idx=defaultdict(dict)
    for a in anchors:idx[int(a['residue']),int(a['r']),int(a['target_t'])][int(a['anchor_t'])]=F(a['Z'])
    out=[]
    for p in pulses:
        s,r,t=int(p['residue']),int(p['r']),int(p['t'])
        v=idx[s,r,t]
        assert all(q in v for q in (t-1,t-2,t-3))
        z=3*v[t-1]-3*v[t-2]+v[t-3]+F(p['terminal_pulse'])
        out.append(dict(residue=s,r=r,depth=t,predicted_Z=str(z)))
    return out

def run():
    p=read(BASE/'results/terminal_pulse_predictions.csv')
    anchors=read(BASE/'data/terminal_pulse_reconstruction_anchors.csv')
    out=predict(p,anchors)
    write(BASE/'results/terminal_Z_reconstruction.csv',out)
    print('Terminal Z predictions',len(out))
    return out
if __name__=='__main__':run()
