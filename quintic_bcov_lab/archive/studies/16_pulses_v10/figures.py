#!/usr/bin/env python3
"""Optional figures; matplotlib only, no external network."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def read(rel):
    with (ROOT/rel).open(newline='') as f:return list(csv.DictReader(f))

def main():
    import matplotlib.pyplot as plt
    old=read('validator/previous_v09_gated_predictions.csv')
    new=read('results/pulse_predictions_Q41_Q53.csv')
    t=list(range(6,11))
    a=[sum(int(z['shell'])==i and bool(z['prediction']) for z in old) for i in t]
    b=[sum(int(z['shell'])==i and bool(z['prediction']) for z in new) for i in t]
    fig,ax=plt.subplots(figsize=(7.2,3.8))
    ax.bar([i-0.18 for i in t],a,width=0.36,label='Frozen v0.9')
    ax.bar([i+0.18 for i in t],b,width=0.36,label='Pulse profiles v1.0')
    ax.set(xticks=t,xlabel='Shell depth t',ylabel='Exact predictions (out of 13)',ylim=(0,14))
    ax.legend();ax.grid(axis='y',alpha=0.2);fig.tight_layout()
    fig.savefig(ROOT/'figures/exact_coverage_by_depth.png',dpi=160);plt.close(fig)
    p=read('validator/front_transport_pairs.csv')
    on=[];off=[]
    from fractions import Fraction as F
    import math
    for z in p:
        v=F(z['pulse_before']);v2=F(z['pulse_after']);r=abs(v+v2)/max(abs(v),abs(v2))
        if z['at_front']=='1':on.append(0.)
        else:off.append(float(r))
    fig,ax=plt.subplots(figsize=(7,3.5))
    ax.scatter([0]*len(on),on,label='Terminal front')
    ax.scatter([1]*len(off),off,label='Still active')
    ax.set(xticks=[0,1],xticklabels=['Last active pulse','Interior active pulse'],ylabel='|D + diagonal D| / max |D|',title='Terminal pulse cancellation')
    ax.grid(alpha=0.15);fig.tight_layout()
    fig.savefig(ROOT/'figures/pulse_transport_residual.png',dpi=160);plt.close(fig)
    print('Saved figures in figures/')

if __name__=='__main__':main()
