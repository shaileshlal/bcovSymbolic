"""Same exact empirical affine pair model, frozen on curve degrees <=40.
No validation import. Generates predictions for degrees 54..70, shells11..15.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,csv,json
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from linear_extension_pairs import rowsys,reduce
from recursive_defects import read,write,phase_fit
from pulse_profile import K
from deep_envelope import W

def predict(sites):
    period,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))
    prior=json.load(open(ROOT/'results/v10_fitted_pulses.json'))
    dat=[z for z in read(ROOT/'data/observed_calibration_Q20_Q40_shell6_15.csv') if K(int(z['degree'])%5,int(z['shell']))-int(z['degree'])//5<=5]
    f=rowsys(dat,phase,period,prior)
    pivots={};bad=[]
    for r in dat:
        d,t=int(r['degree']),int(r['shell']);co,b,_,_=f(d,t,r['pt'])
        co,b=reduce(co,b,pivots)
        if co:
            p=min(co);c=co[p];pivots[p]=({k:v/c for k,v in co.items()},b/c)
        elif b:bad.append((d,t,str(b)))
    assert not bad,bad
    out=[]
    for d,t in sites:
        r=K(d%5,t)-d//5
        if r>5:
            out.append(dict(degree=d,shell=t,prediction='',front_distance=r,status='outside_trained_profile'));continue
        co,b,_,_=f(d,t)
        co,b=reduce(co,b,pivots)
        if not co:
            val=-b*phase[d%5,(d//5)%period]*W(d%5,t,d//5)
            if val.denominator==1:out.append(dict(degree=d,shell=t,prediction=val.numerator,front_distance=r,status='identified'))
            else:out.append(dict(degree=d,shell=t,prediction='',front_distance=r,status='nonintegral'))
        else:out.append(dict(degree=d,shell=t,prediction='',front_distance=r,status='unidentified'))
    return out

def main():
    out=predict([(d,t) for d in range(54,71) for t in range(11,16)])
    write(ROOT/'results'/'prospective_54_70.csv',out)
    print('frozen predictions',sum(bool(z['prediction']) for z in out),'of',len(out))
if __name__=='__main__':main()
