"""Frozen simple baseline families on exact integer calibration (d<=40)."""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import csv
import sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from recursive_defects import phase_fit,read,write
from deep_envelope import W

def lagrange(pts,x):
    out=F(0)
    for i,(k,y) in enumerate(pts):
        t=y
        for j,(m,z) in enumerate(pts):
            if i!=j:t*=F(x-m,k-m)
        out+=t
    return out

def run():
    period,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))
    train=read(ROOT/'data/observed_calibration_Q20_Q40_shell6_15.csv')
    grp=defaultdict(list)
    for r in train:
        d,t=int(r['degree']),int(r['shell']);s=d%5;k=d//5
        if t<11:continue
        sign=phase[s,k%period]
        grp[s,t].append((k,F(int(r['pt']),sign)))
    for z in grp.values():z.sort()
    sites=read(ROOT/'data/withheld_sites_Q41_Q53_shell11_15.csv')
    methods=('quadratic','linear','affine_envelope','constant_envelope','cubic')
    out=[]
    for r in sites:
        d,t=int(r['degree']),int(r['shell']);s=d%5;k=d//5
        a=grp[s,t];sgn=phase[s,k%period];w=W(s,t,k)
        item={'degree':d,'shell':t,'n_train':len(a)}
        for model in methods:
            raw=None
            if model=='quadratic' and len(a)>=3:raw=lagrange(a[-3:],k)
            if model=='cubic' and len(a)>=4:raw=lagrange(a[-4:],k)
            if model=='linear' and len(a)>=2:raw=lagrange(a[-2:],k)
            if model=='constant_envelope' and len(a):
                j,v=a[-1];raw=F(w,W(s,t,j))*v
            if model=='affine_envelope' and len(a)>=2:
                (j,v),(h,u)=a[-2:]
                # y=W*(P)+Q
                P=F(u-v,W(s,t,h)-W(s,t,j));Q=u-P*W(s,t,h)
                raw=P*w+Q
            if raw is None:item[model]='';item[model+'_exact_integer']=''
            else:
                raw*=sgn;item[model]=str(raw);item[model+'_exact_integer']=int(raw.denominator==1)
        out.append(item)
    dest=ROOT/'results';dest.mkdir(exist_ok=True)
    write(dest/'baseline_predictions.csv',out)
    print('rows',len(out),'quadratic_sites',sum(bool(r['quadratic']) for r in out))
    return out
if __name__=='__main__':run()
