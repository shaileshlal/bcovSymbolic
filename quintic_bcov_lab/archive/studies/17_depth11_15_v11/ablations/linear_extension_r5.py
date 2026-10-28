"""Exact, low-capacity affine-in-1/W extension; degree<=50 training only.
Attempts to transfer v1.0 shell residual profiles without interpreting them
through modularity or wall crossing. Profiles not learned before depth11 are
free shared symbols. Rejects inconsistent representations.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json,csv,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(HERE))
from recursive_defects import read,write,phase_fit
from deep_envelope import W
from pulse_profile import K, group

def rowsys(rows,phase,period,prior,fix_prior=True):
    # all unknown B[g,t], P[s,r],Q[s,r]; known priors from v1.0 from t<=10.
    constB={(int(r['class']),int(r['shell'])):F(r['rational']) for r in prior['B']}
    constPQ={(int(r['residue']),int(r['front_distance'])):(F(r['P']),F(r['Q'])) for r in prior['profiles']}
    def equation(d,t,pt=None):
        s=d%5;k=d//5;r=K(s,t)-k;w=W(s,t,k)
        if w==0:return None
        co={};offset=F(0)
        key=group(s),t
        if key in constB and fix_prior: offset+=constB[key]
        else: co[('B',key[0],t)]=F(1)
        if r>0:
            sign=-1 if t%2 else 1
            if (s,r) in constPQ and fix_prior:
                p,q=constPQ[s,r];offset+=sign*(p+q/F(w))
            else:
                co[('P',s,r)]=F(sign)
                co[('Q',s,r)]=F(sign,w)
        value=F(int(pt),phase[s,k%period]*w)-offset if pt is not None else -offset
        return co,value,offset,r
    return equation

def reduce(a,b,pivots):
    # canonicalized sparse elimination (variables lexical order); pivots rows may carry fractions
    for p,(pc,pb) in sorted(pivots.items()):
        v=a.get(p,F(0))
        if v:
            for k,vv in pc.items():
                a[k]=a.get(k,F(0))-v*vv
                if not a[k]:a.pop(k)
            b-=v*pb
    return a,b

def run():
    cal=read(ROOT/'data/shallow_calibration_Q1_Q20.csv')
    period,phase=phase_fit(cal)
    prior=json.load(open(ROOT/'results/v10_fitted_pulses.json'))
    dat=[z for z in read(ROOT/'data/observed_calibration_Q20_Q40_shell6_15.csv') if K(int(z['degree'])%5,int(z['shell']))-int(z['degree'])//5<=5]
    f=rowsys(dat,phase,period,prior)
    pivots={};bad=[]
    for row in dat:
        d,t=int(row['degree']),int(row['shell'])
        co,b,_,r=f(d,t,row['pt'])
        co,b=reduce(co,b,pivots)
        if co:
            p=min(co);c=co[p];co={k:v/c for k,v in co.items()};b=b/c
            pivots[p]=(co,b)
        elif b:bad.append((d,t,str(b),r))
    allvars=set()
    for row in dat:
        eq=f(int(row['degree']),int(row['shell']),row['pt'])
        allvars.update(eq[0])
    sites=read(ROOT/'data/withheld_sites_Q41_Q53_shell11_15.csv'); preds=[]
    for site in sites:
        d,t=int(site['degree']),int(site['shell'])
        if K(d%5,t)-d//5>5:
            preds.append(dict(degree=d,shell=t,prediction='',unfixed_variables='outside_profile_support',r=K(d%5,t)-d//5))
            continue
        co,b,_,r=f(d,t)
        rr,val=reduce(co,b,pivots)
        # value of linear expression = b? eq: A - offset = sum unknown
        # reduce target a - b by eliminating pivots; if no remaining, original A = offset + (-b)?
        # Since rhs b=-offset, row represents sum unknown=-offset + A; identify A from reduced equation.
        # b transform is -offset - linear combination of training rhs -> if 0 remaining, A=-b.
        if not rr:
            A=-val
            w=W(d%5,t,d//5);raw=A*phase[d%5,(d//5)%period]*w
            pred=str(raw.numerator) if raw.denominator==1 else ''
        else:pred=''
        preds.append(dict(degree=d,shell=t,prediction=pred,unfixed_variables=len(rr),r=r))
    print('training rows',len(dat),'pivot rank',len(pivots),'distinct variables',len(allvars),'train conflicts',len(bad),'unique predicted',sum(bool(x['prediction']) for x in preds))
    print('first conflicts',bad[:12]); print('preds', [z for z in preds if z['prediction']][:15])
    out=ROOT/'results';out.mkdir(exist_ok=True)
    write(out/'linear_profile_r5_predictions.csv',preds)
    (out/'linear_profile_r5_audit.json').write_text(json.dumps({'rank':len(pivots),'n_variables':len(allvars),'training_conflicts':bad,'unfixed_site_count':sum(not x['prediction'] for x in preds)},indent=2))
    return preds
if __name__=='__main__':run()
