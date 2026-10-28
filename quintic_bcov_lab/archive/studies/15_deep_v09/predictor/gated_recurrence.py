"""No-oracle, exact-rational PT recurrence: frozen v0.8 diagonal-front hypothesis.

No modularity, DT or wall-crossing dependencies; only shallow PT training,
calibration observations Q<=40 and unlabelled degree/shell sites.

The normalized scale is PT / (phase * W). A shell's scale is *predicted*
to stop changing at K(s,t). Before that point this module abstains.
"""
from fractions import Fraction as Fraction
from pathlib import Path
from recursive_defects import phase_fit,read,write

C={0:10,1:6,2:4,3:4,4:6}
PEER={1:4,4:1,2:3,3:2}

def front(s,t):
    return t-1 if s==0 else t+3-s

def envelope(s,t,k):
    return 5*k*k+(2*s-5)*k+C[s]-2*t

def forecast_one(train,degree,shell,phase,period):
    s=degree%5;t=shell;k=degree//5
    if k<front(s,t):
        return dict(degree=degree,shell=t,prediction='',state='active_front',source_degree='',source_residue='',scale='')
    opts=[]
    for z in (s,PEER.get(s)):
        if z is None:continue
        for r in train:
            d=int(r['degree']);tt=int(r['shell']);kk=d//5
            if not r.get('pt') or d%5!=z or tt!=t or kk<front(z,t):continue
            w=envelope(z,t,kk)
            if not w:continue
            scale=Fraction(int(r['pt']),phase[z,kk%period]*w)
            opts.append((kk,d,z,scale))
    if not opts:
        return dict(degree=degree,shell=t,prediction='',state='needs_mature_anchor',source_degree='',source_residue='',scale='')
    # Latest mature anchor, with own-residue tie preference.
    anchor=max(opts,key=lambda q:(q[0],q[2]==s))
    kk,d,z,a=anchor
    val=a*phase[s,k%period]*envelope(s,t,k)
    if val.denominator!=1:
        # Do not make an uncited nearest-integer guess for a certified exact site.
        return dict(degree=degree,shell=t,prediction='',state='nonintegral_candidate',source_degree=d,source_residue=z,scale=str(a))
    return dict(degree=degree,shell=t,prediction=str(val.numerator),state='mature_recurrence',source_degree=d,source_residue=z,scale=str(a))

def run(root=None):
    root=Path(root) if root is not None else Path(__file__).resolve().parent.parent
    shallow=read(root/'data/shallow_calibration_Q1_Q20.csv')
    period,phase=phase_fit(shallow)
    assert period==4,period
    train=read(root/'data/calibration_deep_Q20_Q40.csv')
    assert train and max(int(r['degree']) for r in train)<=40
    sites=read(root/'data/forecast_sites_Q41_Q53.csv')
    preds=[forecast_one(train,int(r['degree']),int(r['shell']),phase,period) for r in sites]
    out=root/'results/frozen_gated_predictions_Q41_Q53.csv'
    write(out,preds)
    return preds

if __name__=='__main__':
    preds=run()
    from collections import Counter
    print(Counter(x['state'] for x in preds))
