#!/usr/bin/env python3
"""Read withheld labels in a SEPARATE validation program. No changes to predictor."""
from __future__ import annotations
import csv,json,math
from collections import defaultdict,Counter
from pathlib import Path
from recursive_defects import read,groups,phase_fit,write,nearest_integer,evalpoly
ROOT=Path(__file__).resolve().parent

def frac_err(yhat,y):return abs(yhat-y)/max(1,abs(y))

def val5(n):
    n=abs(n)
    if n==0:return None
    v=0
    while n%5==0:n//=5;v+=1
    return v

def check():
    cal=read(ROOT/'data/calibration_Q1_Q20.csv')
    update=read(ROOT/'data/update_Q21_Q25.csv')
    refs=read(ROOT/'data/reference_only_Q26_Q30.csv')
    pred=read(ROOT/'results/predictions_Q26_Q30.csv')
    period,phase=phase_fit(cal)
    reference={(int(x['degree']),int(x['shell'])):int(x['pt']) for x in refs}
    assert len(pred)==len(reference)==24
    scores={col:{'total':0,'exact':0,'sum_relative_error':0.,'max_relative_error':0.,'by_shell':defaultdict(lambda:[0,0])}
            for col in ['quadratic_PT','cubic_PT','linear_times_quadratic_PT','mirror_PT']}
    details=[];bounds={'covered':0,'total':0,'zero_prev_count':0,'zero_prev_correct':0,'persistent_5adic':0,'persistent_count':0,'persistent_attenuation':0,'sign_flips':0}
    for p in pred:
        d=int(p['degree']);t=int(p['shell']);y=reference[(d,t)]
        phase_s=int(p['sign'])
        vals={}
        for col,score in scores.items():
            if p[col] in ('',None):continue
            v=int(p[col]);vals[col]=v
            rel=frac_err(v,y);score['total']+=1;score['exact']+=int(v==y)
            score['sum_relative_error']+=rel;score['max_relative_error']=max(score['max_relative_error'],rel)
            score['by_shell'][t][0]+=int(v==y);score['by_shell'][t][1]+=1
        P=cal+update
        if p['previous_defect'] not in ('',None):
            prev=int(p['previous_defect'])
            y_signed=phase_s*y
            g=groups(P,phase,period)
            own=g[d%5,t]
            e5=y_signed-3*own[-1][1]+3*own[-2][1]-own[-3][1]
            lo=int(p['defect_interval_low']);hi=int(p['defect_interval_high'])
            bounds['total']+=1;bounds['covered']+=int(lo<=e5<=hi)
            if prev==0:
                bounds['zero_prev_count']+=1;bounds['zero_prev_correct']+=int(e5==0)
            if prev!=0 and e5!=0:
                bounds['persistent_count']+=1
                bounds['persistent_5adic']+=int(val5(e5)>=val5(prev))
                bounds['persistent_attenuation']+=int(abs(e5)<=abs(prev))
                bounds['sign_flips']+=int((e5>0)!=(prev>0))
        details.append(dict(degree=d,residue=d%5,shell=t,actual_PT=y,**vals))
    group=groups(cal+update+refs,phase,period)
    support={};pulses={'zero':0,'one':0,'two':0}
    series=[];pair_locks={}
    for (s,t),seq in sorted(group.items()):
        # Use precisely 5 entries per stream starting degree 6.
        cut=[(k,v) for k,v in seq if 5*k+s>=6]
        if len(cut)!=5:raise AssertionError((s,t,len(cut)))
        x=[v for k,v in cut];k0=cut[0][0]
        j1=k0+3;j2=k0+4
        eps1=x[3]-3*x[2]+3*x[1]-x[0]
        eps2=x[4]-3*x[3]+3*x[2]-x[1]
        Q2=evalpoly(cut[:3],cut[-1][0])
        # exact Newton/spline reconstruction without floating-point
        reconstructed=Q2+eps1*math.comb(cut[-1][0]-j1+2,2)+eps2*math.comb(cut[-1][0]-j2+2,2)
        assert reconstructed==cut[-1][1]
        activity=(eps1!=0,eps2!=0)
        pulses['zero' if activity==(False,False) else 'one' if activity==(True,False) else 'two']+=1
        assert not(activity==(False,True)),(s,t)
        series.append(dict(residue=s,shell=t,first_k=k0,defect_at_first_update=eps1,defect_at_next_update=eps2,
                           defect_front_state=('locked' if not eps1 else 'healed' if not eps2 else 'persistent'),
                           val5_first=val5(eps1),val5_second=val5(eps2)))
        for idx,ek in [(3,eps1),(4,eps2)]:
            k=cut[idx][0];support.setdefault((s,k),[])
            if ek:support[s,k].append(t)
    # exact agreement of paired residue curvatures at Q<=25 and Q<=30
    for pair in ((1,4),(2,3)):
        a,b=pair
        for time,label in ((-2,'update'),(-1,'test')):
            equal=[]
            for t in range(min(max(tt for ss,tt in group if ss==a),max(tt for ss,tt in group if ss==b))+1):
                q1=[v for k,v in group[a,t]];q2=[v for k,v in group[b,t]]
                # position -2 is k=4 for each pair; -1 is k=5
                cur1=q1[time]-2*q1[time-1]+q1[time-2]
                cur2=q2[time]-2*q2[time-1]+q2[time-2]
                if cur1==cur2:equal.append(t)
            pair_locks[f'{a}<->{b}:{label}']=equal
    out_scores={}
    for col,ss in scores.items():
        out_scores[col]={k:v for k,v in ss.items() if k!='by_shell'}
        out_scores[col]['mean_relative_error']=ss['sum_relative_error']/ss['total'] if ss['total'] else None
        out_scores[col]['by_shell']={str(k):v for k,v in ss['by_shell'].items()}
    result={'n_streams':len(group),'defect_pulses':pulses,'support':{f'{s},{k}':v for (s,k),v in sorted(support.items())},
            'pair_locks':pair_locks,'defect_interval':bounds,'model_scores':out_scores,
            'validation_has_access_to_reference':True,'predictor_has_access_to_reference':False}
    (ROOT/'results/validated_metrics.json').write_text(json.dumps(result,indent=2))
    write(ROOT/'results/defect_shells.csv',series)
    write(ROOT/'results/predictions_scored_Q26_Q30.csv',details)
    print('Outcome:',json.dumps(result,indent=2))
    return result
if __name__=='__main__':check()
