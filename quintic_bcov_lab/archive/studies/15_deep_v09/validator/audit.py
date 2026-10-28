#!/usr/bin/env python3
"""Independent conditional reference evaluation; never imported by predictor."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
import csv,json,sys,math
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from recursive_defects import phase_fit,read,groups
from gated_recurrence import C,PEER,envelope,front
from deep_reference import evaluate

def polynomial_coeff(points):
    # 3-point exact interpolation; coefficients (a,b,c) in a*k^2+b*k+c.
    (x0,y0),(x1,y1),(x2,y2)=points
    assert x0!=x1 and x0!=x2 and x1!=x2
    def term(i):
        x,y=points[i]; xx=[u for j,(u,_) in enumerate(points) if i!=j]
        v=F(y,(x-xx[0])*(x-xx[1]))
        return v, -v*(xx[0]+xx[1]),v*xx[0]*xx[1]
    parts=[term(i) for i in range(3)]
    return tuple(sum(p[j] for p in parts) for j in range(3))

def build_metrics():
    shallow=read(ROOT/'data/shallow_calibration_Q1_Q20.csv')
    period,phase=phase_fit(shallow);assert period==4
    earlier=read(ROOT/'data/shallow_observed_Q6_Q30.csv')
    G=groups(earlier,phase,period)
    for s in range(5):
        for t in (0,1):
            pts=G[s,t][-3:]; a,b,c=polynomial_coeff(pts)
            assert a!=0
            assert F(5)*b/a==2*s-5,(s,t,'b/a')
            assert F(5)*c/a==C[s]-2*t,(s,t,'c/a')
    reference=read(ROOT/'validator/reference_Q20_Q53.csv')
    old=read(ROOT/'validator/shallow_reference_Q31_Q53.csv')
    assert len(old)==109
    for r in old:
        got=evaluate(int(r['degree']),int(r['shell']))
        assert got['pt']==r['pt'] and got['status']=='certified_within_conditional_framework'
    allsafe=[r for r in reference if r['pt']]
    assert len(allsafe)==116
    ref={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in allsafe}
    assert len([x for x in ref if 31<=x[0]<=53])==106
    assert len([x for x in ref if 41<=x[0]<=53])==65
    gated=read(ROOT/'results/frozen_gated_predictions_Q41_Q53.csv')
    assert len(gated)==65
    states=Counter(r['state'] for r in gated)
    exact=sum(int(r['prediction'])==ref[int(r['degree']),int(r['shell'])] for r in gated if r['prediction'])
    assert states=={'mature_recurrence':32,'active_front':19,'needs_mature_anchor':14}
    assert exact==32,(states,exact)
    ungated=read(ROOT/'results/predictions_paired_repeated.csv')
    exact_u=sum(int(r['prediction'])==ref[int(r['degree']),int(r['shell'])] for r in ungated if r['prediction'])
    assert exact_u==28,exact_u
    group=defaultdict(dict)
    for r in allsafe:
        d=int(r['degree']);t=int(r['shell']);s=d%5;k=d//5
        w=envelope(s,t,k)
        assert w!=0
        group[(s,t)][k]=F(int(r['pt']),phase[s,k%period]*w)
    trans=[];dd=[];known_mature=defaultdict(list)
    for (s,t),v in sorted(group.items()):
        for k in sorted(v):
            if k>=front(s,t):known_mature[s,t].append((k,v[k]))
            if k-1 in v:
                actual=v[k]!=v[k-1]
                predicted=k<=front(s,t)
                trans.append(dict(residue=s,shell=t,k=k,degree=5*k+s,
                                  activity=int(actual),front_prediction=int(predicted),
                                  difference_num=str((v[k]-v[k-1]).numerator),
                                  difference_den=str((v[k]-v[k-1]).denominator)))
                assert actual==predicted,(s,t,k)
            if all(k-i in v for i in range(4)):
                # Normalize back to sign-corrected PT before third difference.
                yy={kk: v[kk]*envelope(s,t,kk) for kk in range(k-3,k+1)}
                third=yy[k]-3*yy[k-1]+3*yy[k-2]-yy[k-3]
                predicted=k <= (t+1 if s==0 else t+5-s)
                actual=third!=0
                assert actual==predicted,(s,t,k,third)
                dd.append(dict(residue=s,shell=t,k=k,degree=5*k+s,
                   activity=int(actual),front_prediction=int(predicted),third_difference=str(third)))
    assert len(trans)==91 and len(dd)==41,(len(trans),len(dd))
    assert sum(r['activity'] for r in trans)==47
    assert sum(r['activity'] for r in dd)==28
    assert all(len({a for k,a in li})==1 for li in known_mature.values())
    assert sum(len(li)>=2 for li in known_mature.values())==18
    pairchecks=[]
    for t in range(6,11):
        for s,z in ((1,4),(2,3)):
            if known_mature.get((s,t)) and known_mature.get((z,t)):
                a=known_mature[s,t][-1][1];b=known_mature[z,t][-1][1]
                assert a==b
                pairchecks.append((s,z,t))
    assert len(pairchecks)==7
    # Persistent quotient defects: shrinkage/phase regularity, tested but not enforced.
    persist=Counter()
    for st,v in group.items():
        diff={k:v[k]-v[k-1] for k in v if k-1 in v}
        for k in diff:
            if k-1 in diff and diff[k] and diff[k-1]:
                persist['count']+=1
                persist['shrink']+=abs(diff[k])<abs(diff[k-1])
                persist['sign_reverse']+=diff[k]*diff[k-1]<0
    assert persist=={'count':24,'shrink':24,'sign_reverse':19}
    for name,rows in [('normalized_scale_defects.csv',trans),('third_difference_front.csv',dd)]:
        out=ROOT/'results'/name
        with out.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(train_scope='shallow models Q<=30; deep scale calibration Q<=40',
      test_scope='Q41..53 and new shell depths 6..10',
      reference_assumptions='conditional rank-zero DT modular ansatz and PT wall-crossing; not predictor inputs',
      shallow_reference_rechecked=len(old),deep_reference_safe_Q31_Q53=106,
      deep_reference_incomplete_Q31_Q40=9,deep_reference_safe_Q41_Q53=65,
      frozen_gated_exact=exact,frozen_gated_predictions=states['mature_recurrence'],
      frozen_gated_abstentions=65-states['mature_recurrence'],
      frozen_ungated_exact=exact_u,frozen_ungated_total=65,
      quotient_defect_activity_accuracy=f'{len(trans)}/{len(trans)}',
      quotient_defect_holdout_accuracy=str(sum(r['degree']>=41 for r in trans))+'/'+str(sum(r['degree']>=41 for r in trans)),
      quotient_defect_nonzero=sum(r['activity'] for r in trans),
      third_defect_activity_accuracy=f'{len(dd)}/{len(dd)}',
      third_defect_holdout_accuracy=str(sum(r['degree']>=41 for r in dd))+'/'+str(sum(r['degree']>=41 for r in dd)),
      third_defect_nonzero=sum(r['activity'] for r in dd),
      mature_scales_total=sum(len(li) for li in known_mature.values()),
      mature_scales_sequences=sum(bool(li) for li in known_mature.values()),
      mature_scales_repeated_sequences=sum(len(li)>=2 for li in known_mature.values()),
      equal_mature_paired_scales=len(pairchecks),persistent_quotient_defects=dict(persist))
    (ROOT/'results/metrics.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=='__main__':build_metrics()
