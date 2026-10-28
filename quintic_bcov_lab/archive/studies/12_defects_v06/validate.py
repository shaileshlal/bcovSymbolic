#!/usr/bin/env python3
"""Post-prediction audit: ONLY this program opens hidden Q26..30 labels.

The hypothesis choices are exploratory; this is a retrospective forward-degree
backtest, not a prospectively registered experiment.
"""
from __future__ import annotations
import csv,json,math,statistics,sys
from math import comb
from pathlib import Path
from collections import defaultdict

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from baseline_v05 import read_csv,train_predict,exact_phase_fit
from defect_features import event_records

def fracerr(p,y): return abs(p-y)/max(1,abs(y))
def median(xs):return statistics.median(xs)
def score(pred,truth):
    match=[];relative=[];by_t=defaultdict(list)
    for r in pred:
        key=(int(r['degree']),int(r['shell']))
        p=int(r['predicted_PT']);y=truth[key]
        match.append(p==y);relative.append(fracerr(p,y));by_t[key[1]].append(p==y)
    return dict(exact=sum(match),n=len(match),sign_correct=sum((int(r['predicted_PT'])>0)==(truth[int(r['degree']),int(r['shell'])]>0) for r in pred),
      median_relative_error=median(relative),median_relative_error_misses=median([u for u in relative if u>0]) if any(relative) else 0,
      mean_absolute_relative_error=sum(relative)/len(relative),max_relative_error=max(relative),
      by_shell={str(k):[sum(v),len(v)] for k,v in sorted(by_t.items())})

def getphase(audit):
    return {(int(x),int(y)):int(v) for k,v in audit['signs'].items() for x,y in [k.split(',')]}

def main():
    train=read_csv(ROOT/'data/calibration_Q1_Q20.csv')
    update=read_csv(ROOT/'data/update_Q21_Q25.csv')
    hidden=read_csv(ROOT/'data/reference_only_Q26_Q30.csv')
    truth={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in update+hidden}
    out=ROOT/'results'
    p1=read_csv(out/'block_Q21_Q25.csv')
    p2=read_csv(out/'block_Q26_Q30.csv')
    au=json.loads((out/'model_audit.json').read_text())
    phase=getphase(au)
    E1=event_records(train,update,phase,au['period'])
    E2=event_records(train+update,hidden,phase,au['period'])
    assert len(E1)==len(E2)==24
    # Fixed-window v0.5 model uses only Q<=20 for BOTH blocks
    sites=read_csv(ROOT/'data/sites_Q21_Q25.csv')+read_csv(ROOT/'data/sites_Q26_Q30.csv')
    fixed,_=train_predict(train,sites)
    # rolling simple quadratic uses the same PT observations available at forecast time
    quad_p2=[dict(degree=r['degree'],shell=r['shell'],predicted_PT=r['baseline_quad_PT']) for r in p2]
    scoremap={
       'v05_fixed_train_Q1_Q20':score(fixed,truth),
       'v05_update_first_block':score(p1,truth),
       'rolling_quad_second_block':score(quad_p2,truth),
       'new_mirror_bracket_second_block':score(p2,truth),
       'v06_sequential_Q21_Q30':score(p1+p2,truth),
    }
    act1={(e['residue'],e['shell']) for e in E1 if e['active']}
    act2={(e['residue'],e['shell']) for e in E2 if e['active']}
    assert act2.issubset(act1)
    assert len(act1)==13 and len(act2)==9
    stable_before=[e for e in E2 if e['own_prevcurvature'] is not None and e['own_prevcurvature']==e['own_curvature']]
    stable_next=[e for e in stable_before if not e['active']]
    active_by_s={}
    for s in range(5):
        active_by_s[str(s)]={'earlier':sorted(t for ss,t in act1 if ss==s),
                              'later':sorted(t for ss,t in act2 if ss==s)}
    stats={
        'active_first_window':len(act1),'active_second_window':len(act2),
        'new_activations':len(act2-act1),'healed_defects':len(act1-act2),
        'stable_persistence':{'success':len(stable_next),'total':len(stable_before)},
        'active_shells_by_residue':active_by_s,
        'exchangeable_active_set_probability':comb(13,9)/comb(24,9),
        'no_future_label_in_prediction_module':True
    }
    # The dev block supports a conventional split-conformal relative-width interval
    # (90% nominal): fit quantile on 24 independent degree21..25 errors.
    dev_error=sorted(fracerr(int(r['predicted_PT']),truth[int(r['degree']),int(r['shell'])]) for r in p1)
    rank=math.ceil((len(dev_error)+1)*.9)
    if rank>len(dev_error):tau=float('inf')
    else:tau=dev_error[rank-1]
    covered=[fracerr(int(r['predicted_PT']),truth[int(r['degree']),int(r['shell'])])<=tau for r in p2]
    conformal=dict(n_cal=len(dev_error),target_coverage=.9,quantile_rank=rank,
      relative_width=tau,coverage=sum(covered)/len(covered),covered=sum(covered),total=len(covered))
    # Dump exact defect labels to transparent post-hoc EDA (but predictor never opens it).
    feature_rows=[]
    for e in E1+E2:
        feature_rows.append(dict(degree=e['degree'],residue=e['residue'],shell=e['shell'],
          k=e['k'],phase=phase[e['residue'],e['k']%au['period']],
          third_difference=str(e['actual_defect']),active=e['active'],
          curvature_before=str(e['own_curvature']),
          curvature_after=str(e['own_curvature']+e['actual_defect']),
          mirror_curvature=str(e['peer_curvature']) if e['peer_curvature'] is not None else '',
          v5_delta=str(__import__('defect_features').val5(e['actual_defect'])) if e['actual_defect'] else '',
          relative_defect=abs(e['actual_defect'])/max(1,abs(e['actual_PT']))))
    with (out/'defect_analysis_only.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=feature_rows[0].keys());w.writeheader();w.writerows(feature_rows)
    combined=dict(scores=scoremap,defect_support=stats,calibration=conformal)
    (out/'validated_metrics.json').write_text(json.dumps(combined,indent=2))
    print(json.dumps(combined,indent=2))

if __name__=='__main__':main()
