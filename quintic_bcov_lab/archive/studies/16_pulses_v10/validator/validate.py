#!/usr/bin/env python3
"""Score predictions after the empirical learner has finished.

The exact labels in this validator are conditional on the earlier independent
rank-zero DT / PT wall-crossing reference construction. The predictor never
imports this file or the reference table.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import csv,json,statistics,sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from recursive_defects import read,write,phase_fit
from deep_envelope import build_amplitudes
from pulse_profile import K

def phase_and_values(rows):
    period,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))
    out=build_amplitudes(rows,phase,period)
    return {(s,t,k):v for (s,t),seq in out.items() for k,v in seq}

def check_front(rows,holdout_min=41):
    a=phase_and_values(rows)
    D={(s,t,k):v-a[s,t,k-1] for (s,t,k),v in a.items() if (s,t,k-1) in a}
    pairs=[]
    for (s,t,k),v in sorted(D.items()):
        nxt=D.get((s,t+1,k+1))
        if nxt is None or (v==0 and nxt==0):continue
        d=5*(k+1)+s
        pairs.append({'degree':d,'residue':s,'shell':t,'k':k,
           'at_front':int(k==K(s,t)), 'pulse_before':str(v),
           'pulse_after':str(nxt),'secondary_residual':str(v+nxt),
           'exact_transport':int(v+nxt==0)})
    return pairs

def score():
    ref={(int(z['degree']),int(z['shell'])):int(z['pt'])
         for z in read(ROOT/'validator/reference_Q20_Q53.csv') if z['pt']}
    pred=read(ROOT/'results/pulse_predictions_Q41_Q53.csv')
    baseline=read(ROOT/'validator/previous_v09_gated_predictions.csv')
    assert len(pred)==len(baseline)==65
    scored=[]
    for z in pred:
        d,t=int(z['degree']),int(z['shell'])
        assert (d,t) in ref,(d,t)
        z=dict(z)
        z['reference']=str(ref[d,t])
        z['exact']=int(bool(z['prediction']) and int(z['prediction'])==ref[d,t])
        z['error']=str(int(z['prediction'])-ref[d,t]) if z['prediction'] else ''
        scored.append(z)
    pcount=sum(bool(z['prediction']) for z in scored)
    exact=sum(z['exact'] for z in scored)
    assert (pcount,exact)==(59,59),(pcount,exact)
    by_shell={str(t):{'predicted':sum(int(z['shell'])==t and bool(z['prediction']) for z in scored),
                      'exact':sum(int(z['shell'])==t and z['exact'] for z in scored),
                      'total':sum(int(z['shell'])==t for z in scored)} for t in range(6,11)}
    assert sum(bool(z['prediction']) for z in baseline)==32
    assert all(not z['prediction'] or int(z['prediction'])==ref[int(z['degree']),int(z['shell'])] for z in baseline)
    base={(int(z['degree']),int(z['shell'])) for z in baseline if z['prediction']}
    new=[z for z in scored if z['prediction'] and (int(z['degree']),int(z['shell'])) not in base]
    assert len(new)==27
    pairs_train=check_front(read(ROOT/'data/calibration_deep_Q20_Q40.csv'))
    pairs_ref=check_front(read(ROOT/'validator/reference_Q20_Q53.csv'))
    assert len(pairs_train)==6 and len(pairs_ref)==32
    assert sum(z['at_front'] for z in pairs_train)==4
    assert sum(z['exact_transport'] for z in pairs_train)==4
    assert sum(z['at_front'] for z in pairs_ref)==15
    assert sum(z['exact_transport'] for z in pairs_ref)==15
    assert all(z['exact_transport']==z['at_front'] for z in pairs_ref)
    write(ROOT/'validator/scored_predictions.csv',scored)
    write(ROOT/'validator/front_transport_pairs.csv',pairs_ref)
    # Retrospective, true-history secondary-pulse diagnostics only. DO NOT
    # feed these residuals to the no-oracle predictor.
    a=phase_and_values(read(ROOT/'validator/reference_Q20_Q53.csv'))
    D={(s,t,k):v-a[s,t,k-1] for (s,t,k),v in a.items() if (s,t,k-1) in a}
    J={(s,t,k):v+D[s,t+1,k+1] for (s,t,k),v in D.items() if (s,t+1,k+1) in D}
    secondary=[]
    for (s,t,k),v in sorted(J.items()):
        nxt=J.get((s,t+1,k+1))
        d=5*(k+2)+s
        if nxt is None or d<46 or not v or not nxt:continue
        guess=-F(k,k+4)*v
        secondary.append({'degree':d,'residue':s,'shell':t,'k':k,'observed_next_J':str(nxt),
            'predicted_next_J':str(guess),'zero_J_error':'1',
            'relative_secondary_error':str(abs(guess-nxt)/abs(nxt))})
    assert len(secondary)==8
    write(ROOT/'validator/secondary_pulse_diagnostics.csv',secondary)
    # The affine-in-inverse-W profile has a closed formula for J. These are
    # retrospective transformed-domain checks, not independent extra labels.
    fitted=json.loads((ROOT/'results/fitted_pulses.json').read_text())
    profiles={(z['residue'],z['front_distance']):(F(z['P']),F(z['Q'])) for z in fitted['profiles']}
    Wfun=lambda s,t,k: 5*k*k+(2*s-5)*k+{0:10,1:6,2:4,3:4,4:6}[s]-2*t
    jchecks=[]
    for (s,t,k),value in sorted(J.items()):
        degree=5*(k+1)+s
        if not value or degree<41:continue
        r=K(s,t)-k
        if r<1 or (s,r) not in profiles or (s,r+1) not in profiles:continue
        Qr=profiles[s,r][1];Qnext=profiles[s,r+1][1]
        reconstructed=(-1 if t%2 else 1)*(Qr*(F(1,Wfun(s,t,k))-F(1,Wfun(s,t+1,k+1)))+
                               Qnext*(F(1,Wfun(s,t+1,k))-F(1,Wfun(s,t,k-1))))
        jchecks.append({'degree':degree,'residue':s,'shell':t,'distance':r,
                        'observed_J':str(value),'model_J':str(reconstructed),
                        'exact':int(value==reconstructed)})
    assert len(jchecks)==14 and all(z['exact'] for z in jchecks)
    write(ROOT/'validator/secondary_closed_form_checks.csv',jchecks)
    secondary_mean=statistics.mean(float(F(z['relative_secondary_error'])) for z in secondary)
    summary={
        'reference_status':'conditional modular DT / wall-crossing, scorer only',
        'hypothesis_status':'exploratory; feature class selected with inspection of extended data',
        'train_cutoff_degree':40,'test_degree_range':[41,53],
        'sites':65,'v09_frozen_exact':32,'v09_frozen_abstentions':33,
        'v10_predicted':pcount,'v10_exact':exact,'v10_abstentions':65-pcount,
        'new_exact_sites':len(new),'by_depth':by_shell,
        'unresolved':[{'degree':int(z['degree']),'shell':int(z['shell'])} for z in scored if not z['prediction']],
        'calibration_diagonal':{'nonzero_pairs':len(pairs_train),'on_front':4,'on_front_exact':4,'off_front_exact':0},
        'full_diagonal':{'nonzero_pairs':len(pairs_ref),'on_front':15,'on_front_exact':15,'off_front_exact':0},
        'closed_secondary_residual':{'exact':len(jchecks),'tested':len(jchecks),'independence':'algebraic consequence of fitted profile'},
        'secondary_pulse_teacher_forced':{'count':len(secondary),'zero_baseline_mean_relative_error':1.0,
                                           'rational_ratio_mean_relative_error':secondary_mean,
                                           'not_an_autonomous_test':True},
        'source_tiers':dict(Counter(z['source'] or 'abstain' for z in scored)),
    }
    (ROOT/'validator/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return summary

if __name__=='__main__':score()
