#!/usr/bin/env python3
"""Score the frozen v0.7 predictor without modifying it.

All predictor calls happen before exposing the target block. In teacher-forced
mode the previous block is revealed only *after* its predictions are written.
Independent reference CSV originates from generate_reference.py.
"""
import csv, json, sys, statistics
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
import recursive_defects as r
from rollout import sites, DEPTH, run as autonomous_rollout

def read(path):return r.read(path)
def write(path,rows):r.write(path,rows)

def metrics(pred,oracle,key):
    arr=[];rat=[];matches=0
    for z in pred:
        d,t=int(z['degree']),int(z['shell']);p=z.get(key,'')
        if p in (None,''):continue
        v=int(oracle[d,t]); y=int(p);matches+=int(v==y)
        arr.append((d,t,y,v,y-v));rat.append(Fraction(abs(y-v),max(1,abs(v))))
    return {'n':len(arr),'exact':matches,'pct_exact':100*matches/len(arr),
            'mean_abs_rel_pct':100*statistics.mean(float(q) for q in rat),
            'median_abs_rel_pct':100*statistics.median(float(q) for q in rat),
            'max_abs_rel_pct':100*max(float(q) for q in rat)}

def block(degree):return f'{31+5*((degree-31)//5)}-{min(35+5*((degree-31)//5),53)}'

def main():
    oracle={(int(z['degree']),int(z['shell'])):int(z['pt']) for z in read(ROOT/'results/reference_Q31_Q53.csv')}
    assert len(oracle)==109
    cal=read(ROOT/'predictor/data/calibration_Q1_Q20.csv')
    train=read(ROOT/'predictor/data/observed_Q6_Q30.csv')
    period,phase=r.phase_fit(cal);assert period==4
    teacher=[]
    for start in range(31,54,5):
        end=min(start+4,53)
        new=r.predict(train,sites(start,end),period,phase)
        teacher+=new
        # This happens AFTER the entire degree block was predicted.
        for z in new:
            d,t=int(z['degree']),int(z['shell'])
            train.append({'degree':str(d),'shell':str(t),'pt':str(oracle[d,t])})
    auto=autonomous_rollout()
    assert len(teacher)==len(auto)==109
    old={(int(z['degree']),int(z['shell'])):z for z in read(ROOT/'results/v07_frozen_Q31_Q35.csv')}
    assert len(old)==24
    for z in auto:
        d,t=int(z['degree']),int(z['shell'])
        if d<=35:assert int(z['mirror_PT'])==int(old[d,t]['mirror_PT'])
    write(ROOT/'results/teacher_forced_predictions_Q31_Q53.csv',teacher)
    report={'reference_assumptions':'rankzero DT modular ansatz plus wall crossing, used only for scoring',
            'calibration_through_degree':30,'no_refits':True,
            'teacher':{},'autonomous':{},'baselines':{},'by_block':[]}
    for key in ('mirror_PT','quadratic_PT','cubic_PT','linear_times_quadratic_PT'):
        report['baselines'][key]=metrics(teacher,oracle,key)
    report['teacher']=metrics(teacher,oracle,'mirror_PT')
    report['autonomous']=metrics(auto,oracle,'mirror_PT')
    scored=[]
    for mode,pred in [('teacher',teacher),('autonomous',auto)]:
        for z in pred:
            d,t=int(z['degree']),int(z['shell']); actual=oracle[d,t]
            choice=int(z['mirror_PT'])
            scored.append(dict(mode=mode,degree=d,shell=t,block=block(d),
                               predicted_pt=str(choice),reference_pt=str(actual),
                               integer_error=str(choice-actual),exact=int(choice==actual),
                               abs_relative_error=str(Fraction(abs(choice-actual),max(1,abs(actual)))),
                               pulse_state=z['pulse_state'],mirror_used=z['mirror_used']))
        for start in range(31,54,5):
            end=min(start+4,53)
            subset=[z for z in pred if start<=int(z['degree'])<=end]
            v=metrics(subset,oracle,'mirror_PT')
            report['by_block'].append({'mode':mode,'degree_start':start,'degree_end':end,**v})
    write(ROOT/'results/individual_predictions_scored.csv',scored)

    # Freeze sign correction and compute exact third-difference defects.
    truth=read(ROOT/'predictor/data/observed_Q6_Q30.csv') + [
        {'degree':str(d),'shell':str(t),'pt':str(v)} for (d,t),v in oracle.items()]
    streams=r.groups(truth,phase,period)
    defect_rows=[]
    for (s,t),seq in sorted(streams.items()):
        for j in range(3,len(seq)):
            k,v=seq[j]
            if seq[j-1][0]!=k-1 or seq[j-2][0]!=k-2 or seq[j-3][0]!=k-3:continue
            d=5*k+s
            eps=v-3*seq[j-1][1]+3*seq[j-2][1]-seq[j-3][1]
            defect_rows.append(dict(degree=d,residue=s,depth=t,k=k,defect=str(eps),active=int(eps!=0)))
    write(ROOT/'results/third_differences.csv',defect_rows)
    supports=[];prev_seen={}
    revivals=[]
    for start in range(26,54,5):
        end=min(start+4,53)
        present=[z for z in defect_rows if start<=int(z['degree'])<=end and int(z['active'])]
        supports.append({'degree_start':start,'degree_end':end,'n_active':len(present),
                         'active_sites':[(int(z['residue']),int(z['depth'])) for z in present]})
    for z in sorted(defect_rows,key=lambda z:(int(z['residue']),int(z['depth']),int(z['k']))):
        if int(z['degree'])<31:continue
        key=(int(z['residue']),int(z['depth']));past=[x for x in defect_rows if (int(x['residue']),int(x['depth']))==key and int(x['degree'])==int(z['degree'])-5]
        if past and not int(past[0]['active']) and int(z['active']):revivals.append((int(z['degree']),key))
    report['defect_support_by_block']=supports
    report['defect_reactivations']=revivals
    assert not revivals, revivals
    assert [x['n_active'] for x in supports]==[9,4,2,1,0,0],supports
    report['frozen_reference_check_Q31_Q35']='24 predictions identical to v0.7 archive'
    (ROOT/'results/metrics.json').write_text(json.dumps(report,indent=2))
    print('SEQUENTIAL',report['teacher'])
    print('AUTONOMOUS',report['autonomous'])
    print('DEFECTS',supports)
    print('ALL OUTPUTS VERIFIED')

if __name__=='__main__':main()
