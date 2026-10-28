#!/usr/bin/env python3
"""No-oracle sequential prediction; produces only forecasts, not validation metrics.

Q1..20 data train the period-four and v0.5 quadratic/mirror rule.
The independent Q21..25 block is an observational update before Q26..30.
No code path in this module opens reference_only_Q26_Q30.csv.
"""
from __future__ import annotations
import sys, csv, json, hashlib, argparse
from pathlib import Path
from fractions import Fraction
from baseline_v05 import read_csv,train_predict,exact_phase_fit,nearest_int
from defect_features import grouped_signed_rows,feat_from_state,DELTA0,event_records,PEER

ROOT=Path(__file__).resolve().parent

def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def predict_second_block(calibration,update,sites,signs,period=4):
    groups=grouped_signed_rows(calibration+update,signs,period)
    out=[]
    for rr in sites:
        d=int(rr['degree']);t=int(rr['shell']);s=d%5;k=d//5
        o=feat_from_state(groups,s,t)
        if o is None:raise ValueError(f'No prior observations for {(s,t)}')
        own=o['state'];p=o['peer_curvature'];C=own['curvature']
        assert own['last_k']==k-1
        bracket=bool(o['peer_bracket'] and own['prev_defect'] not in (None,0))
        mirror_C=p if bracket and p is not None else C
        # This choice is designed for a transient that lies between the most recent
        # target curvature and that of a paired residue class. No holdout used.
        sign=signs[(s,k%period)]
        raw_quad=sign*(2*own['last']-own['penult']+C)
        raw_mirror=sign*(2*own['last']-own['penult']+mirror_C)
        # Calibratable feature-only proxies for uncertainty, not probability bounds.
        defect_prev=own['prev_defect']
        tol = (abs(defect_prev) if defect_prev is not None else abs(C))
        if p is not None:tol += abs(C-p)
        out.append(dict(degree=d,shell=t,residue=s,k=k,period_phase=sign,
          baseline_quad_PT=str(raw_quad),predicted_PT=str(raw_mirror),
          model='mirror_bracket_curvature' if bracket else 'trailing_quadratic',
          observed_curvature=str(C),observed_previous_curvature=str(own['prev_curvature']) if own['prev_curvature'] is not None else '',
          mirror_curvature=str(p) if p is not None else '',
          previous_defect=str(defect_prev) if defect_prev is not None else '',
          uncertainty_proxy=str(tol),
          cert_stable_observed=int(defect_prev==0) if defect_prev is not None else 0))
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--calibration',default=str(ROOT/'data/calibration_Q1_Q20.csv'))
    p.add_argument('--update',default=str(ROOT/'data/update_Q21_Q25.csv'))
    p.add_argument('--outdir',default=str(ROOT/'results'))
    args=p.parse_args()
    train=read_csv(args.calibration);update=read_csv(args.update)
    sites1=read_csv(str(ROOT/'data/sites_Q21_Q25.csv'))
    sites2=read_csv(str(ROOT/'data/sites_Q26_Q30.csv'))
    p1,audit=train_predict(train,sites1)
    signs={(int(a),int(b)):int(v) for key,v in audit['signs'].items() for a,b in [key.split(',')]}
    p2=predict_second_block(train,update,sites2,signs,audit['period'])
    outdir=Path(args.outdir);outdir.mkdir(parents=True,exist_ok=True)
    write_csv(outdir/'block_Q21_Q25.csv',p1)
    write_csv(outdir/'block_Q26_Q30.csv',p2)
    # Train supervised event-risk baselines ONLY on block1 outcomes (Q21..25)
    # and deploy them using only available Q<=25 features to block2.
    E1=event_records(train,update,signs,audit['period'])
    d1=[dict(degree=e['degree'],actual_defect=str(e['actual_defect']),active=e['active'],
       **e['features']) for e in E1]
    write_csv(outdir/'event_training_Q21_Q25.csv',d1)
    audit['feature_engineering']='Signed second differences, third differences, p-adic valuations, paired residue curvature, neighboring shell, genus phase; all computed causally'
    audit['prediction_provenance']={'Q21_Q25':'Q<=20 only','Q26_Q30':'Q<=25 only'}
    audit['sha256']={str(fp.name):hashlib.sha256(fp.read_bytes()).hexdigest() for fp in (Path(args.calibration),Path(args.update))}
    (outdir/'model_audit.json').write_text(json.dumps(audit,indent=2))
    print('Predictions generated: Q21..25',len(p1),' Q26..30',len(p2),'without opening validation labels.')

if __name__=='__main__':main()
