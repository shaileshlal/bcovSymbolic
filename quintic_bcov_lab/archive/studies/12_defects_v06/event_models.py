#!/usr/bin/env python3
"""Fits simple supervised models for presence of *next* third-difference defect.

Training labels from Q21..25 only; test features from Q<=25 only.  This script
NEVER opens Q26..30 reference labels.  The separately stored scored holdout is
only used by validate_event_models.py.
"""
from __future__ import annotations
import csv,math,sys,json
from pathlib import Path
from defect_features import grouped_signed_rows,feat_from_state
from baseline_v05 import read_csv

ROOT=Path(__file__).resolve().parent
BASE=['shell','t_minus_k']
ARITH=BASE+['absolute_C_log','absolute_PT_log','relative_curvature_log',
  'local_ratio_log','v2_C','v5_C','v5_PT','peer_gap_ratio','peer_gap_log',
  'prior_defect_ratio','curvature_stable','mirror_bracket']

def get_input():
    model=json.loads((ROOT/'results/model_audit.json').read_text())
    signs={(int(s),int(k)):int(v) for key,v in model['signs'].items() for s,k in [key.split(',')]}
    train=read_csv(ROOT/'results/event_training_Q21_Q25.csv')
    data=read_csv(ROOT/'data/calibration_Q1_Q20.csv')+read_csv(ROOT/'data/update_Q21_Q25.csv')
    sites=read_csv(ROOT/'data/sites_Q26_Q30.csv')
    groups=grouped_signed_rows(data,signs,model['period'])
    test=[]
    for p in sites:
        d=int(p['degree']);s=d%5;t=int(p['shell']);k=d//5
        features=feat_from_state(groups,s,t)
        assert features['state']['last_k']==k-1
        test.append(dict(degree=d,**features['features']))
    return train,test

def mat(rows,fields):
    import numpy as np
    return np.array([[float(r.get(x,0)) for x in fields]+[float(int(r['residue'])==s) for s in range(5)] for r in rows],dtype=float)

def write(rows,path):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)

def train_models():
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import ExtraTreesClassifier,RandomForestClassifier
    from sklearn.pipeline import make_pipeline
    train,test=get_input()
    y=np.array([int(row['active']) for row in train])
    predictions={(int(a['degree']),int(a['shell'])):dict(degree=a['degree'],residue=a['residue'],shell=a['shell']) for a in test}
    # No holdout hyperparameter tuning; fixed, low-complexity models.
    models={
       'logit_positional':(BASE,make_pipeline(StandardScaler(),LogisticRegression(C=.3,max_iter=2000,random_state=0))),
       'logit_engineered':(ARITH,make_pipeline(StandardScaler(),LogisticRegression(C=.3,max_iter=2000,random_state=0))),
       'extratrees_engineered':(ARITH,ExtraTreesClassifier(n_estimators=200,max_depth=3,min_samples_leaf=2,random_state=0)),
       'randomforest_engineered':(ARITH,RandomForestClassifier(n_estimators=200,max_depth=3,min_samples_leaf=2,random_state=0)),
    }
    for name,(fields,clf) in models.items():
        clf.fit(mat(train,fields),y)
        p=clf.predict_proba(mat(test,fields))[:,1]
        for r,prob in zip(test,p):
            predictions[(int(r['degree']),int(r['shell']))][name+'_p']=round(float(prob),8)
    # Two event-front hypotheses; distinct from predicting the amplitude.
    first_front={}
    for r in train:
        if int(r['active']) and int(r['residue']) not in first_front:
            first_front[int(r['residue'])]=int(r['shell'])
        elif int(r['active']):first_front[int(r['residue'])]=min(first_front[int(r['residue'])],int(r['shell']))
    for r in test:
        s=int(r['residue']);t=int(r['shell']);key=(int(r['degree']),t)
        predictions[key]['persistent_front_p']=int(t>=first_front[s])
        predictions[key]['unit_speed_front_p']=int(t>=first_front[s]+1)
        predictions[key]['stable_absorbing_p']=int(not bool(r['curvature_stable']))
        predictions[key]['first_front']=first_front[s]
    pred=list(predictions.values()); pred.sort(key=lambda r:(int(r['degree']),int(r['shell'])))
    write(pred,ROOT/'results/event_risk_predictions_Q26_Q30.csv')
    (ROOT/'results/event_model_manifest.json').write_text(json.dumps({'train_degrees':'<=25 (24 target events)', 'test_features':'<=25 only','reference_labels_read':False,
       'features':{'positional':BASE,'engineered':ARITH}, 'fronts_from_training':first_front},indent=2))
    print('Event classifiers trained on',len(train),'defect labels; predicted',len(pred),'risk states, no holdout target read')

if __name__=='__main__':train_models()
