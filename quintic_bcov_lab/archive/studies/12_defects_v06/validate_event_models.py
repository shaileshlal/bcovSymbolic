#!/usr/bin/env python3
"""Evaluation ONLY; executed after event_models.py has saved its forecasts."""
import json,csv
from pathlib import Path
from sklearn.metrics import accuracy_score,f1_score,precision_score,recall_score,roc_auc_score,brier_score_loss,balanced_accuracy_score
from baseline_v05 import read_csv
from defect_features import event_records
ROOT=Path(__file__).resolve().parent
train=read_csv(ROOT/'data/calibration_Q1_Q20.csv'); update=read_csv(ROOT/'data/update_Q21_Q25.csv');hold=read_csv(ROOT/'data/reference_only_Q26_Q30.csv')
info=json.loads((ROOT/'results/model_audit.json').read_text()); signs={(int(a),int(b)):int(v) for k,v in info['signs'].items() for a,b in [k.split(',')]}
labels={(e['degree'],e['shell']):e['active'] for e in event_records(train+update,hold,signs)}
pred=read_csv(ROOT/'results/event_risk_predictions_Q26_Q30.csv')
y=[labels[int(r['degree']),int(r['shell'])] for r in pred]
metrics={}
for name in [c[:-2] for c in pred[0] if c.endswith('_p')]:
 p=[float(r[name+'_p']) for r in pred]; yh=[int(v>=.5) for v in p]
 metrics[name]=dict(accuracy=accuracy_score(y,yh),balanced_accuracy=balanced_accuracy_score(y,yh),
  precision=precision_score(y,yh,zero_division=0),recall=recall_score(y,yh,zero_division=0),
  auc=roc_auc_score(y,p),brier=brier_score_loss(y,p),active_predicted=sum(yh),active_actual=sum(y))
(ROOT/'results/event_metrics.json').write_text(json.dumps(metrics,indent=2))
print(json.dumps(metrics,indent=2))
