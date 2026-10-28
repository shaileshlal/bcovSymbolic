"""Degree-ordered ML baselines; no Q>=21 label enters fit.

ML models are deliberately small due to only fifteen independent training rows.
"""
from math import floor
import json, csv
from pathlib import Path
import numpy as np
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, RandomForestClassifier
from sklearn.neighbors import KNeighborsRegressor
from empirical_benchmark import site,SOURCE

train=[(q,v) for q,v in SOURCE if q<=20]
test=[(q,v) for q,v in SOURCE if q>20]

def features(q):
 s=site(q)
 return [q,q%5,q//5,s['genus'],s['chi']]
X=np.array([features(q) for q,_ in train], dtype=float)
Y=np.array([v for _,v in train],dtype=float)
XT=np.array([features(q) for q,_ in test],dtype=float)
actual=np.array([v for _,v in test],dtype=float)

preds={}
preds['Minimum-height zero']=[0]*len(test)
preds['Nearest same-residue']=[next(v for Q,v in reversed(train) if Q%5==q%5) for q,_ in test]
# fit signed degree-quadratic separately in 5 residue classes
sq=[]
for q,_ in test:
 s=q%5
 prev=[(Q,v) for Q,v in train if Q%5==s]
 assert len(prev)==3
 k=(q-prev[0][0])//5
 a,b,c=[v for _,v in prev]
 sq.append(a+k*(b-a)+(k*(k-1)//2)*(c-2*b+a))
preds['Signed quadratic by residue']=sq

for name,model in [
 ('Random forest',RandomForestRegressor(n_estimators=300,max_features=1.0,random_state=42)),
 ('ExtraTrees',ExtraTreesRegressor(n_estimators=300,random_state=42)),
 ('KNN-3',KNeighborsRegressor(n_neighbors=3))]:
 model.fit(X,Y)
 preds[name]=[int(round(x)) for x in model.predict(XT)]

# Physics-informed multiplicative DT-charge features, estimated from first 5 sites only.
first5={s['rankzero_mu']:v//s['kinematic_multiplier'] for q,v in train[:5] for s in [site(q)]}
preds['DT-charge normalization']=[site(q)['kinematic_multiplier']*first5[site(q)['rankzero_mu']] for q,_ in test]

scores={}
for name,p in preds.items():
 errs=[abs(a-b) for a,b in zip(p,actual)]
 scores[name]={'exact':sum(not e for e in errs),'count':len(errs),
               'MAE':float(sum(errs)/len(errs)),
               'median_abs_error':float(np.median(errs)),
               'max_abs_error':float(max(errs))}

results={'training_degrees':'6-20','holdout_degrees':'21-30','scores':scores,
         'warning':'small dataset, feature family was selected retrospectively; published Q21-30 values are comparison-only'}
root=Path(__file__).parent.joinpath('results')
root.joinpath('ablation.json').write_text(json.dumps(results,indent=2)+'\n')
with root.joinpath('model_predictions.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['degree','target_PT']+list(preds))
    for i,(q,target) in enumerate(test):w.writerow([q,target]+[preds[name][i] for name in preds])
print(json.dumps(scores,indent=2))
