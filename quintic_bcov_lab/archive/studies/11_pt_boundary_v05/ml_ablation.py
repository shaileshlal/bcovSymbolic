#!/usr/bin/env python3
"""Non-arithmetic baselines on the same degree split, no labels above degree 20 in fitting.
Score logarithmic magnitudes, not exact integer reconstruction.
"""
import csv,json,math
from pathlib import Path
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor,RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler,PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from run_experiments import lagrange_quad

def read(path):
 with open(path) as f:return list(csv.DictReader(f))

def feat(r):
 d=int(r['degree']);s=d%5;k=d//5;t=int(r['shell']);delta=int(r['delta']);g=int(r['gc'])-delta
 return [k,t,delta,g,d,(g-1)/max(1,d),int(s==0),int(s==1),int(s==2),int(s==3),int(s==4),
         k*t,k*k,t*t, k*(s+1)]

def main():
 a=[r for r in read('data/calibration.csv') if int(r['degree'])>=6]
 b=read('data/reference_only_Q21_Q30.csv')
 X=np.array([feat(r) for r in a]);Y=np.array([math.log10(abs(int(r['pt']))) for r in a]);
 V=np.array([feat(r) for r in b]);Truth=np.array([math.log10(abs(int(r['pt']))) for r in b]);
 m={
  'RandomForest':RandomForestRegressor(n_estimators=300,random_state=17,min_samples_leaf=1,n_jobs=1),
  'ExtraTrees':ExtraTreesRegressor(n_estimators=300,random_state=17,min_samples_leaf=1,n_jobs=1),
  'KNN5':make_pipeline(StandardScaler(),KNeighborsRegressor(n_neighbors=5)),
  'Ridge-polynomial':make_pipeline(StandardScaler(),PolynomialFeatures(2,include_bias=False),Ridge(alpha=100)),
 }
 res={}
 for name,model in m.items():
  model.fit(X,Y);p=model.predict(V)
  res[name]={'median_abs_dex':float(np.median(abs(p-Truth))),'mean_abs_dex':float(np.mean(abs(p-Truth))),
             'median_multiplicative_factor':float(10**np.median(abs(p-Truth))),
             'exact_integer_predictions':0} # rounding floats at 1e15 cannot certify exact; no matches expected
 # Compare with exact recurrence predictor.
 pred=read('results/pt_predictions.csv');pred={(int(r['degree']),int(r['shell'])):int(r['predicted_PT']) for r in pred}
 ep=np.array([math.log10(abs(pred[int(r['degree']),int(r['shell'])])) for r in b]);
 res['Periodic quadratic+mirror']={'median_abs_dex':float(np.median(abs(ep-Truth))),
     'mean_abs_dex':float(np.mean(abs(ep-Truth))),
     'median_multiplicative_factor':float(10**np.median(abs(ep-Truth))),
     'exact_integer_predictions':sum(pred[int(r['degree']),int(r['shell'])]==int(r['pt']) for r in b)}
 Path('results/ml_ablation.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
if __name__=='__main__':main()
