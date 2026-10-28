#!/usr/bin/env python3
"""PT second-difference fingerprint, NO DT input. Compare only in validator.

For each of three calibration residue families s=0,4,3, strip the period-4
sign learned on Q<=20 and measure the second difference across k.
The factor of 5 and alternating sign are an empirical normalization, chosen
from the quintic's arithmetic scale, not imposed via modularity.
"""
import csv,json
from fractions import Fraction
from pathlib import Path

def read(path):
 with open(path) as f:return list(csv.DictReader(f))

FAMILIES=[(0,0,6,1),(4,1,5,1),(3,2,4,-1)]
# tuple(s,starting q-series index,number shell rows,alternating sign at t=0)

def phase(d,info):
 return info['signs'][f'{d%5},{(d//5)%info["period"]}']

def fingerprint(rows,info):
 data={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in rows}
 result=[]
 for s,q0,ns,orient in FAMILIES:
  degrees=sorted([q for q,t in data.keys() if t==0 and q%5==s and q>=6])
  assert len(degrees)>=3
  for t in range(ns):
   if not all((d,t) in data for d in degrees):continue
   vals=[phase(d,info)*data[(d,t)] for d in degrees]
   diffs=[vals[i+2]-2*vals[i+1]+vals[i] for i in range(len(vals)-2)]
   factor=orient*((-1)**t)
   est=Fraction(factor*diffs[0],5)
   stable=(len(diffs)>=3 and diffs[-1]==diffs[-2])
   tail=Fraction(factor*diffs[-1],5) if stable else None
   result.append(dict(residue=s,shell=t,dt_index=q0+t,
                      training_curvature=str(diffs[0]),train_based_DT_candidate=str(est),
                      later_stabilized_curvature=str(diffs[-1]) if stable else '',
                      postholdout_DT_candidate=str(tail) if tail is not None else '',
                      curvature_windows=str(diffs),late_curvature_stabilized=stable))
 return result

def main():
 info=json.loads(Path('results/model_audit.json').read_text())
 train=read('data/calibration.csv')
 T=fingerprint(train,info)
 with open('results/calibration_dt_fingerprint.csv','w',newline='') as f:
  w=csv.DictWriter(f,T[0]);w.writeheader();w.writerows(T)
 # WARNING: validation against holdout *only after saving train fingerprint*.
 # This file deliberately does not use comparison labels to make predictions.
 print('TRAINING-ONLY fingerprint:')
 for x in T:
  print(x['residue'],x['dt_index'],x['train_based_DT_candidate'])

if __name__=='__main__':main()
