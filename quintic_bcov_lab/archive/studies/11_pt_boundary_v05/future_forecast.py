#!/usr/bin/env python3
"""Unverified recurrence forecasts beyond the published Q30 table.
Use all published Q1..30 as calibration, discover where second differences
stabilize, and use *only* those series or validated mirror families.  No
modular form, wall crossing, or DT coefficient is used.
"""
import csv,json
from pathlib import Path
from fractions import Fraction
from run_experiments import lagrange_quad,nearest_int
from fit_empirical_transform import fit

def read(path):
 with open(path) as f:return list(csv.DictReader(f))
DELTA0={0:0,1:1,2:2,3:2,4:1}

def main():
 train=read('data/calibration.csv')+read('data/reference_only_Q21_Q30.csv')
 info=json.loads(Path('results/model_audit.json').read_text())
 d={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in train}
 def phase(q):return info['signs'][f'{q%5},{(q//5)%info["period"]}']
 def y(q,t):return phase(q)*d[(q,t)]
 families={}
 affine={}
 for s in range(5):
  qs=sorted([q for (q,t) in d if t==0 and q%5==s and q>=6])
  assert len(qs)==5
  for t in range(6-DELTA0[s]):
   vals=[y(q,t) for q in qs]
   diffs=[vals[i+2]-2*vals[i+1]+vals[i] for i in range(3)]
   stable=(diffs[-1]==diffs[-2])
   f=lagrange_quad([q//5 for q in qs[-3:]],vals[-3:])
   families[s,t]=dict(stable=stable,curvatures=diffs,polynomial=f)
   # Last 3 samples give a cross-residue affine slope with consistent residues.
   if t>0:
    u=[y(q,0) for q in qs[-3:]];v=vals[-3:]
    A=Fraction(v[1]-v[0],u[1]-u[0]);C=Fraction(v[0])-A*u[0]
    if all(A*x+C==z for x,z in zip(u,v)):affine[s,t]=(A,C)
 output=[]
 for q in range(31,41):
  s=q%5;k=q//5
  gc=(q*q+5*q+10)//10;mc=1-gc
  for t in range(6-DELTA0[s]):
   f=families[s,t]['polynomial'];pred=f(k);method='stabilized_second_difference' if families[s,t]['stable'] else 'no_stable_recurrence'
   if s in (1,2) and t>0:
    A,C=affine.get((5-s,t),(None,None))
    if A is not None:
     latest=sorted([z for z,tt in d if tt==0 and z%5==s])[-2:]
     if all(A*y(z,0)+C==y(z,t) for z in latest):
      pred=A*families[s,0]['polynomial'](k)+C;method='mirror_affine_recent_exact'
   output.append(dict(degree=q,genus=gc-DELTA0[s]-t,PT_m=mc+DELTA0[s]+t,residue=s,shell=t,
                      predicted_PT=str(nearest_int(phase(q)*pred)) if method!='no_stable_recurrence' else '',
                      method=method,verification_status='UNVERIFIED'))
 with Path('results/forward_Q31_Q40.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,output[0]);w.writeheader();w.writerows(output)
 good=[a for a in output if a['predicted_PT']!='']
 print('usable forecasts',len(good),'of',len(output))
 for a in output:
  if a['shell']==0 or (a['degree'] in (31,35,40) and a['shell']<=2):print(a)
if __name__=='__main__':main()
