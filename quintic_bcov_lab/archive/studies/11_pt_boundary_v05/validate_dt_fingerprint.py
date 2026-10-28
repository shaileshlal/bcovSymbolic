#!/usr/bin/env python3
"""Only validation: compare PT finite differences to independently listed DT coefficients.
The latter are reproduced from the publication, never passed to the predictor.
The reference PT tables themselves were generated using related DT information;
this is structural correspondence, NOT independent confirmation of modularity.
"""
from pathlib import Path
import csv,json
from fractions import Fraction
from finite_difference_fingerprint import read,fingerprint
DT_REFERENCE = {
(0,0):5,(0,1):-800,(0,2):58500,(0,3):5817125,(0,4):75474060100,
(4,1):8625,(4,2):-1138500,(4,3):3777474000,(4,4):3102750380125,
(3,2):-1218500,(3,3):441969250,(3,4):953712511250}

def main():
 train=read('data/calibration.csv');withheld=read('data/reference_only_Q21_Q30.csv')
 info=json.loads(Path('results/model_audit.json').read_text())
 calc0=fingerprint(train,info)
 combined=fingerprint(train+withheld,info)
 out=[]
 for a in combined:
  k=(a['residue'],a['dt_index']);truth=DT_REFERENCE.get(k)
  if truth is None:continue
  b=next((r for r in calc0 if (r['residue'],r['dt_index'])==k),None)
  low=int(Fraction(b['train_based_DT_candidate'])) if b is not None else None
  late=int(Fraction(a['postholdout_DT_candidate'])) if a['postholdout_DT_candidate'] else None
  out.append(dict(residue=k[0],q_index=k[1],train_PT_curvature_estimate=low,
                  published_DT_reference=truth,train_error=low-truth if low is not None else '',
                  train_relative_error=abs(low-truth)/abs(truth) if low is not None else '',
                  stabilized_with_extra_PT_data=late if late is not None else '',
                  late_exact=(late==truth) if late is not None else ''))
 with Path('results/dt_fingerprint_comparison.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,out[0]);w.writeheader();w.writerows(out)
 print('PT-second-difference to DT comparison:')
 for a in out: print(a)
 print('train exact',sum(a['train_error']==0 for a in out if a['train_error']!=''),
       'late exact',sum(a['late_exact']==True for a in out if a['late_exact']!=''))
if __name__=='__main__':main()
