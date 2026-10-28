#!/usr/bin/env python3
"""Small deterministic tests of the finite-difference PT/GV experiment."""
import csv,json,subprocess,sys,tempfile
from pathlib import Path
from fractions import Fraction
from run_experiments import exact_phase_fit,train_predict
from fit_empirical_transform import fit
from finite_difference_fingerprint import fingerprint
ROOT=Path(__file__).resolve().parent

def read(p):
 with p.open() as f:return list(csv.DictReader(f))

def main():
 cal=read(ROOT/'data/calibration.csv');ref=read(ROOT/'data/reference_only_Q21_Q30.csv');sites=read(ROOT/'data/sites_unlabeled.csv')
 period,signs=exact_phase_fit(cal)
 assert period==4,period
 pred,info=train_predict(cal,sites)
 assert len(pred)==48
 T={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in ref}
 hits=[r for r in pred if int(r['predicted_PT'])==T[int(r['degree']),int(r['shell'])]]
 assert len(hits)==28,len(hits)
 assert sum(r['model']=='mirror_affine' for r in pred)==6
 assert all(int(r['predicted_PT'])==T[int(r['degree']),int(r['shell'])] for r in pred if r['model']=='mirror_affine')
 assert all((int(r['predicted_PT'])>0)==(T[int(r['degree']),int(r['shell'])]>0) for r in pred)
 known={(int(r['degree']),int(r['shell'])):r for r in ref}
 for d in range(21,31):
  for t in (0,1):
   row=next(r for r in pred if int(r['degree'])==d and int(r['shell'])==t)
   assert int(row['predicted_PT'])==T[d,t]
 r1,r2,diagnostics=fit(str(ROOT/'data/calibration.csv'))
 assert diagnostics['r1_support']=='14/15'
 assert diagnostics['r2_support']=='12/15'
 assert all(r1(h)==2*h-2 for h in range(3,9))
 assert all(r2(h)==(2*h-2)*(2*h-3)//2 for h in range(3,9))
 from collections import defaultdict
 PP={(int(r['degree']),int(r['shell'])):int(r['predicted_PT']) for r in pred}
 hits_gv=0
 for d in range(21,31):
  h=int(known[d,0]['gc'])-int(known[d,0]['delta'])
  g0=PP[d,0]
  g1=PP[d,1]-r1(h)*g0
  g2=PP[d,2]-r1(h-1)*g1-r2(h)*g0
  for t,val in enumerate((g0,g1,g2)):
   if val==int(known[d,t]['gv']):hits_gv+=1
 assert hits_gv==26,hits_gv
 F=fingerprint(cal,info)
 known_extract={(r['residue'],r['dt_index']):Fraction(r['train_based_DT_candidate']) for r in F}
 assert known_extract[0,0]==5 and known_extract[0,1]==-800
 assert known_extract[4,1]==8625 and known_extract[4,2]==-1138500
 assert known_extract[3,2]==-1218500 and known_extract[3,3]==441969250
 future=read(ROOT/'results/forward_Q31_Q40.csv')
 assert len(future)==48
 assert sum(bool(r['predicted_PT']) for r in future)==38
 assert next(int(r['predicted_PT']) for r in future if int(r['degree'])==35 and int(r['shell'])==0)==-550
 assert next(int(r['predicted_PT']) for r in future if int(r['degree'])==40 and int(r['shell'])==0)==725
 # Isolated predictor: copy only source + calibration + degree sites.
 with tempfile.TemporaryDirectory() as z:
  p=Path(z);(p/'data').mkdir();(p/'results').mkdir()
  import shutil
  for fn in ['run_experiments.py']:
   shutil.copyfile(ROOT/fn,p/fn)
  for fn in ['calibration.csv','sites_unlabeled.csv']:
   shutil.copyfile(ROOT/'data'/fn,p/'data'/fn)
  subprocess.run([sys.executable,'run_experiments.py','--train','data/calibration.csv',
                  '--sites','data/sites_unlabeled.csv','--out','results/pt_predictions.csv'],
                 cwd=z,check=True,capture_output=True)
  assert (p/'results/pt_predictions.csv').read_bytes()==(ROOT/'results/pt_predictions.csv').read_bytes()
 print('PASS: 4-period signs, 28/48 PT, 26/30 GV, 12 DT fingerprints, 38 forecasts, isolated predictor.')
if __name__=='__main__':main()
