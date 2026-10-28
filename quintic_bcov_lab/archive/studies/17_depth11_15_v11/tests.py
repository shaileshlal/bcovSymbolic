"""Reproducibility: exact reference, no-oracle predictor, profile, closure, rank checks."""
from __future__ import annotations
import csv,json,sys,subprocess,shutil,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def readcsv(p):
 with open(p) as f:return list(csv.DictReader(f))

def check():
 a=readcsv(ROOT/'validator/reference_Q20_Q53.csv')
 sys.path.insert(0,str(ROOT/'validator'))
 from extended_reference import evaluate
 prev=[z for z in a if z['status']=='certified_within_conditional_framework']
 assert len(prev)==116
 assert all(evaluate(int(z['degree']),int(z['shell']))[4]==z['pt'] for z in prev)
 reference={(int(z['degree']),int(z['shell'])):z['pt'] for z in readcsv(ROOT/'validator/reference_deep_11_15.csv') if z['status']=='certified_conditional'}
 p=readcsv(ROOT/'results/linear_profile_pairs_predictions.csv')
 got=[z for z in p if z['prediction']]
 assert len(p)==65 and len(got)==25
 assert all(z['prediction']==reference[int(z['degree']),int(z['shell'])] for z in got)
 q=readcsv(ROOT/'results/prospective_54_70.csv')
 assert len(q)==85 and sum(bool(z['prediction']) for z in q)==49
 assert all(evaluate(int(z['degree']),int(z['shell']))[4]==z['prediction'] for z in q if z['prediction'])
 a=json.loads((ROOT/'results/anchor_rank.json').read_text());assert a['new_anchor_count']==8
 b=json.loads((ROOT/'results/anchor_reconstruction_audit.json').read_text());assert b['sites_reconstructed']==46 and b['exact']==46
 c=json.loads((ROOT/'results/future_profile_test.json').read_text());assert c['exact']==123 and c['sites']==123 and c['uses_exact_B_oracle']
 d=json.loads((ROOT/'results/closure_tests.json').read_text())
 assert all(x['hankel_ranks']['8']==8 for x in d.values())
 assert all(not any(z['null_dim'] for z in x['recurrence_tests']) for x in d.values())
 # No-oracle test: runnable predictor in an isolated directory with NO validator present.
 with tempfile.TemporaryDirectory() as td:
  root=Path(td);shutil.copytree(ROOT/'predictor',root/'predictor',ignore=shutil.ignore_patterns('__pycache__'))
  (root/'data').mkdir();(root/'results').mkdir()
  for file in ['shallow_calibration_Q1_Q20.csv','observed_calibration_Q20_Q40_shell6_15.csv','withheld_sites_Q41_Q53_shell11_15.csv']:
   shutil.copy2(ROOT/'data'/file,root/'data'/file)
  shutil.copy2(ROOT/'results/v10_fitted_pulses.json',root/'results/v10_fitted_pulses.json')
  assert not (root/'validator').exists()
  for name in ('linear_extension_pairs','prospective','anchor_rank'):
   subprocess.run([sys.executable,str(root/'predictor'/f'{name}.py')],cwd=root,check=True,stdout=subprocess.DEVNULL)
  for fn in ('linear_profile_pairs_predictions.csv','prospective_54_70.csv','anchor_rank.json'):
   assert (ROOT/'results'/fn).read_bytes()==(root/'results'/fn).read_bytes(),fn
 print('PASS: 116 cross-reference checks; 25+49 exact predictions; 8-anchor rank; 123 conditional profile checks; no low-order closure; isolated no-oracle rerun')
if __name__=='__main__':check()
