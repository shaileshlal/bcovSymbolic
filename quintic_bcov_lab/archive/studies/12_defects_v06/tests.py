#!/usr/bin/env python3
"""Unit + no-reference-in-predictor + reproducibility tests."""
from __future__ import annotations
import csv,subprocess,tempfile,shutil,json,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def load(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 scores=json.loads((ROOT/'results/validated_metrics.json').read_text())
 assert scores['scores']['v05_fixed_train_Q1_Q20']['exact']==28
 assert scores['scores']['rolling_quad_second_block']['exact']==15
 assert scores['scores']['new_mirror_bracket_second_block']['exact']==17
 assert scores['scores']['v06_sequential_Q21_Q30']['exact']==31
 assert scores['defect_support']['active_first_window']==13
 assert scores['defect_support']['active_second_window']==9
 assert scores['defect_support']['new_activations']==0
 assert scores['defect_support']['stable_persistence']=={'success':11,'total':11}
 model=json.loads((ROOT/'results/event_metrics.json').read_text())
 assert abs(model['logit_positional']['accuracy']-.875)<1e-12
 assert abs(model['unit_speed_front']['accuracy']-23/24)<1e-12
 future=load(ROOT/'results/future_Q31_Q35.csv')
 leading={(int(x['degree']),int(x['shell'])):int(x['predicted_PT']) for x in future}
 assert leading[(31,0)]==-724500
 assert leading[(35,0)]==-550
 # Copy a self-contained source fixture which deliberately omits the held-out labels
 # entirely. Both prediction scripts must execute without such a file.
 with tempfile.TemporaryDirectory() as td:
  P=Path(td)
  for src in ('baseline_v05.py','defect_features.py','train_and_forecast.py','event_models.py'):
   shutil.copy(ROOT/src,P/src)
  (P/'data').mkdir();(P/'results').mkdir()
  for src in ('calibration_Q1_Q20.csv','update_Q21_Q25.csv','sites_Q21_Q25.csv','sites_Q26_Q30.csv'):
   shutil.copy(ROOT/'data'/src,P/'data'/src)
  assert not (P/'data/reference_only_Q26_Q30.csv').exists()
  subprocess.run([sys.executable,str(P/'train_and_forecast.py')],check=True,capture_output=True,text=True)
  subprocess.run([sys.executable,str(P/'event_models.py')],check=True,capture_output=True,text=True)
  for name in ('block_Q21_Q25.csv','block_Q26_Q30.csv','event_risk_predictions_Q26_Q30.csv'):
   assert digest(ROOT/'results'/name)==digest(P/'results'/name),name
 print('All exact, model, provenance and clean-without-holdout reproduction tests PASS.')

if __name__=='__main__':main()
