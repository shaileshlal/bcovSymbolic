"""True one-step prospective candidate forecasts for Q31..35, no reference labels."""
from pathlib import Path
import csv,json
from train_and_forecast import predict_second_block,write_csv
from baseline_v05 import read_csv
from defect_features import DELTA0
ROOT=Path(__file__).resolve().parent
train=read_csv(ROOT/'data/calibration_Q1_Q20.csv')
update=read_csv(ROOT/'data/update_Q21_Q25.csv')+read_csv(ROOT/'data/reference_only_Q26_Q30.csv')
audit=json.loads((ROOT/'results/model_audit.json').read_text())
signs={(int(a),int(b)):int(v) for k,v in audit['signs'].items() for a,b in [k.split(',')]}
sites=[dict(degree=str(q),shell=str(t)) for q in range(31,36) for t in range(6-DELTA0[q%5])]
preds=predict_second_block(train,update,sites,signs,audit['period'])
for z in preds:z['verification_status']='UNVERIFIED_OUTSIDE_TABLE'
write_csv(ROOT/'results/future_Q31_Q35.csv',preds)
for z in preds:
 if int(z['shell'])<=2:
  print(z['degree'],z['shell'],z['predicted_PT'],z['model'])
