"""Independent scoring. Model has no access to these reference files."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import csv,json
BASE=Path(__file__).resolve().parent.parent

def read(p):return list(csv.DictReader(open(p)))
def idx(a,keys):return {tuple(v[k] for k in keys):v for v in a}

def run():
 res={}
 p=read(BASE/'results/terminal_pulse_predictions.csv');v=idx(read(BASE/'validation/terminal_pulse_holdout_truth.csv'),('residue','r','t'))
 c=sum(F(a['terminal_pulse'])==F(v[a['residue'],a['r'],a['t']]['terminal_pulse']) for a in p)
 res['terminal_pulses']={'predicted':len(p),'exact':c}
 assert c==len(p)==18

 p=read(BASE/'results/diagonal_pulse_predictions.csv');v=idx(read(BASE/'validation/diagonal_holdout_truth.csv'),('residue','offset','r'))
 valid=[];missing=[]
 for a in p:
  key=(a['residue'],a['offset'],a['r'])
  if key not in v:missing.append(key)
  else:valid.append((key,F(a['pulse_prediction']),F(v[key]['epsilon'])))
 c=sum(x==y for key,x,y in valid)
 res['diagonal_pulses']={'attempted':len(p),'scored':len(valid),'exact':c,'unscored':len(missing),
   'failures':[{'site':key,'error':str(x-y),'actual':str(y)} for key,x,y in valid if x!=y]}
 assert len(p)==24 and len(valid)==19 and c==17

 p=read(BASE/'results/r6_Z_predictions.csv');v=idx(read(BASE/'validation/r6_Z_holdout_truth.csv'),('residue','t'))
 c0=sum(F(a['earliest_affine_Z'])==F(v[a['residue'],a['t']]['Z']) for a in p)
 c1=sum(F(a['latest_quadratic_Z'])==F(v[a['residue'],a['t']]['Z']) for a in p)
 c2=sum(F(a['latest_affine_Z'])==F(v[a['residue'],a['t']]['Z']) for a in p)
 res['r6_normalized']={'scored':len(p),'old_affine_exact':c0,'latest_two_affine_exact':c2,'latest_three_exact':c1,
    'nonzero_corrections':sum(F(a['integer_correction'])!=0 for a in p)}
 assert len(p)==32 and c0==0 and c1==32 and c2==32

 # Reconstruct exactly the terminal coefficient after three given prior normalized values.
 p=read(BASE/'results/terminal_Z_reconstruction.csv')
 atlas=idx(read(BASE/'validation/normalized_PT_observations.csv'),('residue','r','depth'))
 c=sum(F(a['predicted_Z'])==F(atlas[a['residue'],a['r'],a['depth']]['Z']) for a in p)
 res['terminal_Z']={'scored':len(p),'exact':c};assert len(p)==18 and c==18
 (BASE/'results/benchmark.json').write_text(json.dumps(res,indent=2)+'\n')
 print(json.dumps(res,indent=2))
if __name__=='__main__':run()
