"""Score previously saved predictions against a separate conditional reference.
NOT imported by any predictor. Exact fractions; full integer equality is the metric.
"""
from __future__ import annotations
import csv,json,sys
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
from extended_reference import evaluate

def table(path):return list(csv.DictReader(open(path)))
ref={(int(z['degree']),int(z['shell'])):int(z['pt']) for z in table(ROOT/'validator/reference_deep_11_15.csv') if z['status']=='certified_conditional' and 41<=int(z['degree'])<=53}
assert len(ref)==65

results={}
for file,model in [('baseline_predictions.csv','linear'),('baseline_predictions.csv','quadratic'),('baseline_predictions.csv','cubic'),('baseline_predictions.csv','affine_envelope'),('baseline_predictions.csv','constant_envelope'),('linear_profile_pairs_predictions.csv','prediction')]:
 rows=table(ROOT/'results'/file)
 exact=0;available=0;integral=0;errs=[];by=defaultdict(lambda:dict(exact=0,predicted=0,abstained=0))
 for r in rows:
  d,t=int(r['degree']),int(r['shell'])
  raw=r.get(model,'')
  if not raw:
   by[t]['abstained']+=1
   continue
  available+=1;by[t]['predicted']+=1
  v=F(raw);integral+=int(v.denominator==1)
  same=v==ref[d,t];exact+=int(same);by[t]['exact']+=int(same)
  errs.append(float(abs(F(v-ref[d,t],ref[d,t]))))
 errs.sort()
 results[model]={'predicted':available,'exact':exact,'integral_candidates':integral,'abstained':65-available,
   'median_relative_error':errs[len(errs)//2] if errs else None,
   'max_relative_error':max(errs) if errs else None,'by_depth':dict(by)}

probe=table(ROOT/'results/prospective_54_70.csv')
by=defaultdict(lambda:dict(exact=0,predicted=0,abstained=0)); wrong=[]
for r in probe:
 d,t=int(r['degree']),int(r['shell']);p=r['prediction']
 if not p:by[t]['abstained']+=1;continue
 val=evaluate(d,t)
 assert val[3]=='certified_conditional'
 by[t]['predicted']+=1
 if int(p)==int(val[4]):by[t]['exact']+=1
 else:wrong.append((d,t,p,val[4]))
results['forward_54_70']={'predicted':sum(x['predicted'] for x in by.values()),'exact':sum(x['exact'] for x in by.values()),
 'abstained':sum(x['abstained'] for x in by.values()),'by_depth':dict(by),'wrong':wrong}
assert not wrong
(root:=ROOT/'results').mkdir(exist_ok=True)
(root/'benchmark.json').write_text(json.dumps(results,indent=2))
for k,v in results.items():print(k,v['exact'],'/',v['predicted'],'abstained',v['abstained'],'median_rel',v.get('median_relative_error'))
