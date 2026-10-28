"""Exploratory scan of residue periods for the PT/DT charge-normalized data.

The charge multiplier (-1)^chi chi is supplied as a kinematic feature;
the 5-fold residue structure itself is selected using ONLY training Q<=20.
"""
from collections import defaultdict
from fractions import Fraction as F
import json
from pathlib import Path
from empirical_benchmark import SOURCE,site

train=[(q,v) for q,v in SOURCE if q<=20]
test=[(q,v) for q,v in SOURCE if q>20]

# Rolling-origin tuning: for Q=16..20, train on all earlier degrees Q6..Q-1.
def predict(Q,earlier,period):
 ratios=defaultdict(list)
 global_ratios=[]
 for q,v in earlier:
  ratio=F(v,site(q)['kinematic_multiplier'])
  ratios[q%period].append(ratio)
  global_ratios.append(ratio)
 group=ratios[Q%period]
 learned=sum(group,F(0))/len(group) if group else sum(global_ratios,F(0))/len(global_ratios)
 return learned*site(Q)['kinematic_multiplier']

scores=[]
for period in range(1,16):
 valid=[]
 for Q,v in train:
  if Q<16:continue
  prior=[(q,u) for q,u in train if q<Q]
  pred=predict(Q,prior,period)
  valid.append((Q,pred==v,abs(pred-v)))
 out=[]
 for Q,v in test:
  pred=predict(Q,train,period)
  out.append((Q,pred==v,abs(pred-v)))
 scores.append({'period':period,'train_rolling_exact':sum(r[1] for r in valid),'train_rolling_total':len(valid),
                'holdout_exact':sum(r[1] for r in out),'holdout_total':len(out),
                'holdout_integer_abs_error_sum':str(sum((r[2] for r in out),F(0)))})

# p=5 is chosen as the smallest training-consistent period with perfect roll-forward.
assert next(a['period'] for a in scores if a['train_rolling_exact']==a['train_rolling_total'])==5
assert next(a for a in scores if a['period']==5)['holdout_exact']==10
Path(__file__).parent.joinpath('results/period_scan.json').write_text(json.dumps(scores,indent=2)+'\n')
for r in scores:print(r)
