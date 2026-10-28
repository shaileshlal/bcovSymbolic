"""Separate no-oracle empirical predictor: reads Q<=20 labels, Q>=21 METADATA ONLY.

Does not open the published Q>=21 reference labels or the modular-series table.
The physical feature (-1)^chi*chi is made explicit; this is a restricted
arithmetic learner, not an unconstrained model discovery.
"""
import csv,hashlib,json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from math import floor

ROOT=Path(__file__).resolve().parent
CALIB=ROOT/'data'/'calibration_Q6_Q20.csv'
SITES=ROOT/'data'/'withheld_sites_Q21_Q40.csv'

def chargetuple(Q):
 s=Q%5
 d={0:0,1:1,2:2,3:2,4:1}[s]
 m= -(Q*(Q+5)//10)+d
 chi=m+Q-5
 sign=1 if chi%2==0 else -1
 return min(s,5-s),m,sign*chi

if not CALIB.exists() or not SITES.exists():
 raise FileNotFoundError('Run prepare_blind_splits.py once to create isolated inputs.')
train=list(csv.DictReader(CALIB.open()))
site=list(csv.DictReader(SITES.open()))
assert all(int(r['Q'])<=20 for r in train)
assert all(int(r['Q'])>=21 for r in site)
groups=defaultdict(list)
for r in train:
 Q=int(r['Q']);val=int(r['PT_integer'])
 mu,m,factor=chargetuple(Q)
 groups[mu].append(Fraction(val,factor))
learned={mu:sum(vals)/len(vals) for mu,vals in groups.items()}
assert all(v.denominator==1 for v in learned.values())
assert set(learned)=={0,1,2}

rows=[]
for r in site:
 Q=int(r['Q']);mu,m,mult=chargetuple(Q)
 prediction=learned[mu]*mult
 assert prediction.denominator==1
 rows.append({'Q':Q,'genus':1-m,'m':m,'mu':mu,'charge_factor':mult,'predicted_PT':int(prediction)})
output=ROOT/'results'/'isolated_predictions_Q21_Q40.csv'
output.parent.mkdir(exist_ok=True)
with output.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
manifest={'training_samples':len(train),'withheld_sites':len(site),
          'learned_small_integers':{str(k):str(v) for k,v in learned.items()},
          'used_inputs':{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in (CALIB,SITES)},
          'output_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
          'note':'No held-out integer labels read by this program.'}
(ROOT/'results'/'isolated_run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('wrote',output)
print(json.dumps(manifest,indent=2))
