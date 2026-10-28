"""Post-hoc *conditional* shape test at shells16..20. Uses exact reference B_t.
The B_t values are ORACLES here, so this is not a PT integer prediction.
Fits P,Q from data through depth 15 and validates on unseen deeper shells.
"""
from __future__ import annotations
import csv,json,sys
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from deep_envelope import W
from recursive_defects import read,phase_fit
from pulse_profile import group,K
from pt_from_dt_solver import pt_mC,simple_vacuum_pt
from extended_reference import evaluate,DELTA
per,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))

def B(c,t):
 d={0:50,1:49,2:48}[c];s=d%5;k=d//5
 return F(simple_vacuum_pt(d,pt_mC(d)+DELTA[s]+t),phase[s,k%per]*W(s,t,k))

def signed_profile(d,t,pt):
 s=d%5;k=d//5;r=K(s,t)-k;w=W(s,t,k)
 a=F(int(pt),phase[s,k%per]*w)
 return (-1 if t%2 else 1)*(a-B(group(s),t)),w,r

rows=[]
for fn in ('reference_Q20_Q53.csv','reference_deep_11_15.csv'):
 for row in csv.DictReader(open(ROOT/'validator'/fn)):
  if row['pt'] and row['status'] in ('certified_within_conditional_framework','certified_conditional'):
   rows.append(row)
profs={}
for s in range(5):
 for r in range(1,6):
  pts=[]
  for z in rows:
   d,t=int(z['degree']),int(z['shell'])
   if d%5!=s:continue
   y,w,rr=signed_profile(d,t,z['pt'])
   if rr==r:pts.append((w,y,t,d))
  distinct={w:y for w,y,t,d in pts}
  if len(distinct)<3:raise AssertionError(('insufficient calibration',s,r,len(distinct)))
  (w0,y0),(w1,y1)=sorted(distinct.items())[:2]
  Q=(y1-y0)/F(1,w1)-F(0) if False else (y1-y0)/(F(1,w1)-F(1,w0))
  P=y0-Q/F(w0)
  assert all(y==P+Q/F(w) for w,y in distinct.items()),('failed calibration',s,r)
  profs[s,r]=(P,Q)

hits=0;fail=[];by=Counter()
for d in range(55,111):
 for t in range(16,21):
  r=K(d%5,t)-d//5
  if not 1<=r<=5:continue
  z=evaluate(d,t)
  if z[3]!='certified_conditional':continue
  y,w,rr=signed_profile(d,t,z[4])
  P,Q=profs[d%5,r]
  if y!=P+Q/F(w):fail.append((d,t,r,str(y),str(P+Q/F(w))))
  else:hits+=1
  by[t]+=1
print('Future profile test depth16-20:',hits,'/',hits+len(fail),'by depth',dict(by),'fail sample',fail[:8])
(ROOT/'results/future_profile_test.json').write_text(json.dumps({'conditional_reference_only':True,'uses_exact_B_oracle':True,'training_depth_max':15,
  'test_depths':[16,17,18,19,20],'sites':hits+len(fail),'exact':hits,'errors':fail[:10],'by_depth':dict(by)},indent=2))
