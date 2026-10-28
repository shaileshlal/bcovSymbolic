"""Oracle-assisted finite-window rank audit (not counted as blind prediction)."""
import csv,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from linear_extension_pairs import rowsys,reduce
from pulse_profile import K
from recursive_defects import read,phase_fit
from deep_envelope import W

per,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))
prior=json.load(open(ROOT/'results/v10_fitted_pulses.json'))
train=[z for z in read(ROOT/'data/observed_calibration_Q20_Q40_shell6_15.csv') if K(int(z['degree'])%5,int(z['shell']))-int(z['degree'])//5<=5]
f=rowsys(train,phase,per,prior)
pivots={}
def add(co,b):
 co,b=reduce(co,b,pivots)
 if co:
  p=min(co);c=co[p];pivots[p]=({k:v/c for k,v in co.items()},b/c)
 else:assert not b, ('inconsistent profile',b)
for z in train:
 co,b,_,_=f(int(z['degree']),int(z['shell']),z['pt']);add(co,b)
calibrated=len(pivots)
anchors=[tuple(x) for x in json.load(open(ROOT/'results/anchor_rank.json'))['greedy_anchor_sites']]
ref={(int(z['degree']),int(z['shell'])):int(z['pt']) for z in csv.DictReader(open(ROOT/'validator/reference_deep_11_15.csv')) if z['status']=='certified_conditional'}
for d,t in anchors:
 co,b,_,_=f(d,t,str(ref[d,t]));add(co,b)
assert len(pivots)-calibrated==len(anchors)
checked=0;fails=[]
for d in range(41,54):
 for t in range(11,16):
  if K(d%5,t)-d//5>5:continue
  co,b,_,_=f(d,t)
  co,b=reduce(co,b,pivots)
  assert not co,(d,t,co)
  v=-b*phase[d%5,(d//5)%per]*W(d%5,t,d//5)
  if v!=ref[d,t]:fails.append((d,t,str(v),str(ref[d,t])))
  checked+=1
print('Oracle-assisted anchor check:', checked, 'exact matches',checked-len(fails),'anchors',len(anchors),'failures',fails)
assert not fails
(ROOT/'results/anchor_reconstruction_audit.json').write_text(json.dumps({'oracle_assisted':True,'anchor_count':len(anchors),'sites_reconstructed':checked,'exact':checked,'anchor_sites':anchors},indent=2))
