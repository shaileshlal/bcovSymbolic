"""Count independent empirical integer anchors by exact rational row rank.
No target values are read: only unlabelled (degree,shell) indices.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from linear_extension_pairs import rowsys,reduce
from recursive_defects import phase_fit,read
from pulse_profile import K

def run():
 per,ph=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))
 prior=json.load(open(ROOT/'results/v10_fitted_pulses.json'))
 train=[z for z in read(ROOT/'data/observed_calibration_Q20_Q40_shell6_15.csv') if K(int(z['degree'])%5,int(z['shell']))-int(z['degree'])//5<=5]
 f=rowsys(train,ph,per,prior)
 pivots={};bad=[]
 for item in train:
  co,b,_,_=f(int(item['degree']),int(item['shell']),item['pt'])
  co,b=reduce(co,b,pivots)
  if co:
   p=min(co); c=co[p];pivots[p]=({k:v/c for k,v in co.items()},b/c)
  elif b:bad.append(item)
 assert not bad
 initial=len(pivots)
 stats={}; anchors=[]
 for t in range(11,16):
  sites=[d for d in range(41,54) if K(d%5,t)-d//5<=5]
  inc=0
  for d in sites:
   co,b,_,_=f(d,t)
   co,b=reduce(co,b,pivots)
   if co:
    p=min(co);c=co[p]
    pivots[p]=({k:v/c for k,v in co.items()},b/c)
    inc+=1;anchors.append((d,t))
  stats[t]={'eligible':len(sites),'new_independent_anchors':inc}
 print('training rank',initial,'hypothetical target ranks',stats,'total new anchors',len(anchors),'chosen anchors',anchors)
 (ROOT/'results/anchor_rank.json').write_text(json.dumps({'training_rank':initial,'by_depth':stats,'new_anchor_count':len(anchors),'greedy_anchor_sites':anchors},indent=2))
if __name__=='__main__':run()
