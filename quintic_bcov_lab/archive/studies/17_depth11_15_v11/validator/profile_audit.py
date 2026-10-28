"""Post-hoc examination of empirical profile complexity. Uses validation-only DT
vacuum scale to subtract unknown base, NOT a predictor or independent data.
"""
import sys,csv,json
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from deep_envelope import W
from recursive_defects import read,phase_fit
from pulse_profile import K,group
from pt_from_dt_solver import pt_mC,simple_vacuum_pt
DELTA={0:0,1:1,2:2,3:2,4:1}
period,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))

def Bclass(g,t):
 d={0:50,1:49,2:48}[g];s=d%5;k=d//5
 m=pt_mC(d)+DELTA[s]+t
 return F(simple_vacuum_pt(d,m),phase[s,k%period]*W(s,t,k))

def lagrange(pts,z):
 out=F(0)
 for i,(x,y) in enumerate(pts):
  val=y
  for j,(xx,yy) in enumerate(pts):
   if i!=j:val*=F(z-xx,x-xx)
  out+=val
 return out

ref=[]
for fn in ('reference_Q20_Q53.csv','reference_deep_11_15.csv'):
 for r in csv.DictReader(open(ROOT/'validator'/fn)):
  if r.get('pt') and r['status'] in ('certified_within_conditional_framework','certified_conditional'):
   ref.append(r)
series=defaultdict(list)
for row in ref:
 d,t=int(row['degree']),int(row['shell']);s=d%5;k=d//5;r=K(s,t)-k
 if r<=0:continue
 w=W(s,t,k); A=F(int(row['pt']),phase[s,k%period]*w)
 y=(-1 if t%2 else 1)*(A-Bclass(group(s),t))
 series[s,r].append((F(1,w),y,t,d,w))

out=[]
for (s,r),items in sorted(series.items()):
 # distinct W only
 uniq={z:y for z,y,t,d,w in items}
 if len(uniq)<2:continue
 points=sorted(uniq.items())
 fits=[]
 for deg in range(0,min(6,len(points))):
  xs=points[:deg+1]
  if all(lagrange(xs,z)==y for z,y in points):
   fits.append(deg);break
 dmin=fits[0] if fits else None
 out.append({'s':s,'r':r,'points':len(items),'unique_W':len(points),'min_degree':dmin,
 'overdetermined':dmin is not None and len(points)>=dmin+2,
 'max_depth':max(z[2] for z in items),'min_depth':min(z[2] for z in items)})
print('profiles with >=2 unique W',len(out))
for r in range(1,11):
 st=[(z['s'],z['unique_W'],z['min_degree'],z['overdetermined']) for z in out if z['r']==r]
 if st: print('r',r,st)
print('overdetermined certified by degree',Counter(z['min_degree'] for z in out if z['overdetermined']))
print('nonzero min degree at r>=4', [z for z in out if z['r']>=4 and z['overdetermined']][:10])
(ROOT/'results'/'posthoc_profile_audit.json').write_text(json.dumps(out,indent=2))
