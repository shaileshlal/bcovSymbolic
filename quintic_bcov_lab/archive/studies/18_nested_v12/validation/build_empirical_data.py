"""Validation-only construction of calibration labels from observed PT + mature PT scales.
All data files under data/ use only earlier-r or lower-degree observations.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
import csv,sys
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE))
from predict_terminal import write,terminal_t
z=defaultdict(dict);details=defaultdict(dict)
for row in csv.DictReader(open(BASE/'validation/normalized_PT_observations.csv')):
 s=int(row['residue']);r=int(row['r']);t=int(row['depth'])
 z[s,r][t]=F(row['Z']);details[s,r][t]=row

def eps(s,r,t):
 a=z[s,r]
 return a[t]-3*a[t-1]+3*a[t-2]-a[t-3] if all(j in a for j in range(t-3,t+1)) else None

calib={0:(6,7,8),1:(6,7,8),2:(4,5,6),3:(4,5),4:(5,6)}
held={0:(9,10,11),1:(9,10,11,12),2:(7,8,9),3:(6,7,8,9),4:(7,8,9,10)}
tr=[];hold=[];truth=[];rec=[]
for s,rs in calib.items():
 for r in rs:
  t=terminal_t(s,r);e=eps(s,r,t)
  assert e is not None and e.denominator==1 and e and eps(s,r,t+1)==0
  tr.append(dict(residue=s,r=r,t=t,terminal_pulse=str(e.numerator)))
for s,rs in held.items():
 for r in rs:
  t=terminal_t(s,r);e=eps(s,r,t)
  assert e is not None and e.denominator==1 and e and eps(s,r,t+1)==0
  hold.append(dict(residue=s,r=r,t=t))
  truth.append(dict(residue=s,r=r,t=t,terminal_pulse=str(e.numerator)))
  for u in range(t-3,t):
   rec.append(dict(residue=s,r=r,target_t=t,anchor_t=u,Z=str(z[s,r][u])))
write(BASE/'data/terminal_pulse_calibration.csv',tr)
write(BASE/'data/terminal_pulse_holdout_sites.csv',hold)
write(BASE/'data/terminal_pulse_reconstruction_anchors.csv',rec)
write(BASE/'validation/terminal_pulse_holdout_truth.csv',truth)

# Secondary pulse diagonal data: calibration r<=9, and a prespecified target grid r=10,11.
diag=defaultdict(dict)
for (s,r),seq in z.items():
 for t in seq:
  v=eps(s,r,t)
  if v is not None:diag[s,r][t]=v
calib=[];sites=[];targets=[]
for s in range(5):
 for ell in range(4):
  seq=[]
  for r in range(4,10):
   t=terminal_t(s,r)-ell;e=diag[s,r].get(t)
   if e is not None and e:seq.append(dict(residue=s,offset=ell,r=r,depth=t,epsilon=str(e)))
  if len(seq)<3:continue
  calib.extend(seq)
  for r in (10,11):
   t=terminal_t(s,r)-ell
   sites.append(dict(residue=s,offset=ell,r=r,depth=t))
   e=diag[s,r].get(t)
   if e is not None:targets.append(dict(residue=s,offset=ell,r=r,depth=t,epsilon=str(e)))
write(BASE/'data/diagonal_pulse_calibration.csv',calib)
write(BASE/'data/diagonal_holdout_sites.csv',sites)
write(BASE/'validation/diagonal_holdout_truth.csv',targets)

# First failure r6, train d<=53; test d>53; select t<=22 for reference quality.
r6tr=[];r6sites=[];r6truth=[]
for s in range(5):
 for t,v in sorted(z[s,6].items()):
  data=details[s,6][t];d=int(data['degree'])
  if t>22:continue
  if d<=53:r6tr.append(dict(residue=s,t=t,degree=d,W=data['W'],Z=str(v)))
  else:
   r6sites.append(dict(residue=s,t=t,degree=d,W=data['W']))
   r6truth.append(dict(residue=s,t=t,degree=d,W=data['W'],Z=str(v)))
write(BASE/'data/r6_Z_calibration.csv',r6tr)
write(BASE/'data/r6_Z_holdout_sites.csv',r6sites)
write(BASE/'validation/r6_Z_holdout_truth.csv',r6truth)
print('terminal',len(tr),'calibration',len(hold),'targets; diag',len(calib),'calibration',len(sites),'sites',len(targets),'scorable; r6',len(r6tr),'train',len(r6sites),'holdouts')
