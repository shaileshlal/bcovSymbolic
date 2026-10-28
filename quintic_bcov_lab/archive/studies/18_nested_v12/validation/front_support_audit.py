"""Post hoc empirical front-support audit from PT observations; no conjecture proof."""
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import csv,json,sys
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE))
from empirical_coordinates import terminal_t
z=defaultdict(dict)
for a in csv.DictReader(open(BASE/'validation/normalized_PT_observations.csv')):
 z[int(a['residue']),int(a['r'])][int(a['depth'])]=F(a['Z'])
counts=Counter();errors=[]
for (s,r),seq in z.items():
 if r<9:continue
 for t in seq:
  if not all(u in seq for u in range(t-3,t+1)):continue
  e=seq[t]-3*seq[t-1]+3*seq[t-2]-seq[t-3]
  if t==terminal_t(s,r):
   counts['terminal_scored']+=1;counts['terminal_active']+=(e!=0)
   if e==0:errors.append((s,r,t,'unexpected_zero_at_terminal'))
  if t>terminal_t(s,r):
   counts['postterminal_scored']+=1;counts['postterminal_zero']+=(e==0)
   if e!=0:errors.append((s,r,t,str(e)))
assert not errors and counts['terminal_scored']==13 and counts['postterminal_scored']==33
(BASE/'results/front_support_audit.json').write_text(json.dumps({'counts':dict(counts),'anomalies':errors},indent=2)+'\n')
print(dict(counts))
