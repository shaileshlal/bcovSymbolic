from pathlib import Path
import csv
from observations_raw import parse,PT,GV
D=Path(__file__).resolve().parent
GV_D=parse(GV);PT_D=parse(PT)
delta0={0:0,1:1,2:2,3:2,4:1}
cal=[];cmp=[];sites=[]
for d in range(1,31):
 for t in range(6-delta0[d%5]):
  delta=delta0[d%5]+t
  pt=PT_D[d]['values'][delta];gv=GV_D[d]['values'][delta]
  if pt is None or gv is None: continue
  rr=dict(degree=d,residue=d%5,k=d//5,shell=t,delta=delta,gc=GV_D[d]['base'],mc=PT_D[d]['base'],pt=pt,gv=gv)
  if d<=20:cal.append(rr)
  else:
   cmp.append(rr)
   sites.append(dict(degree=d,shell=t))
(D/'data').mkdir(exist_ok=True)
for path,vals in [(D/'data/calibration.csv',cal),(D/'data/reference_only_Q21_Q30.csv',cmp),(D/'data/sites_unlabeled.csv',sites)]:
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,vals[0]);w.writeheader();w.writerows(vals)
print('train',len(cal),'test',len(cmp))
