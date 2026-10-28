#!/usr/bin/env python3
"""Creates exact empirical Z solely from PT integer observations + mature PT anchors.
No DT, modular, geometric operations imported into this transformation.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
import csv,sys
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE))
from empirical_coordinates import W,group,Z_from_PT,sigma

def read(path):return list(csv.DictReader(open(path)))
def write(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def run():
 anchors=read(BASE/'validation/mature_PT_anchor_observations.csv')
 scales={}
 for a in anchors:
  g,t=int(a['group']),int(a['depth']);d=int(a['degree']);s=d%5;k=d//5
  B=F(int(a['pt']),sigma(s,k)*W(s,t,k))
  assert (g,t) not in scales
  scales[g,t]=B
 out=[]
 for z in read(BASE/'validation/conditional_PT_observations.csv'):
  d,t=int(z['degree']),int(z['depth']);s=d%5;k=d//5;g=group(s)
  if (g,t) not in scales:continue
  Z=Z_from_PT(int(z['pt']),d,t,scales[g,t])
  out.append(dict(degree=d,depth=t,residue=s,r=int(z['r']),k=k,W=W(s,t,k),pt=z['pt'],Z=str(Z),scale_source_group=g,scale_source_depth=t))
 write(BASE/'validation/normalized_PT_observations.csv',out)
 print('Normalized PT samples',len(out),'excluded for missing mature anchor:',len(read(BASE/'validation/conditional_PT_observations.csv'))-len(out))
 return out
if __name__=='__main__':run()
