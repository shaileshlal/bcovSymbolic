#!/usr/bin/env python3
"""Validator-only generation of a PT observations table; predictor sees no DT code.
Exact PT integers are conditional on the validator's modular/wall assumptions.
"""
from pathlib import Path
import sys,csv
BASE=Path(__file__).resolve().parent.parent
SOURCE=BASE/'validation'/'reference_backend'
sys.path.insert(0,str(SOURCE))
sys.path.insert(0,str(BASE))
from extended_reference import evaluate
from empirical_coordinates import K,W,group

def write(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def main():
 reference=[]
 for s in range(5):
  for r in range(3,14):
   for t in range(8,29):
    k=K(s,t)-r;d=5*k+s
    if d<20:continue
    z=evaluate(d,t)
    if z[3]!='certified_conditional':continue
    reference.append(dict(degree=d,depth=t,residue=s,r=r,k=k,W=W(s,t,k),pt=z[4]))
 anchors=[]
 for g,s in ((0,0),(1,4),(2,3)):
  for t in range(8,29):
   k=K(s,t);d=5*k+s
   z=evaluate(d,t)
   if z[3]!='certified_conditional':continue
   anchors.append(dict(group=g,degree=d,depth=t,residue=s,k=k,W=W(s,t,k),pt=z[4]))
 write(BASE/'validation/conditional_PT_observations.csv',reference)
 write(BASE/'validation/mature_PT_anchor_observations.csv',anchors)
 print('PT observations',len(reference),'PT mature anchors',len(anchors))
if __name__=='__main__':main()
