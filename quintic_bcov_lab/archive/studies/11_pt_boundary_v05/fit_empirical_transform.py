#!/usr/bin/env python3
"""Learn a local PT -> GV triangular transform from the integers alone.
No analytic PT/GV generating-function formula is used in training.
Find maximum-consensus low-degree polynomials on Q6..20 data.
"""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations
import csv,json
from pathlib import Path
from run_experiments import lagrange_quad

def read(path):
 with open(path) as f: return list(csv.DictReader(f))

def interp_fit(xs,ys):
 # generic Lagrange interpolation, fraction only
 def f(z):
  out=Fraction(0)
  for i in range(len(xs)):
   v=Fraction(ys[i])
   for j in range(len(xs)):
    if i!=j: v*=Fraction(z-xs[j],xs[i]-xs[j])
   out+=v
  return out
 return f

def ransac_exact(points,deg):
 best=None
 for ii in combinations(range(len(points)),deg+1):
  xs=[points[i][0] for i in ii];ys=[points[i][1] for i in ii]
  if len(set(xs)) !=len(xs):continue
  p=interp_fit(xs,ys)
  count=sum(p(x)==v for x,v in points)
  # Prefer smaller fitted integer coefficients for tie-break.
  c=[p(i) for i in range(deg+1)]
  complexity=sum(abs(a.numerator).bit_length()+a.denominator.bit_length() for a in c)
  key=(count,-complexity)
  if best is None or key>best[0]:best=(key,p,ii,c)
 return best[1],best[0][0],best[2],best[3]

def fit(train_path):
 data=read(train_path);d={}
 for row in data:
  q=int(row['degree']);t=int(row['shell'])
  d[(q,t)]=(int(row['pt']),int(row['gv']),int(row['gc'])-int(row['delta'])+t)
 pts=[]
 for q in range(6,21):
  if (q,0) not in d or(q,1) not in d: continue
  p0,g0,h=d[(q,0)];p1,g1,_=d[(q,1)]
  assert p0==g0
  pts.append((h,Fraction(p1-g1,g0)))
 r1,num1,i1,c1=ransac_exact(pts,1)
 pts2=[]
 for q in range(6,21):
  if (q,2) not in d:continue
  p0,g0,h=d[(q,0)];p1,g1,_=d[(q,1)];p2,g2,_=d[(q,2)]
  coeff=Fraction(p2-g2-r1(h-1)*g1,g0)
  pts2.append((h,coeff))
 r2,num2,i2,c2=ransac_exact(pts2,2)
 assert num1>=12 and num2>=10,(num1,num2)
 return r1,r2,dict(r1_support=f'{num1}/{len(pts)}',r2_support=f'{num2}/{len(pts2)}',
                    r1_values=[str(r1(x)) for x in [3,4,5,6]],r2_values=[str(r2(x)) for x in [3,4,5,6]],
                    r1_outliers=[(q, str(Fraction(d[(q,1)][0]-d[(q,1)][1],d[(q,0)][1]) - r1(d[(q,0)][2]))) for q in range(6,21) if (q,1) in d and Fraction(d[(q,1)][0]-d[(q,1)][1],d[(q,0)][1])!=r1(d[(q,0)][2])],
                    r2_outliers=[(q,str(Fraction(d[(q,2)][0]-d[(q,2)][1]-r1(d[(q,0)][2]-1)*d[(q,1)][1],d[(q,0)][1])-r2(d[(q,0)][2]))) for q in range(6,21) if (q,2) in d and Fraction(d[(q,2)][0]-d[(q,2)][1]-r1(d[(q,0)][2]-1)*d[(q,1)][1],d[(q,0)][1]) !=r2(d[(q,0)][2])])

if __name__=='__main__':
 r1,r2,a=fit('data/calibration.csv')
 print('R1 h=3,4,5,6:',a['r1_values'],'support',a['r1_support']);print('R2:',a['r2_values'],'support',a['r2_support'])
 print('outliers',a['r1_outliers'],a['r2_outliers'])
 Path('results/empirical_pt_gv_transform.json').write_text(json.dumps(a,indent=2))
