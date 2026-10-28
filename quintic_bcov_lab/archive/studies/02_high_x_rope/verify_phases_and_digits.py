"""Exact GV finite-character identity and lossless rational numeral encoding."""
from fractions import Fraction
from math import comb
from pathlib import Path
import sys
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P/'quintic_ai_predictor'))
from hybrid.digits import encode_fraction,decode_fraction
from high_degree_phase_study import reference, gv_phase,gv_character

def multiply(a,b):
 out={}
 for j,v in a.items():
  for k,w in b.items():out[j+k]=out.get(j+k,0)+v*w
 return {j:v for j,v in out.items() if v}

def main():
 # Laurent identity (2-z-z^-1)^h = sum_r (-1)^(h+r) C(2h,r) z^(h-r)
 for h in range(1,15):
  lhs={0:1}
  for _ in range(h):lhs=multiply(lhs,{0:2,1:-1,-1:-1})
  rhs={h-r:(-1)**(h+r)*comb(2*h,r) for r in range(2*h+1)}
  assert lhs==rhs,(h,lhs,rhs)
 A=reference()
 n=0
 for g,a in A.items():
  for q in a:
   assert decode_fraction(encode_fraction(q))==q
   n+=1
 print('PASS exact Laurent identity for h=1..14 and exact digit round-trip for',n,'rational coefficients')
if __name__=='__main__':main()
