"""Exact finite-sequence closure checks. No extrapolation uses these post-hoc coefficients.
Tests small homogeneous polynomial-coefficient linear recurrences and Hankel rank.
Finite data can reject candidate orders but never prove finite generation.
"""
import json,sys,csv
from pathlib import Path
from fractions import Fraction as F
import sympy as sp
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from deep_envelope import W
from recursive_defects import read,phase_fit
from pt_from_dt_solver import pt_mC,simple_vacuum_pt
period,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))
DELTA={0:0,1:1,2:2,3:2,4:1}

def Bseq(cls,N=16):
 d={0:50,1:49,2:48}[cls];s=d%5;k=d//5;out=[]
 for t in range(N):
  m=pt_mC(d)+DELTA[s]+t
  out.append(F(simple_vacuum_pt(d,m),phase[s,k%period]*W(s,t,k)))
 return out

def matrix_rank(rows):
 return sp.Matrix(rows).rank() if rows else 0

out={}
for cls in range(3):
 b=Bseq(cls,16)
 ranks={n:matrix_rank([[sp.Rational(b[i+j].numerator,b[i+j].denominator) for j in range(n)] for i in range(n)]) for n in range(2,9)}
 # search exact nonzero null vector over polynomial-coefficient recurrences:
 # sum_{j=0..order}sum_{p=0..deg} u_{j,p} (n)^p B[n+j]=0.
 laws=[]
 for order,deg in [(1,0),(2,0),(3,0),(4,0),(5,0),(6,0),(7,0),(1,1),(2,1),(3,1),(1,2),(2,2),(3,2)]:
  rows=[]
  for n in range(len(b)-order):
   rows.append([sp.Rational(b[n+j].numerator,b[n+j].denominator)*(n**p) for j in range(order+1) for p in range(deg+1)])
  rnk=matrix_rank(rows)
  null=(order+1)*(deg+1)-rnk
  laws.append({'order':order,'degree':deg,'equations':len(rows),'null_dim':null,'rank':rnk})
 out[cls]={'hankel_ranks':ranks,'recurrence_tests':laws,'digits_t15':len(str(abs(b[15].numerator))),'B_0':str(b[0]),'B_15':str(b[15])}
 print('cls',cls,'hankel',ranks,'finite-order nulls',[(z['order'],z['degree'],z['null_dim']) for z in laws if z['null_dim']>0])
(ROOT/'results'/'closure_tests.json').write_text(json.dumps(out,indent=2))
