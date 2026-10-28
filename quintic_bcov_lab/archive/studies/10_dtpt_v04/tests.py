"""Verify exact rational DT modular reconstruction and quintic charge predictions."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
from modular_core import A,vec,h_series,verify_postfit
from pt_from_dt_solver import pt_from_rankzero_dt,omega_rankzero,qhat,applicable,wall_candidates,pt_mC
import empirical_benchmark as bench

assert A.rank()==9
assert verify_postfit()
assert omega_rankzero(1,F(83,120))==8625
assert omega_rankzero(2,F(-13,120))==-1218500
assert bench.learned=={1:8625,2:-1218500,0:5}
assert len(bench.ROWS)==25
assert len([r for r in bench.ROWS if r['degree']>10])==20
assert all(r['prediction_PT']==r['reference_PT_leading'] for r in bench.ROWS)

PT={-53:319125,-52:40986000,-51:132211590000,
    -50:-105493513413000,-49:19064998024136500}
for m,ref in PT.items():
    assert applicable(21,m)
    pred,baseline,walls=pt_from_rankzero_dt(21,m)
    assert pred==ref,(m,pred,ref)
assert wall_candidates(21,-53)==[]
assert pt_from_rankzero_dt(21,-50)[0]-pt_from_rankzero_dt(21,-50)[1]==-488750
assert pt_from_rankzero_dt(21,-49)[0]-pt_from_rankzero_dt(21,-49)[1]==-74922500

# PT->GV conversion with disconnected genus-51 degree-20 * genus-zero degree-1 pair.
gv={}
for m in sorted(PT):
    g=1-m
    known=sum(n*comb(2*h-2,h-g) for h,n in gv.items() if h>=g)
    gv[g]=PT[m]-known-(175*2875 if m==-49 else 0)
assert gv=={54:319125,53:7158750,52:129691149375,
            51:-118821918509250,50:30276261813046500}

# Frontier coefficients are not taken from published rows by the DT/PT predictors;
# these comparisons happen after the modular fit and wall-crossing calculation.
assert all(r['wall_inequality_pass'] and not r['nonvacuum_candidate_count']
           for r in bench.ROWS if r['degree']>=21)
assert all(r['theorem_single_wall_certified'] for r in bench.FORECAST)

from bcov_exact_response import scale_5,D,m,R,vals
assert scale_5==107 and D==21 and m==32
assert vals[:3]==[F(1),F(-3450),F(4790625)]
assert vals[3]==-3685910000
assert vals[4]==1751393550000
print('PASS: exact rank-zero DT polar reconstruction, 25 leading PT fits (20 withheld),')
print('      five degree-21 PT values, five PT->GV integers, wall optimality,')
print('      first disconnected product, BCOV genus-54 homogeneous response.')
