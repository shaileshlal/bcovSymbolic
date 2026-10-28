"""Check exact q-series support of X ambiguity via the physical period frame."""
import sys,json
from pathlib import Path
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parent/'quintic_ai_predictor'
sys.path.insert(0,str(ROOT/'engine'))
from physical_boundaries import q_coefficients,Jet,frame
D=[]
for g in (6,9,10,11,12):
    values,k,_=frame('mum',14)
    xp=Jet({0:Q(1)},14)
    cols=[]
    for j in range(3*g-2):
        cols.append(q_coefficients(xp*(k**(g-1)),4))
        xp=xp*values[4]
    diag=[]
    for d in range(5):
        for j in range(d+1,3*g-2):assert cols[j][d]==0,(g,j,d)
        want=Q((-1)**(g-1+d)*5**(5*d),625**(g-1))
        assert cols[d][d]==want,(g,d,cols[d][d],want)
        diag.append(str(want))
    D.append({'genus':g,'MUM_q_degree':4,'response_jacobian_rank':5,
              'zero_columns_for_j_greater_than':4,
              'number_high_tail_columns_with_exact_zero_response':3*g-2-5,
              'diagonal':[str(z) for z in diag]})
P=Path(__file__).resolve().parent
(P/'exact_response_audit.json').write_text(json.dumps(D,indent=2)+'\n')
print(json.dumps(D,indent=2))
