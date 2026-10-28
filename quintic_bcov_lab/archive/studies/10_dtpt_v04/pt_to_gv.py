"""Extract near-Castelnuovo GV integers from the predicted degree-21 PT jet.

Includes the first disconnected contribution (degree 20 + degree 1)
which becomes visible at m=-49, avoiding an incorrect 503125 offset.
"""
from math import comb
import json
from pathlib import Path
from pt_from_dt_solver import pt_from_rankzero_dt

TOP=54
PT={m:pt_from_rankzero_dt(21,m)[0] for m in range(-53,-48)}
DISCONNECTED={-49:175*2875}
GV={}
for m in sorted(PT):
    g=1-m
    known=sum(v*comb(2*h-2,h-g) for h,v in GV.items() if h>=g)
    disconnected=DISCONNECTED.get(m,0)
    GV[g]=PT[m]-known-disconnected

REF={54:319125,53:7158750,52:129691149375,51:-118821918509250,50:30276261813046500}
assert GV==REF,(GV,REF)
assert DISCONNECTED[-49]==503125

result={'DT_PT_predictions':{str(m):str(v) for m,v in PT.items()},
        'GV_predictions_after_plethystic_subtraction':{str(g):str(v) for g,v in GV.items()},
        'degree20_degree1_first_disconnected_at_m_minus49':'175*2875=503125',
        'status':'Exact GV back-transformation conditional on the rank-zero modular ansatz and wall-crossing hypotheses.'}
Path(__file__).parent.joinpath('results').mkdir(exist_ok=True)
Path(__file__).parent.joinpath('results/pt_to_gv_g54.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
