"""Exact B-model high-X reconstruction with all unknown positive-degree GV withheld."""
import sys,json,time, argparse
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parent/'quintic_ai_predictor'
sys.path.insert(0,str(ROOT))
from hybrid.boundary import boundary_system
from engine.bcov_yy import Poly
D={r['genus']:r for r in json.loads((ROOT/'reference_amplitudes.json').read_text())['amplitudes']}
rows=[]
for g in [10,11,12]:
    t=time.monotonic()
    model=boundary_system(g, Poly.from_data(D[g]['anomaly']),mum_degrees=[0]);truth=[Fraction(a) for a in D[g]['ambiguity']]
    assert model.verify(truth)
    free=set(model.free_indices)
    fixed=[j for j in range(len(truth)) if all(v[j]==0 for v in model.kernel)]
    assert all(model.particular[j]==truth[j] for j in fixed)
    data={'genus':g,'exact_boundary_rank':model.rank,'total_coefficients':len(truth),'free_indices':list(model.free_indices), 'high_coefficients_fixed':list(range(max(free)+1,3*g-2)), 'all_fixed_match':True,'seconds':round(time.monotonic()-t,1)}
    rows.append(data)
    Path(__file__).resolve().parent/'boundary_verification.json'.write_text(json.dumps(rows,indent=2))
    print(data,flush=True)
