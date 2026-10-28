import json,csv
from pathlib import Path
from fractions import Fraction as Q
from experiment import dimensions, kernel_vector, section_midpoint, solve, matmul
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results'/'exact_lattice_audit.json').read_text())
assert r['n_slots']==22 and r['n_independent']==14
for obj in r['records']:
 g=obj['g']; K,m,D=dimensions(g)
 assert D==obj['D'] and obj['m']==m
 H=[[Q(x) for x in line] for line in obj['H']]
 n0=[Q(v) for v in obj['n_anchor']]
 c=[Q(v) for v in obj['physical_c']]
 chosen=solve(H,[-v for v in n0]); assert chosen==c
 assert matmul(H,c)==[-x for x in n0]
 assert all(H[i][j]==0 for i in range(D) for j in range(i+1,D))
 for j in range(D):
  step=solve(H,[Q(int(k==j)) for k in range(D)])
  assert matmul(H,step)==[int(k==j) for k in range(D)]
unit=json.loads((ROOT/'results'/'unit_charge_steps.json').read_text())
assert len(unit)==22
for u in unit:
 assert section_midpoint(u['g'],[Q(x) for x in u['delta_c']])==Q(u['delta_midpoint'])
transport=json.loads((ROOT/'results'/'integral_transport_degree26.json').read_text())
assert sum(x['n_integrality_entries'] for x in transport)==572
assert all(x['all_integer'] for x in transport)
assert all(all(isinstance(v,int) for v in vs) for x in transport for vs in x['unit_charge_response_selected'].values())
p=json.loads((ROOT/'results'/'ml_lattice_projection.json').read_text())
x=[r for r in p if r['split']=='train_to_8' and r['model']=='mixed_plus_genus__majority']
assert sum(row['n_exact'] for row in x)==6
assert sum(row['n_exact'] for row in x if row['independent_particular'])==2
print('All exact audits and held-out projection reproducibility tests pass.')
