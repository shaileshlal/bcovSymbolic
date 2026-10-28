"""A universal nearest-neighbor GV response for a shift of the LAST residual charge.

Derivation: for a unit top residual-site perturbation at d=D_g,
  delta F_g(q)= omega_0(z(q))**(2g-2) * z(q)**D,
with omega_0(z)=1+120z+..., q=z+770z**2+... .
Then [q^(D+1)]delta F_g=240*(g-1)-770*D.
When D=1, subtract the degree-1 twice-covered GV contribution 2**(2g-3).
Independent exact transport values are in previous published/archived tests.
"""
from __future__ import annotations
import json
from fractions import Fraction as Q
from pathlib import Path
from build_data import residual_dimension
ROOT=Path(__file__).resolve().parent


def period_mirror_first_coeffs():
    a=[Q(j,5) for j in range(1,5)]
    first_period=Q(1)
    for v in a:first_period*=v
    log_period=first_period*(sum((1/x for x in a),Q(0))-4)
    return {'period_omega0_z1':first_period*3125,
            'mirror_q_z2':log_period*3125}


def terminal_neighbor_transport(g,d=None):
    if d is None:d=residual_dimension(g)
    assert d==residual_dimension(g),'Formula requires terminal residual site'
    return 240*(g-1)-770*d-(2**(2*g-3) if d==1 else 0)


def main():
    source=ROOT/'source'/'integral_transport_degree26.json'
    data=json.loads(source.read_text())
    assert period_mirror_first_coeffs()=={'period_omega0_z1':120,'mirror_q_z2':770}
    records=[]
    for e in data:
        g,d=e['g'],e['D']
        true=e['unit_charge_response_selected'][str(d)][d]
        pred=terminal_neighbor_transport(g,d)
        assert pred==true,(g,d,pred,true)
        records.append({'g':g,'d':d,'observed_exact':true,'predicted_formula':pred,
                        'exception_double_cover':d==1})
    high=[]
    for g in [51,54,55,56,76]:
        d=residual_dimension(g)
        high.append({'g':g,'d_terminal':d,'coefficient':terminal_neighbor_transport(g),
                     'predicted_change_per_unit_charge':'exact structural response, not physical charge prediction'})
    out={'period_mirror_first_coeffs':{k:int(v) for k,v in period_mirror_first_coeffs().items()},
         'observations':records,'checked':len(records),
         'formula':'240(g-1)-770D-1_(D=1)*2^(2g-3)',
         'higher_genus_structural_responses':high,
         'first_actual_nonzero_case':{'g':51,'d':20,'physical_charge':175,
                'exact_change_next_degree':175*terminal_neighbor_transport(51)},
         'frontier_unknown_template':{'g':54,'d':21,
                'delta_n_54_22_per_delta_n_54_21':terminal_neighbor_transport(54)}}
    (ROOT/'results'/'neighbor_regulator.json').write_text(json.dumps(out,indent=2)+'\n')
    print('neighbor responses verified',len(records),'of',len(records))
    print('g51 actual charge shift at degree21:',175*terminal_neighbor_transport(51))
    print('g54 response per unit:',terminal_neighbor_transport(54))

if __name__=='__main__':main()
