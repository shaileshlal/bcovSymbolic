"""One-command, fail-fast tests for data provenance, integer predictions and boundary transport."""
from __future__ import annotations
from pathlib import Path
import csv, gzip, json, math
from fractions import Fraction as Q
from build_data import top_extremal, residual_dimension, margin10
from frontier_solver import classify_charge,first_nonzero_section_taylor
from neighbor_regulator import terminal_neighbor_transport,period_mirror_first_coeffs
from exact_solver import lower_triangular_solve,matvec,section_bernstein_coeffs
ROOT=Path(__file__).resolve().parent


def run():
    manifest=json.loads((ROOT/'data'/'dataset_manifest.json').read_text())
    assert manifest['total']==199400
    assert (manifest['certified_zero'],manifest['extremal_nonzero'],manifest['unresolved'])==(63544,16,135840)
    # Leakage audit: unresolved labels have empty displacement strings, even when in genus 54.
    with gzip.open(ROOT/'data'/'all_lattice_sites_g2_g1000.csv.gz','rt') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==199400
    assert sum(r['label_status']=='unresolved' for r in rows)==135840
    assert all(r['displacement']=='' for r in rows if r['label_status']=='unresolved')
    assert any(r['g']=='54' and r['d']=='21' and r['displacement']=='' for r in rows)
    assert classify_charge(51,20)['delta_n']==175
    assert classify_charge(52,20)['delta_n']==0
    assert classify_charge(53,20)['delta_n']==0
    assert classify_charge(54,21)['delta_n'] is None
    assert first_nonzero_section_taylor(51,20,175)['coefficient_exact']==str(175*5**100)
    assert first_nonzero_section_taylor(54,21,1)['coefficient_exact']==str(5**107)
    assert top_extremal(4)==(51,20,175)
    assert top_extremal(5)==(76,25,275)
    assert top_extremal(6)==(106,30,-400)
    # The exact symbolic rule was fitted only k=4..11, and tests all eight k12..19.
    models=json.loads((ROOT/'results'/'model_benchmarks.json').read_text())
    assert models['learned']['sign_period']==4
    assert models['learned']['minimum_polynomial_degree']==2
    assert models['metrics']['mod4_plus_quadratic_abs']['exact_integers']==8
    assert models['metrics']['rotary_phase4_poly2']['exact_integers']==8
    for name in ('ridge_quadratic_raw','random_forest_raw','extra_trees_raw','kNN_k3_raw'):
        assert models['metrics'][name]['exact_integers']==0
    # Exact general nearest-neighbor regulation with the D=1 multicover exception.
    assert period_mirror_first_coeffs()=={'period_omega0_z1':120,'mirror_q_z2':770}
    historical=json.loads((ROOT/'source'/'integral_transport_degree26.json').read_text())
    for e in historical:
        g,d=e['g'],e['D']
        assert terminal_neighbor_transport(g,d)==e['unit_charge_response_selected'][str(d)][d]
    assert terminal_neighbor_transport(51)==-3400
    assert terminal_neighbor_transport(54)==-3450
    assert 175*terminal_neighbor_transport(51)==-595000
    # All 22 old exact residual coordinates give n_min=0 and 0 physical displacement.
    audit=json.loads((ROOT/'source'/'exact_lattice_audit.json').read_text())
    assert sum(e['D'] for e in audit['records'])==22
    for e in audit['records']:
        g,D=e['g'],e['D']
        H=[[Q(t) for t in row] for row in e['H']]
        prior=[Q(q) for q in e['n_anchor']]
        c=lower_triangular_solve(H,[-p for p in prior])
        assert c==[Q(q) for q in e['physical_c']]
        assert all(Q(n)==0 for n in e['physical_GV'])
        for d in range(1,D+1):
            v=[0]*D;v[d-1]=1
            step=lower_triangular_solve(H,v)
            assert matvec(H,step)==v
            # Leading section coefficient of the perturbation is universal.
            b=section_bernstein_coeffs(g,step)
            assert b[d-1]==Q(5)**(4*g-4-5*d)
            assert all(b[i]==0 for i in range(d-1))
        last=[0]*(D-1)+[1]
        b=section_bernstein_coeffs(g,lower_triangular_solve(H,last))
        assert b==[Q(0)]*(D-1)+[Q(5)**(4*g-4-5*D)]
    print('PASS: data partition, no unknown/zero leakage, 8/8 extremal holdouts, 22 exact rational calibration, 9 neighborhood laws.')

if __name__=='__main__':run()
