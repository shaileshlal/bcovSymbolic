"""Exact triangular GV-lattice solver on the bundled residual calibration cases.

This is *not* an all-genus BCOV engine: H and affine anchor for g=4..12 were
computed by the previous exact BCOV release, with g11..12 retrospective.
"""
from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
import json, csv
ROOT=Path(__file__).resolve().parent


def lower_triangular_solve(H, rhs):
    assert len(H)==len(rhs)
    result=[]
    for i in range(len(H)):
        assert all(H[i][j]==0 for j in range(i+1,len(H)))
        assert H[i][i]
        result.append((rhs[i]-sum((H[i][j]*result[j] for j in range(i)),Q(0)))/H[i][i])
    return result


def matvec(H,v):
    return [sum((a*b for a,b in zip(row,v)),Q(0)) for row in H]


def physical_c_from_charges(entry,charges):
    H=[[Q(v) for v in row] for row in entry['H']]
    n0=[Q(v) for v in entry['n_anchor']]
    return lower_triangular_solve(H,[Q(n)-a for n,a in zip(charges,n0)])


def section_bernstein_coeffs(g,c):
    return [Q((-1)**(g-2-r))*a for r,a in enumerate(c)]


def main():
    audit=json.loads((ROOT/'source'/'exact_lattice_audit.json').read_text())
    records=[]
    for e in audit['records']:
        g=e['g']; D=e['D']
        H=[[Q(v) for v in row] for row in e['H']]
        assert all(H[i][j]==Q(0) for i in range(D) for j in range(i+1,D))
        for d in range(1,D+1):
            assert H[d-1][d-1] == Q((-1)**(g-1+d))*Q(5)**(5*d-4*g+4)
        c0=physical_c_from_charges(e,[0]*D)
        assert c0==[Q(v) for v in e['physical_c']]
        for i in range(D):
            charge=[0]*D;charge[i]=1
            cs=lower_triangular_solve(H,charge)
            assert matvec(H,cs)==[Q(v) for v in charge]
            new_c=[a+b for a,b in zip(c0,cs)]
            # Correct affine check: n_anchor + H*(c0+cs) = unit charge.
            assert [u+v for u,v in zip([Q(q) for q in e['n_anchor']],matvec(H,new_c))] == charge
        records.append({'g':g,'D':D,'physical_charge':[0]*D,
                        'exact_minimum_match':True,
                        'unit_shift_directions_verified':D,
                        'boundary_particular_independent':bool(e['independent_particular_matches_reference_split']),
                        'first_bernstein_coefficient':str(section_bernstein_coeffs(g,c0)[0])})
    # Prior held-out ML errors are retrospective charge mismatches, not physical nonzero labels.
    predictions=[{'g':9,'model_charges':[0,-1,3400],'true_charges':[0,0,0]},
                 {'g':10,'model_charges':[0,-1,3231],'true_charges':[0,0,0]}]
    shifts=[]
    for p in predictions:
        e=next(z for z in audit['records'] if z['g']==p['g'])
        dc=lower_triangular_solve([[Q(v) for v in row] for row in e['H']],
                                  [Q(a-b) for a,b in zip(p['true_charges'],p['model_charges'])])
        shifts.append({'g':p['g'],'correction_in_GV':[a-b for a,b in zip(p['true_charges'],p['model_charges'])],
                       'exact_shift_in_Xbasis':[str(x) for x in dc]})
    out={'calibration':records,'model_error_corrections':shifts,
         'validation':'exact rational triangular inversion, historical g4..12 archive'}
    (ROOT/'results'/'exact_solver_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Exact residual checks:',sum(e['D'] for e in records),'of',sum(e['D'] for e in records))
    for item in shifts:print(item)

if __name__=='__main__':main()
