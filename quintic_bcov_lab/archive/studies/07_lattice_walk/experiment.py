"""Quintic GV-lattice walking on the exact residual x(1-x) ambiguity.

Reads the previously released corrected genus-2..12 amplitudes. The affine
boundary particular and response matrix are calculated independently of target
genus pure-X labels. All arithmetic in Fraction. Reference used only for
validation and retrospective denominator-cap audit.
"""
from __future__ import annotations
import sys, json, csv, math, time
from pathlib import Path
from fractions import Fraction as Q
from functools import reduce
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'Quintic_Intermediate_ML'/'source'
sys.path[:0]=[str(SOURCE),str(SOURCE/'quintic_ai_predictor'/'engine'),str(SOURCE/'quintic_ai_predictor'/'hybrid'),str(ROOT/'quintic_ai_predictor'/'hybrid')]
import bcov_yy as alg
from exact_series import evaluate_poly
from period_frames import frame
from physical_boundaries import seed_gv, record_gv, q_coefficients, required_gv
from boundary_fast import boundary_system
from intermediate import dimensions, kernel_vector, reference_split
from lattice import integral_coset

def s(a):return str(a)
def lcm(a,b):return math.lcm(a,b)
def poly(coeffs):return sum((v*alg.X**j for j,v in enumerate(coeffs) if v),alg.Poly())
def basis_poly(g,r):return poly(kernel_vector(g,r))

def section_midpoint(g, c):
    K,m,D=dimensions(g)
    return sum((Q((-1)**(g-2-r),2**D)*v for r,v in enumerate(c)),Q(0))

def M_matrix(g):
    """Exact d-by-r response GV(g,d) / residual c_r. Lower triangular."""
    K,m,D=dimensions(g)
    vals,k,_=frame('mum',D+9)
    H=[[Q(0) for _ in range(D)] for d in range(D)]
    for r in range(D):
        gw=q_coefficients(evaluate_poly(basis_poly(g,r),vals)*k**(g-1),D)
        dn=[Q(0)]*(D+1)
        for d in range(1,D+1):
            dn[d]=gw[d]-sum((Q(a)**(2*g-3)*dn[d//a] for a in range(2,d+1) if d%a==0),Q(0))
            H[d-1][r]=dn[d]
    assert all(H[i][j]==0 for i in range(D) for j in range(i+1,D))
    for d in range(1,D+1):
        response=Q((-1)**(g-1+d))*Q(5)**(5*d-4*g+4)
        assert H[d-1][d-1]==response,(g,d,H[d-1][d-1],response)
    return H

def solve(H, v):
    n=len(H);c=[]
    for i in range(n):c.append((v[i]-sum((H[i][j]*c[j] for j in range(i)),Q(0)))/H[i][i])
    return c

def matmul(H,c):return [sum((H[i][j]*c[j] for j in range(len(c))),Q(0)) for i in range(len(c))]

def do():
    t0=time.monotonic()
    data=json.loads((SOURCE/'reference_amplitudes.json').read_text())['amplitudes']
    gv=seed_gv(4)
    records=[];step_rows=[];unit_rows=[]
    for item in data:
        g=item['genus'];K,m,D=dimensions(g)
        a=[Q(v) for v in item['ambiguity']]
        anomaly=alg.Poly.from_data(item['anomaly'])
        full=anomaly+poly(a)
        if D:
            t=time.monotonic()
            split=reference_split(g,a)
            independent=(g<=10)
            if independent:
                boundary=boundary_system(g,anomaly,mum_degrees=[0])
                assert list(boundary.free_indices)==list(range(m+1,g))
                assert boundary.verify(a)
                assert list(boundary.particular)==split['anchor'],('independent anchor mismatch',g)
                assert all(boundary.verify([z+y for z,y in zip(boundary.particular,kernel_vector(g,r))]) for r in range(D))
                anchor=boundary.particular
                rank=boundary.rank
            else:
                anchor=split['anchor'] # RETROSPECTIVE ONLY: dependent on target genus reference
                rank=3*g-2-D
            H=M_matrix(g)
            pre=dict(gv)
            gv0=dict(pre)
            record_gv(g,anomaly+poly(anchor),D,gv0)
            baseline=[gv0[g,d] for d in range(1,D+1)]
            c_phys=split['c']
            assert [baseline[i]+val for i,val in enumerate(matmul(H,c_phys))]==[Q(0)]*D,('GV relative check',g)
            # No same-genus target labels used to compute c_min: only independent
            # BCOV affine particular, lower-genus GV, and minimum |n|=0.
            c_min=solve(H,[-v for v in baseline])
            assert c_min==c_phys,('minimum reconstruction error',g)
            cap=reduce(lcm,(z.denominator for z in c_phys),1) # ORACLE audit, never predictive
            chosen=[]
            for r in range(D):
                b=baseline[r]+sum((H[r][j]*chosen[j] for j in range(r)),Q(0))
                coset=integral_coset(b,H[r][r],cap)
                guess,ngv=coset.closest_index(0)
                assert guess==c_phys[r] and ngv==0,(g,r,'cap did not contain physical minimum')
                chosen.append(guess)
                step_rows.append({'g':g,'d':r+1,'r':r,'b_x_sign':(-1)**(g-2-r),
                    'c_exact':s(c_phys[r]),'response':s(H[r][r]),'n_anchor_after_previous':s(b),
                    'selected_GV':ngv,'minheight_c':s(guess),'retrospective_denominator_cap':cap,
                    'coset_n_step':s(coset.n_step),'coset_c_step':s(Q(coset.m_step,cap)),
                    'offset_in_coset_steps':0})
            unit_effect=[]
            for d in range(1,D+1):
                e=[Q(int(i+1==d)) for i in range(D)]
                step=solve(H,e) # exact residual coordinate translation for GV unit charge
                assert matmul(H,step)==e
                trial=[p+q for p,q in zip(c_phys,step)]
                if independent:
                    assert boundary.verify([p+z for p,z in zip(anchor,
                        [sum((step[r]*kernel_vector(g,r)[j] for r in range(D)),Q(0))+
                         sum((c_phys[r]*kernel_vector(g,r)[j] for r in range(D)),Q(0)) for j in range(K+1)])])
                unit_rows.append({'g':g,'d':d,'delta_c':[s(v) for v in step],
                    'delta_midpoint':s(section_midpoint(g,step)),
                    'section_basis':'sum (-1)^(g-2-r) dc_r x^(r+1)(1-x)^(D-r-1)'})
            records.append({'g':g,'D':D,'K':K,'m':m,'boundary_rank':rank,
                'independent_particular_matches_reference_split':independent,
                'lattice_triangular':True,'diagonal':[s(H[r][r]) for r in range(D)],
                'H':[[s(z) for z in row] for row in H],
                'n_anchor':[s(z) for z in baseline],
                'minimum_c':[s(v) for v in c_min],
                'physical_c':[s(v) for v in c_phys],
                'minimum_GV':[s(z) for z in [baseline[i]+z for i,z in enumerate(matmul(H,c_min))]],
                'physical_GV':[s(gv.get((g,d),'not_yet_evaluated')) for d in range(1,D+1)],
                'reference_denominator_cap':cap,
                'elapsed_boundary_and_matrix_s':round(time.monotonic()-t,2)})
            print(f'g={g:2} D={D} {"INDEPENDENT" if independent else "REFERENCE-DERIVED"} minheight={D}/{D}, {time.monotonic()-t:.2f}s',flush=True)
        record_gv(g,full,4,gv)
        if D:assert all(gv[g,d]==0 for d in range(1,D+1))
    # Correct hindsight set for n[g,d] to valid exact physical 0
    for rec in records:rec['physical_GV']=[s(gv[rec['g'],d]) for d in range(1,rec['D']+1)]
    output={'protocol':{'prediction_data':'known genus-g anomaly, independent conifold/orbifold/MUM-constant-map affine particular, lower-genus physically fixed amplitudes, q/mirror map',
       'hidden_data':'all target-genus pure-X coefficients; target-genus positive degree GV values until scoring',
       'selection':'minimize absolute integer GV degree by degree (0); test with retrospective physical coefficient denominator cap separately',
       'caution':'all observed predicted GV=0 geometrically; these are control cases, not test of nonzero lattice correction. g11-12 particular derived from held-out physical labels and are NOT independent predictions'},
       'n_genera':len(records),'n_slots':len(step_rows),'n_exact':len(step_rows),'n_independent':sum(r['D'] for r in records if r['independent_particular_matches_reference_split']),
       'records':records,'unit_GV_steps':unit_rows}
    (ROOT/'results').mkdir(exist_ok=True)
    (ROOT/'results'/'exact_lattice_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    with (ROOT/'results'/'lattice_coordinates.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=step_rows[0].keys());w.writeheader();w.writerows(step_rows)
    with (ROOT/'results'/'unit_charge_steps.json').open('w') as f:json.dump(unit_rows,f,indent=2)
    print('TOTAL:',len(step_rows),'of',len(step_rows),'min-height recoveries. Time',round(time.monotonic()-t0,1),'s')

if __name__=='__main__':do()
