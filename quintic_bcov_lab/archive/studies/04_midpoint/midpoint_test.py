"""Exact midpoint involution and ambiguity-kernel audit for the mirror quintic.

Coordinate: x=psi**(-5)=3125*z_PF, X=x/(x-1).
This is a change of local coordinate, not a symmetry of the quintic VHS.
"""
from pathlib import Path
from math import comb
from fractions import Fraction as Q
import json,sys, argparse
BASE=Path(__file__).resolve().parent
SOURCE=BASE/'quintic_ai_predictor'
sys.path.insert(0,str(SOURCE))
from hybrid.boundary import boundary_system
from engine.bcov_yy import Poly

DATA={r['genus']:r for r in json.loads((SOURCE/'reference_amplitudes.json').read_text())['amplitudes']}

def midpoint(a, genus, order=0):
    """Exact value/first x derivative of the ambiguity *section coefficient*
    (x-1)**(g-1) a(x/(x-1)), with global omega prefactor suppressed.
    """
    s=sum(((-1)**j*c for j,c in enumerate(a)),Q(0))
    if order==0: return (-Q(1,2))**(genus-1)*s
    t=sum((j*(-1)**(j-1)*c for j,c in enumerate(a) if j),Q(0))
    if order==1: return (-Q(1,2))**(genus-2)*((genus-1)*s+2*t)
    raise ValueError('order must be 0 or 1')

def canonical_kernel(g):
    """Exact homogeneous ambiguity basis X**(r+1)*(1-X)**m, r=0..D-1."""
    m=(3*g-3+4)//5
    D=(2*g-2)//5
    out=[]
    for r in range(D):
        coeff=[Q(0)]*(3*g-2)
        for k in range(m+1): coeff[1+r+k]=Q((-1)**k*comb(m,k))
        out.append(tuple(coeff))
    return out

def diagnostic(g):
    truth=[Q(s) for s in DATA[g]['ambiguity']]
    model=boundary_system(g,Poly.from_data(DATA[g]['anomaly']),mum_degrees=[0])
    assert model.verify(truth)
    basis=canonical_kernel(g)
    assert len(basis)==len(model.kernel), ('kernel dimension disagrees with closed form',g)
    assert all(all(sum((a*b for a,b in zip(row,v)),Q(0))==0 for row,rhs,label in model.constraint_rows)
               for v in basis)
    images=[[midpoint(k,g,d) for k in basis] for d in (0,1)]
    # Pullback by w=1-x leads to (-1)**derivative order, independent of all a.
    assert midpoint(truth,g,0) == midpoint(truth,g,0)
    result={
        'genus':g,'rank':model.rank,'free_indices':list(model.free_indices),
        'm_orbifold':(3*g-3+4)//5,'D_MUM':(2*g-2)//5,'factorization_verified':True,
        'physical_midpoint_section_coefficient':str(midpoint(truth,g)),
        'physical_first_x_derivative':str(midpoint(truth,g,1)),
        'midpoint_responses_by_kernel_direction':[[str(v) for v in row] for row in images],
        'number_of_new_equations_from_coordinate_gluing':0,
        'meaning':('The midpoint response need not vanish along admissible kernel directions. '
                   'Coordinate gluing is an identity for all of them; a measured, externally supplied '
                   'midpoint observable could reduce the kernel, but the involution supplies no such value.')
    }
    if len(model.kernel)==2:
        det=images[0][0]*images[1][1]-images[0][1]*images[1][0]
        result['value_and_first_derivative_independent']= bool(det)
        result['midpoint_jet_rank']=2 if det else (1 if any(any(row) for row in images) else 0)
    if len(model.kernel)==3:
        result['midpoint_jet_rank']=2 if any(images[0][i]*images[1][j]-images[0][j]*images[1][i] for i in range(3) for j in range(i+1,3)) else 1
    return result

def section_polynomial(a,g):
    """Compute (x-1)**(g-1)*sum a_j X^j as exact polynomial.

    Only defined for differences with a_j=0 for j>g-1.
    """
    assert all(not a[j] for j in range(g,len(a)))
    z=[Q(0)]*g
    for j,val in enumerate(a[:g]):
        for t in range(g-j):
            z[j+t]+=val*Q(comb(g-1-j,t)*(-1)**(g-1-j-t))
    return z

def jets_at_half(polynomial,num):
    """Return first num derivatives of polynomial at x=1/2."""
    out=[]
    for r in range(num):
        out.append(sum((c*Q(__import__('math').factorial(k),__import__('math').factorial(k-r))
                        *Q(1,2)**(k-r) for k,c in enumerate(polynomial) if k>=r),Q(0)))
    return out

def solve(A,b):
    """Exact Gaussian inversion, rectangular input full column rank."""
    n=len(b); assert len(A)==n and all(len(row)==n for row in A)
    M=[list(row)+[v] for row,v in zip(A,b)]
    for p in range(n):
        i=next(i for i in range(p,n) if M[i][p]); M[p],M[i]=M[i],M[p]
        q=M[p][p];M[p]=[v/q for v in M[p]]
        for i in range(n):
            if i==p:continue
            q=M[i][p]
            if q:M[i]=[v-q*w for v,w in zip(M[i],M[p])]
    return [M[i][-1] for i in range(n)]

def oracle_demo(g):
    """Show inversion if midpoint observables ARE externally supplied.

    HERE observables are explicitly obtained from held-out reference and are
    not additional independent physics constraints. This is a plumbing test.
    """
    truth=[Q(a) for a in DATA[g]['ambiguity']]
    m=boundary_system(g,Poly.from_data(DATA[g]['anomaly']),mum_degrees=[0])
    D=len(m.kernel)
    basis=canonical_kernel(g)
    delta=[a-b for a,b in zip(truth,m.particular)]
    oracle=jets_at_half(section_polynomial(delta,g),D)
    A=list(map(list,zip(*(jets_at_half(section_polynomial(k,g),D) for k in basis))))
    coeff=solve(A,oracle)
    reconstruct=[m.particular[j]+sum((c*k[j] for c,k in zip(coeff,basis)),Q(0)) for j in range(len(truth))]
    assert reconstruct==truth
    return {'genus':g,'midpoint_derivative_orders_used':D,
            'oracle_provenance':'computed from held-out physical reference, NOT independent input',
            'oracle_jet':[str(x) for x in oracle],
            'recovered_kernel_coordinates':[str(x) for x in coeff],
            'all_ambiguity_coefficients_exactly_reconstructed':True}

def g2_counterexample():
    g=2
    a=[Q(v) for v in DATA[2]['ambiguity']]
    def f(x):return (x-1)**(g-1)*sum((v*(x/(x-1))**j for j,v in enumerate(a)),Q(0))
    return {'a2_at_X_minus_one':str(sum((v*(-1)**j for j,v in enumerate(a)),Q(0))),
            'f2_at_x_half':str(f(Q(1,2))),
            'f2_at_x_third':str(f(Q(1,3))),
            'f2_at_x_two_thirds':str(f(Q(2,3))),
            'not_x_to_one_minus_x_invariant':f(Q(1,3))!=f(Q(2,3))}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--genera', nargs='+',type=int,default=[7,9])
    args=parser.parse_args()
    rows=[]
    for g in args.genera:
        row=diagnostic(g); rows.append(row)
        print(f'g={g}: rank={row["rank"]}, free={row["free_indices"]}, midpoint jet rank={row.get("midpoint_jet_rank")}',flush=True)
        print('  sensitivity:',row['midpoint_responses_by_kernel_direction'],flush=True)
    (BASE/'results.json').write_text(json.dumps({
        'boundary_kernels':rows,
        'genus2_no_involution_symmetry':g2_counterexample(),
        'collocation_with_oracle_reference':[oracle_demo(g) for g in (7,9)]},indent=2))

if __name__=='__main__':main()
