"""Fresh small-genus integral-coset reconstruction; no Airy input.

The genus-three denominator cap is a retrospective empirical input, NOT
deduced from GV integrality. The minimum-absolute-GV rule is a heuristic.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import json,sys
sys.path.insert(0,str(Path(__file__).parent/'engine'))
import bcov_yy as alg
from physical_boundaries import seed_gv,record_gv
from leading_only import bernoulli_abs,constant_coefficients,leading_airy

def select(intercept,response,D):
    # m/D in the coefficient lattice, N = intercept + response*m/D integral.
    A=response.numerator*intercept.denominator
    B=-intercept.numerator*D*response.denominator
    modulus=intercept.denominator*D*response.denominator
    common=gcd(A,modulus)
    assert B%common==0
    period=modulus//common
    m0=(B//common*pow(A//common,-1,period))%period if period>1 else 0
    n0=intercept+response*Q(m0,D);step=response*Q(period,D)
    assert n0.denominator==step.denominator==1
    z=-n0/step;floor=z.numerator//z.denominator
    candidates=[(abs(n0+step*j),j) for j in [floor,floor+1]]
    _,jump=min(candidates)
    m=m0+period*jump
    return Q(m,D),n0+step*jump,step

def main():
    state=alg.BCOV();gv=seed_gv(9);rows=[];bb=bernoulli_abs(3);c=constant_coefficients(bb)
    for g,D in [(2,144),(3,36288)]:
        anomaly=state.integrate(g).anomaly;coeff=[c[g]]
        for d in range(1,3*g-2):
            partial=anomaly+sum((v*alg.X**j for j,v in enumerate(coeff)),alg.Poly())
            trial=dict(gv);record_gv(g,partial,d,trial);intercept=trial[g,d]
            response=(-1)**(g-1+d)*Q(5)**(5*d-4*g+4)
            chosen,ngd,step=select(intercept,response,D)
            # Check the triangular response directly, without relying on its formula.
            trial2=dict(gv);record_gv(g,partial+alg.X**d,d,trial2)
            assert trial2[g,d]-intercept==response
            coeff.append(chosen)
            rows.append({'g':g,'degree':d,'cap':D,'intercept':str(intercept),
                         'response':str(response),'chosen_coefficient':str(chosen),
                         'chosen_GV':str(ngd),'GV_lattice_step':str(step)})
        state.set_ambiguity(g,coeff);assert not state.residual(g)
        record_gv(g,state.amplitudes[g].polynomial(),9,gv)
    expected={(2,4):534750,(2,5):75478987900,(3,4):8625,(3,5):-15663750,(3,6):3156446162875}
    for key,v in expected.items():assert gv[key]==v,(key,gv[key],v)
    reference=json.loads((Path(__file__).parent/'reference_amplitudes.json').read_text())
    for a in reference['amplitudes']:
        if a['genus'] in [2,3]:
            assert list(map(Q,a['ambiguity']))==list(state.amplitudes[a['genus']].ambiguity)
    A,B=128340000,-217350000
    assert A%240==0 and (10*B+21*A)%60480==0
    result={'inputs':{'seed':'genus-one vertex and PF/MUM period data','denominator_caps':{'2':144,'3':36288},
                     'selection':'minimum absolute GV representative; empirical, no uniqueness theorem'},
            'smith_invariants':[24,60480],'rows':rows,
            'coefficients':{str(g):list(map(str,state.amplitudes[g].ambiguity)) for g in [2,3]},
            'GV':{str(g):{str(d):str(gv[g,d]) for d in range(1,10)} for g in [2,3]},
            'normalization_invariance':'Multiplying the entire candidate lattice by a known nonzero scalar does not select a different minimum-absolute-GV member.'}
    (Path(__file__).parent/'results'/'lattice_audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
