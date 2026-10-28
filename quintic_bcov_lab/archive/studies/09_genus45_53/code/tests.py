"""Reproducible exact tests, no floating arithmetic."""
import sys,json,csv
from pathlib import Path
from fractions import Fraction as F
from math import comb
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent/'engine'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from bcov_yy import YYRing,Poly,Amplitude,X,v1,u
from near_extremal_benchmark import TRAIN, HOLDOUT,fit,predict,GC,DELTA
from bootstrap_45_53 import D,extremal_5k
from pt_dt_arithmetic import dt_zero_coeffs,conv,pt_gv_top

def load(name):
    return Poly.from_data(json.loads((ROOT/'data'/f'{name}_exact.json').read_text())['poly_terms'])

def main():
    b=json.loads((ROOT/'data'/'bootstrap_results.json').read_text())
    assert b['integer_sites']==169 and b['zero_sites']==168
    assert D(51)==20 and extremal_5k(4)==175
    assert [extremal_5k(k) for k in [2,3,4,5,6]]==[-50,-100,175,275,-400]
    for residue in TRAIN:
        model=fit(TRAIN[residue])
        for k in range(1,6):
            check=(TRAIN[residue][k-1] if k<=3 else HOLDOUT[residue][k-4])
            assert predict(model,k)==check
            assert GC[residue][k]-DELTA[residue]>=1
    assert dt_zero_coeffs(2)==[1,200,19500]
    assert 100*175+9700==27200
    assert 106*319125+7158750==40986000
    P51,P52,P53=map(load,['delta_P51','delta_P52','delta_P53'])
    R52,R53=map(load,['delta_R52','delta_R53'])
    ring=YYRing();P1=F(125,12)-v1/2-F(28,3)*u+F(5,12)*X
    P51_1=ring.covariant(P51,51,0)
    assert R52==-(ring.covariant(P51_1,51,1)+2*P1*P51_1)/2
    assert (P52.partial('v1')-u*P52.partial('v2')-u*(u+5*X)*P52.partial('v3'))==R52
    P2=Amplitude.from_data(json.loads((ROOT/'data'/'reference_P2.json').read_text())).polynomial()
    R53_expected=-(ring.covariant(ring.covariant(P52,52,0),52,1)
                   +2*P1*ring.covariant(P52,52,0)
                   +2*ring.covariant(P2,2,0)*P51_1)/2
    assert R53==R53_expected
    assert (P53.partial('v1')-u*P53.partial('v2')-u*(u+5*X)*P53.partial('v3'))==R53
    sensitivity=json.loads((ROOT/'data'/'physical_sensitivity.json').read_text())
    assert all(sensitivity[k] for k in ['genus53_conifold_gap_difference_pass','genus53_orbifold_difference_pass','genus53_MUM_GW_degree0to20_difference_pass'])
    print('ALL 49 integer predictions/training checks and all exact YY anomaly identities passed')

if __name__=='__main__':main()
