"""Regression tests: reference split, integral basis and exact boundary kernel."""
from fractions import Fraction as F
from pathlib import Path
import json
from intermediate import load_reference, reference_split, dimensions, kernel_vector
from quintic_ai_predictor.hybrid.boundary import boundary_system
from quintic_ai_predictor.engine.bcov_yy import Poly

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'data'/'reference_amplitudes.json'

def fast_tests():
    n=0; targets=0
    for g,a,anomaly in load_reference(SOURCE):
        K,m,D=dimensions(g)
        p=reference_split(g,a)
        assert len(p['c'])==D
        assert p['physical']==a
        assert p['D']==D
        assert p['pivot_indices']==list(range(m+1,g))
        assert all(a[j]==p['anchor'][j]+p['residual'][j] for j in range(K+1))
        assert all(not p['anchor'][j] for j in p['pivot_indices'])
        assert all(not p['residual'][j] for j in range(g,K+1))
        assert sum(p['bernstein_x'])==p['xQ'][1] if D==1 else True
        # Re-evaluate x-section for rational x by both bases.
        for num,den in [(1,2),(1,3),(3,7)]:
            x=F(num,den);X=x/(x-1)
            section=(x-1)**(g-1)*sum((aj*X**j for j,aj in enumerate(p['residual'])),F(0))
            poly=sum((aj*x**j for j,aj in enumerate(p['xQ'])),F(0))
            assert section==poly,(g,x)
        n+=1; targets+=D
    assert (n,targets)==(11,22)
    print('fast tests passed: 11 genera, 22 residual targets, exact section identities')

def slow_boundary_tests(max_genus=9):
    checked=[]
    for g,a,anomaly in load_reference(SOURCE):
        if g>max_genus:break
        p=reference_split(g,a)
        bm=boundary_system(g,Poly.from_data(anomaly),mum_degrees=[0])
        assert bm.verify(a)
        assert bm.free_indices==tuple(p['pivot_indices'])
        assert tuple(p['anchor'])==bm.particular
        for k in range(p['D']):
            v=kernel_vector(g,k)
            assert all(sum((coef*b for coef,b in zip(row,v)),F(0))==0 for row,_,_ in bm.constraint_rows)
        checked.append(g)
    print('exact independently solved boundary systems:',checked)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--boundary',action='store_true');ap.add_argument('--max-genus',type=int,default=9)
    a=ap.parse_args()
    fast_tests()
    if a.boundary:slow_boundary_tests(a.max_genus)
