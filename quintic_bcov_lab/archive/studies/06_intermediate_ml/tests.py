"""Reproducibility and safety checks for the intermediate-only benchmark."""
from pathlib import Path
from fractions import Fraction
import json
import run_experiment as e
r=e.load_rows()
assert len(r)==22
assert sorted({v['genus'] for v in r})==list(range(4,13))
for g in range(4,13):
    a=[z for z in r if z['genus']==g]
    assert len(a)==(2*g-2)//5
    assert len({x['sign'] for x in a})==1
    for x in a:
        assert Fraction(x['exact_b'])==(-1)**(g-2-x['r'])*Fraction(x['exact_c'])
        assert x['j']==(3*g-3+4)//5+1+x['r']
# No reference-derived true physical coefficients allowed in model features.
import inspect
src=inspect.getsource(e.feat)
for bad in ('a_pivot_phys','physical','anchor','true_exact_c','exact_c'):
    assert bad not in src
for tag in ('genus','mixed','mixed_plus_genus','mixed_residual','v01_trees'):
    p,_=e.run(r,8,tag,'combined_v01' if tag=='v01_trees' else 'majority')
    assert len(p)==14 and all(z['genus']>=9 for z in p)
    assert all(0<=z['error_dex']<100 for z in p)
res=json.loads((Path(e.__file__).parent/'results/benchmark.json').read_text())
assert res['dataset']['sign_uniform_within_genus']
assert res['shape_tests']['8']['n_train_shape']==3
print('All exact label/basis, feature-leakage, and chronological holdout checks passed.')
