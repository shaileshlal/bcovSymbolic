"""Blinded-by-time symbolic/ML benchmark on the proven extremal n[g,5k].

Train extremal k=4..11; holdout k=12..19. Proven formula is only used to
create ground-truth labels (in build_data), not to choose model's rule.
This experiment does *not* predict unresolved interior GV charges.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction
import collections, csv, json, math
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Ridge
from sympy import Symbol, Rational, interpolate, simplify

ROOT=Path(__file__).resolve().parent
rows=list(csv.DictReader((ROOT/'data'/'extremal_displacements.csv').open()))
tr=[r for r in rows if r['split']=='train']
te=[r for r in rows if r['split']=='test']
assert len(tr)==len(te)==8
ktr=np.array([int(r['k']) for r in tr],dtype=np.float64)
kts=np.array([int(r['k']) for r in te],dtype=np.float64)
ytr=np.array([int(r['n']) for r in tr],dtype=np.float64)
yts=np.array([int(r['n']) for r in te],dtype=np.float64)


def find_shortest_sign_period(train):
    """Smallest residue period supported by *at least two instances per residue*."""
    for p in range(1,9):
        d=collections.defaultdict(list)
        for a in train:d[int(a['k'])%p].append(1 if int(a['n'])>0 else -1)
        if len(d)!=p or any(len(v)<2 or len(set(v))>1 for v in d.values()):continue
        return p,{str(k):v[0] for k,v in sorted(d.items())}
    return None,{}


def smallest_exact_abs_polynomial(train, max_degree=5):
    k=Symbol('k')
    for degree in range(max_degree+1):
        points=[(int(r['k']),abs(int(r['n']))) for r in train[:degree+1]]
        candidate=simplify(interpolate(points,k))
        if all(candidate.subs(k,int(r['k']))==abs(int(r['n'])) for r in train):
            return candidate,degree
    return None,None

p, signs=find_shortest_sign_period(tr)
P,degree=smallest_exact_abs_polynomial(tr)
assert p==4 and degree==2 and str(P)==str(25*(Symbol('k')**2-Symbol('k')+2)/2)
signpred=np.array([signs[str(int(k)%p)] for k in kts],int)
pred_symbolic=np.array([int(P.subs(Symbol('k'),int(k)))*sg for k,sg in zip(kts,signpred)],dtype=float)

# Compare to plausible small-data black-box predictors with no oracle target-genus input.
models={
    'zero_origin':np.zeros_like(yts),
    'last_observed_integer':np.full_like(yts,ytr[-1]),
    'ridge_quadratic_raw':None,
    'kNN_k3_raw':None,
    'random_forest_raw':None,
    'extra_trees_raw':None,
    'mod4_plus_quadratic_abs':pred_symbolic,
}
Xtr=np.stack([ktr,ktr**2],axis=1);Xts=np.stack([kts,kts**2],axis=1)
models['ridge_quadratic_raw']=Ridge(alpha=1e-6).fit(Xtr,ytr).predict(Xts)
models['kNN_k3_raw']=KNeighborsRegressor(n_neighbors=3).fit(ktr.reshape(-1,1),ytr).predict(kts.reshape(-1,1))
models['random_forest_raw']=RandomForestRegressor(n_estimators=256,random_state=2026,min_samples_leaf=1).fit(ktr.reshape(-1,1),ytr).predict(kts.reshape(-1,1))
models['extra_trees_raw']=ExtraTreesRegressor(n_estimators=256,random_state=2026,min_samples_leaf=1).fit(ktr.reshape(-1,1),ytr).predict(kts.reshape(-1,1))


def model_metric(pred):
    return dict(mae=float(np.mean(np.abs(pred-yts))),
        median_ae=float(np.median(np.abs(pred-yts))),
        exact_integers=int(np.sum(np.rint(pred).astype(np.int64)==yts)),
        sign_accuracy=float(np.mean(np.sign(pred)==np.sign(yts))),
        max_abs_error=float(np.max(np.abs(pred-yts))))

metrics={name:model_metric(pred) for name,pred in models.items()}

# Continuous RoPE-style phase representation: 2D rotation at theta=pi*k/2.
# Compare its learned signed polynomial envelope to the symbolic sign+poly route.
def phase_design(ks,period,degree):
    phi=2*np.pi*np.asarray(ks)/period
    return np.stack([f*ks**j for j in range(degree+1)
                     for f in (np.cos(phi),np.sin(phi))],axis=1)
phase=phase_design(ktr,4,2)
phase_coef, *_=np.linalg.lstsq(phase,ytr,rcond=None)
phase_pred=phase_design(kts,4,2)@phase_coef
models['rotary_phase4_poly2']=phase_pred
metrics['rotary_phase4_poly2']=model_metric(phase_pred)

# Pre-holdout phase-frequency selection: fit on k=4..9, choose with k=10..11.
# Keep only models with 6 or fewer effective dimensions, rank sufficient, and no hidden holdout labels.
phase_candidates=[]
for period in range(2,9):
    for dg in range(0,3):
        X=phase_design(ktr[:6],period,dg)
        if np.linalg.matrix_rank(X,tol=1e-7)<X.shape[1]:continue
        z,_r,_rank,_s=np.linalg.lstsq(X,ytr[:6],rcond=None)
        pp=phase_design(ktr[6:],period,dg)@z
        err=float(np.mean(np.abs(pp-ytr[6:])))
        phase_candidates.append({'period':period,'degree':dg,'validation_mae':err,
                                  'parameter_count':X.shape[1]})
phase_candidates.sort(key=lambda a:(round(a['validation_mae'],9),a['parameter_count'],a['period']))

with (ROOT/'results'/'heldout_predictions.csv').open('w',newline='') as f:
    cols=['k','g','d','physical_integer']+list(models)
    w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
    for i,row in enumerate(te):
        w.writerow({'k':row['k'],'g':row['g'],'d':row['d'],
            'physical_integer':row['n'],**{m:float(v[i]) for m,v in models.items()}})
report={
    'protocol':{
        'purpose':'Arithmetic trend discovery on analytically known extremal integer displacements, not unknown interior',
        'origin':'unrestricted integral GV charge, minimum height n=0',
        'train_k':'4..11','test_k':'12..19','test_genus':[int(r['g']) for r in te],
        'synthetic_truth_warning':'All extremal labels are generated from the same independent geometric theorem; model success is rediscovery and extrapolation of that theorem, not new GV computation.',
    },
    'learned':{'sign_period':p,'sign_map_by_k_mod4':signs,'minimum_polynomial_degree':degree,
        'absolute_magnitude_polynomial':str(P),'rotary_phase_angle':'pi*k/2',
        'fourier_identity_on_integers':'cos(pi*k/2)+sin(pi*k/2) equals the period-four sign pattern',
        'phase_coefficients':[float(q) for q in phase_coef]},
    'preholdout_phase_selection_top5':phase_candidates[:5],
    'metrics':metrics,
    'heldout':[{**{k:row[k] for k in ('k','g','d','n')},'symbolic_prediction':int(pr),
                 'rotary_phase_prediction':float(pa)}
               for row,pr,pa in zip(te,pred_symbolic,phase_pred)]
}
(ROOT/'results'/'model_benchmarks.json').write_text(json.dumps(report,indent=2)+'\n')
print('learned period',p,'sign_map',signs,'abs_poly',P)
for k,v in metrics.items():print(k,v)
print('phase selection top5:',phase_candidates[:5])
