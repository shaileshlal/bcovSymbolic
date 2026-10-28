"""Small-sample baseline: map mixed antiholomorphic polynomial features to a_g(X).

The genus holdout is strictly chronological. This experiment is a diagnostic,
not evidence that BCOV determines the ambiguity without boundary data.
"""
from fractions import Fraction as Q
import json,math,numpy as np,gzip
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.dummy import DummyRegressor

def signed_log(q):
    q=Q(q)
    if q==0:return 0.
    return math.copysign(math.log10(1+abs(float(q))),float(q))

def rows_from_checkpoint(path,airy_data=None):
    data=json.load(open(path))['amplitudes']
    alphas={}
    if airy_data:
        with gzip.open(airy_data,'rt') as fp:
            alphas={int(g):Q(v['alpha']) for g,v in json.load(fp).items()}
    x=[];y=[];meta=[]
    for amp in data:
        g=amp['genus'];a=list(map(Q,amp['ambiguity']));a0=a[0]
        poly=amp['anomaly']
        buckets={}
        for m,n,d in poly:
            buckets.setdefault(m[4],[]).append((m,signed_log(Q(int(n),int(d))/a0)))
        maxj=len(a)-1
        for j in range(1,maxj+1):
            features=[g,j,j/(3*g-3),3*g-3-j,len(poly)]
            if airy_data:features.extend([signed_log(alphas[g]/a0),signed_log(alphas[g])])
            for dj in [-2,-1,0,1,2]:
                bucket=buckets.get(j+dj,[])
                values=[v for _,v in bucket]
                features.extend([len(bucket),sum(values),sum(map(abs,values)),max(map(abs,values),default=0)])
                for v in [1,2,3]:
                    vv=[c for m,c in bucket if m[v]]
                    features.extend([len(vv),sum(vv)])
            x.append(features);y.append(signed_log(a[j]/a0));meta.append((g,j,str(a[j])))
    return np.array(x),np.array(y),meta

def benchmark(path,airy_data=None):
    x,y,meta=rows_from_checkpoint(path,airy_data)
    gg=np.array([g for g,j,a in meta]);tr=(gg<=8);te=(gg>=9)
    model=ExtraTreesRegressor(n_estimators=180,max_features=.8,min_samples_leaf=2,random_state=7,n_jobs=1)
    model.fit(x[tr],y[tr]);base=DummyRegressor(strategy='mean').fit(x[tr],y[tr])
    p=model.predict(x[te]);q=base.predict(x[te]);real=y[te]
    return {'airy_features':bool(airy_data),'train_genera':[2,8],'test_genera':[9,12],
            'train_coefficients':int(tr.sum()),'test_coefficients':int(te.sum()),
            'target':'signed log10(1+|a_gj/a_g0|)',
            'model_mae':float(np.abs(p-real).mean()),'constant_baseline_mae':float(np.abs(q-real).mean()),
            'model_sign_accuracy':float(np.mean(np.sign(p)==np.sign(real))),
            'constant_sign_accuracy':float(np.mean(np.sign(q)==np.sign(real))),
            'important_note':'The held-out anomaly polynomial already contains lower-genus physically fixed data; this is a conditional, not a generative BCOV predictor.'}

def experimental_full_projection(path,target_genus=9):
    """Predict all pure-X slots from mixed polynomial, then solve exact partial boundary.

    Uses only g<=8 amplitude labels to train. A constant-map a0 is permitted.
    All other target genus ambiguity coefficients are held out.
    This is intended primarily as a stringent negative control.
    """
    from .boundary import boundary_system
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'engine'))
    from bcov_yy import Poly
    from leading_only import bernoulli_abs,constant_coefficients
    x,y,meta=rows_from_checkpoint(path)
    g=np.array([z[0] for z in meta]);tr=(g<=8);test=g==target_genus
    if target_genus<=8:raise ValueError('use a held-out target genus >=9')
    mdl=ExtraTreesRegressor(n_estimators=180,max_features=.8,min_samples_leaf=2,random_state=7,n_jobs=1)
    mdl.fit(x[tr],y[tr]);pr=mdl.predict(x[test]);indexes=[j for gg,j,a in meta if gg==target_genus]
    a0=constant_coefficients(bernoulli_abs(target_genus))[target_genus]
    predicted=[a0]
    for value in pr:
        z=math.copysign(math.expm1(abs(float(value))*math.log(10.)),value)
        predicted.append(Q(str(z))*a0)
    if indexes!=list(range(1,3*target_genus-2)):raise ValueError('missing sequence indices')
    amps={a['genus']:a for a in json.load(open(path))['amplitudes']}
    truth=tuple(map(Q,amps[target_genus]['ambiguity']))
    affine=boundary_system(target_genus,Poly.from_data(amps[target_genus]['anomaly']),mum_degrees=[0])
    candidate=affine.project_by_free_coordinates(predicted)
    assert affine.verify(candidate)
    # Relative error with numerically huge rational pivots measured using normalized ratios.
    errors=[abs(float((c-t)/t)) for c,t in zip(candidate[1:],truth[1:]) if t]
    return {'target_genus':target_genus,'train_genera':[2,8],
            'constraint_type':'conifold + orbifold + MUM constant map, not degree>=1 GV',
            'partial_rank':affine.rank,'free_slots':list(affine.free_indices),
            'relative_error_median':float(np.median(errors)),
            'relative_error_max':float(np.max(errors)),
            'exact_true_coefficients':sum(int(a==b) for a,b in zip(candidate,truth)),
            'total_coefficients':len(truth),
            'verified_exact_boundary_constraints':True}
