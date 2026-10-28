#!/usr/bin/env python3
"""Small-data, chronology-safe quintic intermediate ambiguity benchmarks.

All labels originate in *known reference amplitudes*; never use YY physical
pivots, boundary anchors, or other target genus ambiguity coefficients as features.
The mixed-YY data at held-out genus presuppose physically chosen lower genera.
"""
from __future__ import annotations
import argparse,csv,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor,ExtraTreesClassifier
from sklearn.linear_model import LinearRegression
from sklearn.dummy import DummyClassifier

ROOT=Path(__file__).resolve().parent
RAW=ROOT/'data'/'intermediate_training.csv'

def log_abs_exact(s:str)->float:
    f=abs(Fraction(s))
    if not f:return float('-inf')
    def lg(n):
        x=str(n)
        return len(x)-1+math.log10(int(x[:15])/10**(len(x[:15])-1))
    return lg(f.numerator)-lg(f.denominator)

def load_rows(path=RAW):
    with open(path,newline='') as fp: a=list(csv.DictReader(fp))
    rows=[]
    for x in a:
        q=Fraction(x['c_X_kernel']);b=Fraction(x['b_x_bernstein']); g=int(x['genus']);r=int(x['residual_index'])
        assert b==(-1)**(g-2-r)*q
        y=log_abs_exact(x['c_X_kernel'])
        assert abs(y-abs(float(x['label_signed_log10'])))<1e-6
        z={k:float(x[k]) for k in ('mixed_log10_max_at_Xj','mixed_terms_at_Xj','mixed_v1_terms_at_Xj','mixed_v2_terms_at_Xj','mixed_v3_terms_at_Xj','relative_pivot_degree')}
        z.update({'genus':g,'r':r,'j':int(x['X_pivot_degree']),'D':int(x['D']),
                  'sign':int(x['label_sign']), 'log_mag':y,'signed_log':float(x['label_signed_log10']),
                  'M':float(x['mixed_log10_max_at_Xj']), 'exact_c':x['c_X_kernel'], 'exact_b':x['b_x_bernstein']})
        assert z['sign']==(1 if q>0 else -1)
        rows.append(z)
    assert len(rows)==22
    return rows

def feat(rows,tag):
    g=np.array([v['genus'] for v in rows],float)
    r=np.array([v['r'] for v in rows],float)
    M=np.array([v['M'] for v in rows],float)
    if tag=='genus':return np.column_stack([g,r])
    if tag=='mixed':return np.column_stack([M,r])
    if tag=='mixed_plus_genus':return np.column_stack([M,r,g])
    if tag=='mixed_residual':return np.column_stack([g,r])
    if tag=='v01':return np.column_stack([g,r, [v['j'] for v in rows], [3*v['genus']-3 for v in rows], [v['D'] for v in rows],
       *([v[k] for v in rows] for k in ('mixed_terms_at_Xj','mixed_v1_terms_at_Xj','mixed_v2_terms_at_Xj','mixed_v3_terms_at_Xj','M'))])
    raise ValueError(tag)

def fit_magnitude(train,test,tag):
    tr_y=np.array([a['log_mag'] for a in train]); X=feat(train,tag);Xt=feat(test,tag)
    if tag=='mixed_residual':
        mdl=LinearRegression().fit(X,tr_y-np.array([a['M'] for a in train]))
        pred=mdl.predict(Xt)+np.array([a['M'] for a in test])
    else:
        mdl=LinearRegression().fit(X,tr_y)
        pred=mdl.predict(Xt)
    return pred, {'intercept':float(mdl.intercept_),'coefficients':list(map(float,mdl.coef_))}

def fit_sign(train,test,tag):
    if tag=='majority':
        mdl=DummyClassifier(strategy='most_frequent').fit(feat(train,'genus'),[a['sign'] for a in train])
        return mdl.predict(feat(test,'genus'))
    if tag=='v01_trees':
        mdl=ExtraTreesClassifier(n_estimators=180,max_features=.8,min_samples_leaf=2,random_state=7,n_jobs=1).fit(feat(train,'v01'),[a['sign'] for a in train]);return mdl.predict(feat(test,'v01'))
    if tag=='mod3_diagnostic':
        # Not prespecified, hypothesis suggested by inspecting five training genera.
        votes={}
        for a in train:votes.setdefault(a['genus']%3,[]).append(a['sign'])
        fallback=1 if sum(a['sign'] for a in train)>=0 else -1
        return np.array([1 if sum(votes.get(a['genus']%3,[fallback]))>=0 else -1 for a in test])
    raise ValueError(tag)

def summarize_predictions(agg):
    v=np.array([abs(x['pred_log_mag']-x['actual_log_mag']) for x in agg])
    sg=np.array([x['pred_sign']==x['actual_sign'] for x in agg])
    sh=np.array([x['pred_r'] is not None and x['pred_r'] < x['train_max_D'] for x in agg])
    return {'n_coefficients':len(agg),'n_genera':len(set(x['genus'] for x in agg)),
       'MAE_log10_magnitude':float(v.mean()),'median_abs_log10_magnitude':float(np.median(v)),
       'p90_abs_log10_magnitude':float(np.quantile(v,.90)),
       'median_magnitude_factor':float(10**np.median(v)),
       'sign_accuracy_by_coefficient':float(np.mean(sg)),
       'sign_accuracy_by_genus':float(np.mean([next(x for x in agg if x['genus']==g)['pred_sign']==next(x for x in agg if x['genus']==g)['actual_sign'] for g in sorted(set(x['genus'] for x in agg))])),
       'MAE_unseen_position':float(v[~sh].mean()) if (~sh).any() else None,
       'n_unseen_position':int((~sh).sum())}

def run(rows, upto,mag_tag,sign_tag, rolling=False):
    buckets = [( [a for a in rows if a['genus']<g], [a for a in rows if a['genus']==g]) for g in range(9,13)] if rolling else [([a for a in rows if a['genus']<=upto],[a for a in rows if a['genus']>upto])]
    out=[]; params=[]
    for train,test in buckets:
        if not train or not test:continue
        if mag_tag=='v01_trees':
            m=ExtraTreesRegressor(n_estimators=180,max_features=.8,min_samples_leaf=2,random_state=7,n_jobs=1).fit(feat(train,'v01'),np.array([a['signed_log'] for a in train])); p=m.predict(feat(test,'v01'))
            # Predict sign and log mag from a single signed regressor as in v0.1
            sign=np.where(p>=0,1,-1)
            mag=np.abs(p)
            pars={}
        else:
            mag,pars=fit_magnitude(train,test,mag_tag)
            sign=fit_sign(train,test,sign_tag)
        params.append({'training_max_genus':max(t['genus'] for t in train),'parameters':pars})
        maxD=max(v['D'] for v in train)
        for a,m,s in zip(test,mag,sign):
            out.append({'genus':a['genus'],'r':a['r'],'pred_r':a['r'],'train_max_D':maxD,
             'actual_sign':a['sign'],'pred_sign':int(s),
             'actual_log_mag':a['log_mag'],'pred_log_mag':float(m),
             'error_dex':abs(a['log_mag']-float(m)),
             'true_exact_c':a['exact_c'],
             'training_max_genus':max(t['genus'] for t in train)})
    return out,params

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default=str(RAW));ap.add_argument('--output',default=str(ROOT/'results'))
    args=ap.parse_args();outdir=Path(args.output);outdir.mkdir(parents=True,exist_ok=True)
    rows=load_rows(Path(args.input))
    # Pre-specified magnitudes. Sign predictions are deliberately separate.
    maglist=['genus','mixed','mixed_plus_genus','mixed_residual','v01_trees']
    experiments={}
    all_preds=[]
    for split in ['train_to_8','train_to_9','expanding_9_to_12']:
        x={}
        for mag in maglist:
            signs=['combined_v01'] if mag=='v01_trees' else ['majority','v01_trees','mod3_diagnostic']
            for s in signs:
                pred,params=run(rows,8 if split!='train_to_9' else 9,mag,s, rolling=split.startswith('expanding'))
                name=mag+'__'+s
                x[name]={'scores':summarize_predictions(pred),'fit_parameters':params}
                if s in ('majority','combined_v01','mod3_diagnostic'):
                    all_preds += [{'split':split,'model':name,**a} for a in pred]
        experiments[split]=x
    # Oracle-sign upper bound on signed-log MAE: magnitude error only.
    # Shape-only within-genus oracle c0: conditional diagnostic, never used in blind model.
    by_g={}
    for a in rows:by_g.setdefault(a['genus'],[]).append(a)
    shapes=[]
    for g,group in by_g.items():
        c0=next(a['log_mag'] for a in group if a['r']==0)
        for a in group:
            if a['r']>0: shapes.append({'genus':g,'r':a['r'],'log10_abs_c_r_over_c0':a['log_mag']-c0})
    shape_tests={}
    for cutoff in (8,9):
        tr=[a for a in shapes if a['genus']<=cutoff]
        te=[a for a in shapes if a['genus']>cutoff]
        t=np.array([v['log10_abs_c_r_over_c0'] for v in tr]);r=np.array([v['r'] for v in tr])
        slope=float(np.dot(t,r)/np.dot(r,r))
        predictions=[{'genus':v['genus'],'r':v['r'],'true_shape_log10':v['log10_abs_c_r_over_c0'],
                     'pred_shape_log10':slope*v['r'],
                     'abs_error_dex':abs(slope*v['r']-v['log10_abs_c_r_over_c0'])} for v in te]
        errors=[z['abs_error_dex'] for z in predictions]
        shape_tests[str(cutoff)]={'n_train_shape':len(tr),'n_test_shape':len(te),
             'slope_per_index':slope,'MAE_dex':float(np.mean(errors)),
             'median_abs_error_dex':float(np.median(errors)),
             'requires_oracle_c0_at_test_genus':True,'predictions':predictions}
    report={'schema':'quintic-yy-restricted-ml-v1','provenance':'reference derived canonical YY residual labels, mixed terms from exact anomaly',
      'features_allowed':['genus','residual index','kernel dimensions','mixed anomaly features at same X degree'],
      'forbidden_features':['physical YY ambiguity pivot coefficients','reference-derived anchor','reference-derived midpoint value','target-genus GV invariant'],
      'conditioning':'held-out genus anomaly is known, hence its lower-genus physical ambiguities are available; not a sequential bootstrap',
      'dataset':{'n_reference_genera':11,'n_nonzero_genera':9,'n_coefficients':22,'g_range':[4,12],
         'sign_uniform_within_genus':all(len({a['sign'] for a in rr})==1 for rr in by_g.values())},
      'experiments':experiments,'shape_log_ratios':shapes,'shape_tests':shape_tests,
      'warning':'Sign period-three is a post-hoc exploratory diagnostic, not an established law; only four held-out genera.'}
    (outdir/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
    with (outdir/'predictions.csv').open('w',newline='') as fp:
        keys=list(all_preds[0]);w=csv.DictWriter(fp,fieldnames=keys);w.writeheader();w.writerows(all_preds)
    print('N',len(rows));
    for split in experiments:
        print('\n',split)
        for m in ['genus__majority','mixed__majority','mixed_plus_genus__majority','mixed_residual__majority','v01_trees__combined_v01','mixed__mod3_diagnostic']:
            z=experiments[split][m]['scores'];print(m, f"MAE={z['MAE_log10_magnitude']:.3f}dex median={z['median_abs_log10_magnitude']:.3f} factor={z['median_magnitude_factor']:.2f} sign={z['sign_accuracy_by_coefficient']:.3f} signG={z['sign_accuracy_by_genus']:.3f}")
if __name__=='__main__': main()
