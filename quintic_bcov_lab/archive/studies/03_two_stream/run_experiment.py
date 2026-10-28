"""Strict forward-genus two-stream YY/GV ablation.

Truth: audited exact quintic amplitudes g2..12. Only lower-genus reference
information (g<=4) is used to populate the GV charge stream. Genus 10..12
are never used to select predictor hyperparameters. All q-series/MUM arithmetic
is performed over fractions. Predictors are deliberately small, ridge-based
surrogates to avoid pretending a six-genus dataset fits a deep transformer.
"""
from __future__ import annotations
import json,sys,math
from pathlib import Path
from fractions import Fraction as F
from math import comb
import numpy as np
P=Path(__file__).resolve().parent
SRC=P/'quintic_ai_predictor'
sys.path.insert(0,str(SRC/'engine'))
sys.path.insert(0,str(P))
from bcov_yy import Poly
from period_frames import frame
from exact_series import evaluate_poly
from physical_boundaries import seed_gv,record_gv,gv_kernel,gw_from_gv,q_coefficients,constant_map
from high_degree_phase_study import reverse_cumulants

DATA={x['genus']:x for x in json.loads((SRC/'reference_amplitudes.json').read_text())['amplitudes']}
G=sorted(DATA)
A={g:[F(z) for z in DATA[g]['ambiguity']] for g in G}

def signature(q):
    q=F(q)
    if not q:return 0.
    # Prevent overflow: no float of numerator/denominator until separate logs
    return math.copysign(math.log10(abs(q.numerator))-math.log10(q.denominator)+1,q.numerator)

def log_coefficients(values,n):
    """Formal logarithm 1+sum b_i z**i with b0=1."""
    cs=[F(0)]
    for j in range(1,n+1):
        cs.append(values[j]-sum((F(k,j)*cs[k]*values[j-k] for k in range(1,j)),F(0)))
    return cs

def anom_cumulant(g,sector,rmax=6):
    k=3*g-3
    e_1,e_2,e_3=sector
    j0=k-e_1-2*e_2-3*e_3
    coeff={m[4]:F(int(n),int(d)) for m,n,d in DATA[g]['anomaly']
           if m[:4]==[0,e_1,e_2,e_3]}
    top=coeff.get(j0,F(0))
    if not top: return [0.]*(rmax+1)
    b=[F(1)]+[coeff.get(j0-j,F(0))/top for j in range(1,rmax+1)]
    return [float(v)/g for v in log_coefficients(b,rmax)]

def make_gv():
    gv=seed_gv(4)
    for g in [2,3,4]:
        p=Poly.from_data(DATA[g]['anomaly'])+Poly({(0,0,0,0,j):v for j,v in enumerate(A[g])})
        record_gv(g,p,4,gv)
    assert gv[2,4]==534750 and gv[3,4]==8625
    for g in range(5,13):
        for d in range(1,5):gv[g,d]=F(0) # independently known Castelnuovo vanishing
    return gv

def create_data():
    gv=make_gv()
    features={}; physical_checks={}; records={}
    vals,k,_=frame('mum',14)
    for g in range(4,13):
        p_an=Poly.from_data(DATA[g]['anomaly'])
        an_q=q_coefficients(evaluate_poly(p_an,vals)*k**(g-1),4)
        f_actual=Poly({(0,0,0,0,j):v for j,v in enumerate(A[g])})
        full_q=q_coefficients(evaluate_poly(p_an+f_actual,vals)*k**(g-1),4)
        target=[constant_map(g)]+[gw_from_gv(g,d,gv) for d in range(1,5)]
        errors=[full_q[d]-target[d] for d in range(5)]
        assert not any(errors),(g,errors)
        rhs=[target[i]-an_q[i] for i in range(5)]
        # Construct finite character representation of the h>=1 part, using
        # physical BPS spin h, multicover m, and curve degree d. No X-index here.
        phase=[]
        for lam in [0.35/math.sqrt(g),0.8/math.sqrt(g),1.6/math.sqrt(g)]:
            acc=0.
            for d in range(1,5):
                for m in range(1,d+1):
                    if d%m:continue
                    for h in range(1,5):
                        acc+=float(gv[h,d//m])/m*(2*math.sin(m*lam/2))**(2*h-2)
            phase.append(math.copysign(math.log10(1+abs(acc)),acc))
        # One numerically-stable coordinate from the exact MUM physical RHS.
        # These data are known without the held-out high-X coefficients.
        charge=[signature(v) for v in target[1:]]
        boundary=[signature(v) for v in rhs[1:]]
        ann1=anom_cumulant(g,(1,0,0));ann2=anom_cumulant(g,(0,1,0));ann3=anom_cumulant(g,(0,0,1))
        b,c=reverse_cumulants(A[g],6)
        features[g]={'yy1':ann1,'yy2':ann2,'yy3':ann3,'charge':charge,
                     'boundary':boundary,'phases':phase,
                     'target': [float(x)/g for x in c],
                     'b':[float(x) for x in b]}
        physical_checks[g]={'mum_degree_0_to_4_exact':True,
                             'known_gv_constraints':[str(q) for q in target],
                             'low_ambiguity_rhs':[str(q) for q in rhs],
                             'charge':[str(q) for q in charge]}
        records[g]={'genus':g,'high_degree_start':3*g-9,
                    'max_GV_directly_visible_X':4,
                    'high_degree_direct_GV_response_zero':3*g-9>4}
        print('exact GV/MUM check',g, 'good',flush=True)
    return features,physical_checks,records,gv

def stack_x(g,r, mode,fb):
    """One normalized sample per (g, tail distance r)."""
    base=[1.,1./g,1./g**2,1./g**3]
    if mode=='base2': return base[:3]
    if mode=='base3': return base
    if mode=='yy':return base+[fb[g]['yy1'][r]]
    if mode=='gv':return base+[fb[g]['boundary'][3]]
    if mode=='yy+gv':return base+[fb[g]['yy1'][r],fb[g]['boundary'][3]]
    if mode=='yy+phase':return base+[fb[g]['yy1'][r],fb[g]['phases'][1]]
    if mode=='gv+phase':return base+[fb[g]['boundary'][3],fb[g]['phases'][1]]
    if mode=='yy+gv+phase':return base+[fb[g]['yy1'][r],fb[g]['boundary'][3],fb[g]['phases'][1]]
    raise ValueError(mode)

def fit_ridge_predict(mode,r,train,test,fb,alpha=1.0):
    mat=np.array([stack_x(g,r,mode,fb) for g in train],float)
    fut=np.array([stack_x(g,r,mode,fb) for g in test],float)
    y=np.array([fb[g]['target'][r] for g in train],float)
    # Stable unpenalized inverse-genus basis (very accurate physics baseline)
    nbase=3 if mode=='base2' else 4
    base=mat[:,:nbase];fbase=fut[:,:nbase]
    coef=np.linalg.lstsq(base,y,rcond=1e-13)[0]
    pred=fbase@coef
    if mat.shape[1]==nbase:return pred
    rtrain=y-base@coef
    extra=mat[:,nbase:]; extra_future=fut[:,nbase:]
    # Residualize extra channel against genus-polynomial features.
    # Result measures additional cross-genus information conditional on baseline.
    w=np.linalg.lstsq(base,extra,rcond=1e-13)[0]
    ex=extra-base@w
    exf=extra_future-fbase@w
    scales=np.sqrt(np.mean(ex*ex,axis=0))
    good=scales>1e-11
    if not np.any(good):return pred
    ex=ex[:,good]/scales[good]
    exf=exf[:,good]/scales[good]
    update=np.linalg.solve(ex.T@ex+alpha*np.eye(ex.shape[1]),ex.T@rtrain)
    return pred+exf@update

def pred_tail(mode,train,test,fb,alpha=1.0):
    C={g:[0.]*7 for g in test}
    for r in range(1,7):
        pp=fit_ridge_predict(mode,r,train,test,fb,alpha)
        for g,p in zip(test,pp): C[g][r]=p*g
    out={}
    for g,c in C.items():
        b=[1.]
        for j in range(1,7):b.append(sum((k*c[k]*b[j-k]/j for k in range(1,j+1))))
        out[g]=b
    return out

def eval_preds(p,fb,test):
    err=[abs(p[g][j]/fb[g]['b'][j]-1) if fb[g]['b'][j] else float('inf')
         for g in test for j in range(1,7)]
    return {'median_relative_error':float(np.median(err)),
            'mean_relative_error':float(np.mean(err)),
            'p90_relative_error':float(np.quantile(err,.9)),
            'max_relative_error':float(np.max(err)),
            'all_relative_errors':[float(z) for z in err]}

def run():
    fb,checks,records,gv=create_data()
    modes=['base2','base3','yy','gv','yy+gv','yy+phase','gv+phase','yy+gv+phase']
    alphas=[0.001,0.01,0.1,1.,10.,100.]
    # Tune every method ONLY on earlier future-genus validation g8,9,
    # then refit from g4..9, hold out g10..12 entirely.
    tuning={}
    for mode in modes:
        scores={}
        for alpha in (alphas if '+' in mode or mode in ('yy','gv') else [1.]):
            v=[]
            for end in [7,8]:
                pr=pred_tail(mode,list(range(4,end+1)),[end+1],fb,alpha)
                v.extend(eval_preds(pr,fb,[end+1])['all_relative_errors'])
            scores[str(alpha)]=float(np.median(v))
        best=min(scores,key=scores.get)
        tuning[mode]={'best_alpha':float(best),'val_median_relative_error':scores[best], 'candidates':scores}
    out={}
    for mode in modes:
        p=pred_tail(mode,list(range(4,10)),[10,11,12],fb,tuning[mode]['best_alpha'])
        out[mode]={**eval_preds(p,fb,[10,11,12]),
                   'by_genus':{str(g):eval_preds(p,fb,[g])['median_relative_error'] for g in [10,11,12]},
                   'forecast_b':{str(g):[float(x) for x in p[g]] for g in [10,11,12]}}
    # Explicit physical character identity check with +/- Fourier phases.
    phase_check=[]
    for h in range(1,5):
        for m in range(1,5):
            for lam in (0.17,0.61,1.03):
                z=sum(((-1)**(h-1+j))*comb(2*h-2,j)*np.exp(1j*m*(h-1-j)*lam) for j in range(2*h-1))
                phase_check.append(float(abs(z-(2*np.sin(m*lam/2))**(2*h-2))))
    results={'reference_genera':G,'training_genera':list(range(4,10)),
             'tuning_validation_genera':[8,9], 'strict_heldout_genera':[10,11,12],
             'target':'b[g,r]=a[g,K-r]/a[g,K] for r=1..6; top a[g,K] is an oracle in all models',
             'features':'YY: reversed cumulant of v1X^j terms; GV: exact BPS numbers of h<=4,d<=4, exact MUM/GV RHS, GV spin sine characters on physical m,h,d (no YY X-index substitution)',
             'tuning':tuning,'results':out,'physical_verification':checks,
             'visible_degree_audit':records, 'character_identity_max_error':max(phase_check),
             'gv_known':{f'{g},{d}':str(v) for (g,d),v in sorted(gv.items()) if g<=4}}
    (P/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    print('\nSTRICT HELDOUT:',[10,11,12])
    for m in modes:
        print(f"{m:16} median={out[m]['median_relative_error']:.8g} p90={out[m]['p90_relative_error']:.8g} val={tuning[m]['val_median_relative_error']:.8g} alpha={tuning[m]['best_alpha']}")
    print('GV/Fourier character identity max error',max(phase_check))
    print('wrote',P/'results.json')

if __name__=='__main__':run()
