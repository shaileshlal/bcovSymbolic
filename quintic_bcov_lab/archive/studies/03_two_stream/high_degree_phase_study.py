"""Exact high-X cumulants and strict forward-genus ablations (Python/NumPy)."""
from __future__ import annotations
import json, math, argparse
from fractions import Fraction
from pathlib import Path
import numpy as np
from math import comb
ROOT=Path(__file__).resolve().parent
REF=ROOT/'quintic_ai_predictor/reference_amplitudes.json'

def reference():
    dat=json.loads(REF.read_text())['amplitudes']
    return {row['genus']: [Fraction(q) for q in row['ambiguity']] for row in dat}

def reverse_cumulants(a, n=6):
    K=len(a)-1
    b=[a[K-j]/a[K] for j in range(min(K,n)+1)]
    c=[Fraction(0)]
    for i in range(1,len(b)):
        c.append(b[i]-sum((Fraction(k,i)*c[k]*b[i-k] for k in range(1,i)),Fraction(0)))
    check=[Fraction(1)]
    for i in range(1,len(b)):
        check.append(sum((Fraction(k,i)*c[k]*check[i-k] for k in range(1,i+1)),Fraction(0)))
    assert b==check
    return b,c

def gv_phase(h,k,lam):
    """Exact GV character as a finite sum of complex phases, h>=0."""
    return sum(((-1)**(h+r))*comb(2*h,r)*np.exp(1j*k*(h-r)*lam)
               for r in range(2*h+1))

def gv_character(h,k,lam):
    return (2*np.sin(k*lam/2))**(2*h)

def feature_basis(js,name):
    js=np.asarray(js,dtype=float)
    if name=='onehot':
        return np.eye(6)[js.astype(int)-1]
    if name=='ordinary_polynomial':
        x=js/6.
        return np.stack([np.ones_like(x),x,x*x,x*x*x],axis=-1)
    if name=='rope':
        # Frequency bank, rather than a trained attention model: RoPE-like phase features.
        t=2*np.pi*js/13
        return np.stack([np.ones_like(t),np.cos(t),np.sin(t),np.cos(2*t),np.sin(2*t)],axis=-1)
    if name=='rope_full':
        t=2*np.pi*js/13
        return np.stack([np.cos(t),np.sin(t),np.cos(2*t),np.sin(2*t),np.cos(3*t),np.sin(3*t)],axis=-1)
    if name=='gv_character_full':
        t=np.pi*js/9
        return np.stack([(2*np.sin(t/2))**(2*h) for h in range(6)],axis=-1)
    if name=='gv_character':
        # Fixed character dictionary; this is a proposed analogy on X-indices,
        # NOT a derivation identifying X-power with BPS curve degree.
        t=np.pi*js/9
        return np.stack([np.ones_like(t),
              (2*np.sin(t/2))**2,
              (2*np.sin(t/2))**4,
              (2*np.sin(t/2))**6,
              (2*np.sin(t/2))**8],axis=-1)
    raise ValueError(name)

def genus_basis(gg,kind):
    gg=np.asarray(gg,dtype=float)
    if kind=='inverse2':return np.stack([np.ones_like(gg),1/gg,1/gg**2],axis=-1)
    if kind=='inverse3':return np.stack([np.ones_like(gg),1/gg,1/gg**2,1/gg**3],axis=-1)
    if kind=='linear':return np.stack([np.ones_like(gg),gg],axis=-1)
    raise ValueError(kind)

def fit_predict_cumulants(C,train,test,basis='onehot',gbasis='inverse2',ridge=0):
    num=6
    rows=[];ys=[]
    for g in train:
        for j in range(1,num+1):
            rows.append((g,j));ys.append(float(C[g][j]/g))
    xt=np.array([np.outer(genus_basis([g],gbasis)[0],feature_basis([j],basis)[0]).ravel() for g,j in rows])
    yt=np.array(ys)
    rows_test=[(g,j) for g in test for j in range(1,num+1)]
    xs=np.array([np.outer(genus_basis([g],gbasis)[0],feature_basis([j],basis)[0]).ravel() for g,j in rows_test])
    if ridge>0:
        scale=np.maximum(np.linalg.norm(xt,axis=0),1.e-12)
        xt=xt/scale; xs=xs/scale
        w=np.linalg.solve(xt.T@xt+ridge*np.eye(xt.shape[1]),xt.T@yt)
    else:w=np.linalg.lstsq(xt,yt,rcond=1e-12)[0]
    return {(g,j):p*g for (g,j),p in zip(rows_test,xs@w)}

def predict_tail(C,train,test,basis,gbasis,ridge=0):
    cp=fit_predict_cumulants(C,train,test,basis,gbasis,ridge)
    pred={}
    for g in test:
        c=[0]+[cp[g,j] for j in range(1,7)]
        b=[1.]
        for i in range(1,7):b.append(sum((k/i*c[k]*b[i-k] for k in range(1,i+1))))
        pred[g]=b
    return pred

def metrics(pred,B,test):
    errs=[];sign=[]
    for g in test:
        for j in range(1,7):
            v=float(B[g][j]);x=pred[g][j]
            errs.append(abs(x-v)/abs(v) if v else float('inf'))
            sign.append(np.sign(v)==np.sign(x))
    return {'median_relative_error':float(np.median(errs)),
            'mean_relative_error':float(np.mean(errs)),
            'p90_relative_error':float(np.quantile(errs,.90)),
            'max_relative_error':float(np.max(errs)),
            'sign_accuracy':float(np.mean(sign))}

def top_scale_forecast(A,train,test):
    gg=np.array(train,dtype=float)
    yy=np.array([math.log10(float(A[g][-1]/A[g][0])) for g in train])
    # Prespecified two-dimensional model; no test information used for fitting.
    M=np.stack([np.ones_like(gg),gg],axis=1)
    w=np.linalg.lstsq(M,yy,rcond=None)[0]
    return {g:float(10**(w[0]+w[1]*g)) for g in test}

def main():
    A=reference()
    B={};C={}
    for g,a in A.items():
        b,c=reverse_cumulants(a,6)
        B[g]=b;C[g]=c
    train=list(range(6,10));test=[10,11,12]
    config=[('linear (index-aligned)', 'onehot','linear',0),
            ('inverse-g^2 (index-aligned)', 'onehot','inverse2',0),
            ('inverse-g^3 (index-aligned)', 'onehot','inverse3',0),
            ('polynomial index / inverse-g^2','ordinary_polynomial','inverse2',1e-8),
            ('RoPE-like phases / inverse-g^2','rope','inverse2',1e-8),
            ('GV-character phases / inverse-g^2','gv_character','inverse2',1e-8),
            ('RoPE full rank / inverse-g^3','rope_full','inverse3',0),
            ('GV character full rank / inverse-g^3','gv_character_full','inverse3',0),
            ('RoPE-like phases / inverse-g^3','rope','inverse3',1e-8),
            ('GV-character phases / inverse-g^3','gv_character','inverse3',1e-8)]
    output=[]
    for name,index,gshape,ridge in config:
        pred=predict_tail(C,train,test,index,gshape,ridge)
        output.append({'model':name,**metrics(pred,B,test)})
    scale=top_scale_forecast(A,train,test)
    pred=predict_tail(C,train,test,'onehot','inverse3')
    absolute_errors=[];by_genus={}
    for g in test:
        per=[]
        for k in range(7):
            # a_g,top = a_g,0 * alpha; a0 is the exact, separately known constant map.
            observed=float(A[g][-1-k]/A[g][0])
            predicted=scale[g]*pred[g][k]
            err=abs(predicted-observed)/abs(observed) if observed else float('inf')
            absolute_errors.append(err);per.append(err)
        by_genus[str(g)]={'median_relative_error':float(np.median(per)), 'max_relative_error':float(max(per)),
                          'top_scale_relative_error':float(abs(scale[g]/float(A[g][-1]/A[g][0])-1))}
    # g9 exact high tail is uniquely fixed with conifold orbifold MUM constant as independently audited
    phase_checks={}
    for h in range(1,6):
        differences=[abs(gv_phase(h,k,x)-gv_character(h,k,x)) for k in range(1,6)
                     for x in (.12,.39,1.12)]
        phase_checks[str(h)]=float(max(differences))
    result={'schema':'quintic-high-X-phase-v1','reference_genera':sorted(A),
            'train_genera':train,'heldout_genera':test,'tail_coefficients_k':list(range(1,7)),
            'normalization':'b[g,k]=a[g,3g-3-k]/a[g,3g-3]; c[g,k]=[y^k] log sum b[g,k]y^k',
            'predicted_target':'b[g,k] given exact a[g,top] for shape only',
            'models':output,
            'full_absolute_forecast':{'a0_oracle':True,'a_top_predicted':True,'median_relative_error':float(np.median(absolute_errors)),
                                      'p90_relative_error':float(np.quantile(absolute_errors,.9)), 'by_genus':by_genus},
            'gv_fourier_identity_max_errors':phase_checks,
            'note':'GV-character dictionary on X-position is heuristic; no identification of curve degree with X exponent. Small g sample, phase results not proof of optimality.'}
    out=ROOT/'results_high_degree_phase.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print('TRAIN',train,'TEST',test)
    for r in output:print(r['model'], 'median',round(r['median_relative_error'],7),'p90',round(r['p90_relative_error'],7),'sign',r['sign_accuracy'])
    print('FULL using known a0, learned atop and tail',result['full_absolute_forecast'])
    print('GV identity residual',phase_checks)
    print('wrote',out)
if __name__=='__main__':main()
