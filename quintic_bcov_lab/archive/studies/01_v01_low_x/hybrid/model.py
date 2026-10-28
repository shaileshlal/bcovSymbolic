"""Quintic low-X predictive model. Uses exact ground truth, no test-set labels in fitting.

The model predicts (a_{g,1},...,a_{g,4})/a_{g,0}, NOT exact rationals.
It fits log-series cumulants per genus using a small inverse-genus feature dictionary.
Physical boundary fixing and exact rational reconstruction are separate stages.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
import gzip, json, math, sys
import numpy as np

sys.set_int_max_str_digits(0)


def log_jet(normalized: list[Q]):
    """Logarithm of polynomial normalized[0]=1, through X^4, exactly."""
    assert normalized[0] == 1
    c = [Q(0)] * len(normalized)
    for n in range(1, len(c)):
        c[n] = normalized[n] - sum(
            (Q(k,n)*c[k]*normalized[n-k] for k in range(1,n)), Q(0))
    return c


def exp_jet(log_coeff: np.ndarray):
    """Floating predicted log-series to polynomial series (not rational reconstruction)."""
    log_coeff=np.asarray(log_coeff,dtype=float)
    if log_coeff.ndim==1:log_coeff=log_coeff[None,:]
    out=np.zeros((len(log_coeff),log_coeff.shape[1]+1));out[:,0]=1
    for n in range(1,out.shape[1]):
        for k in range(1,n+1):
            out[:,n]+=k*log_coeff[:,k-1]*out[:,n-k]/n
    return out


def load_exact_dataset(path):
    with gzip.open(path,'rt') as f: rows=json.load(f)['rows']
    gs=[]; logs=[]; a0=[]; exact=[]
    for row in rows:
        g=int(row['genus']); aa=tuple(map(Q,row['coefficients']))
        if g<5 or len(aa)!=5:continue
        c=log_jet([Q(1)]+[x/aa[0] for x in aa[1:]])
        gs.append(g);logs.append([float(c[j]/g) for j in range(1,5)])
        a0.append(str(aa[0]));exact.append(tuple(map(str,aa)))
    return dict(g=np.asarray(gs),cum=np.asarray(logs),a0=a0,exact=exact)


@dataclass
class InverseGenusRegressor:
    degree: int = 2
    min_train_genus: int = 40
    max_train_genus: int = 53
    coefficients: np.ndarray | None = None

    def basis(self,genus):
        g=np.asarray(genus,dtype=float).reshape(-1)
        if np.any(g<=0):raise ValueError('g must be positive')
        x=1./g
        return np.stack([x**d for d in range(self.degree+1)],axis=-1)

    def fit(self,data):
        gg=data['g'];sel=(gg>=self.min_train_genus)&(gg<=self.max_train_genus)
        if np.count_nonzero(sel)<2*(self.degree+1):raise ValueError('Too few training genera')
        # Joint multi-output LS; reuse the same design matrix for four log coefficients.
        self.coefficients=np.linalg.lstsq(self.basis(gg[sel]),data['cum'][sel],rcond=None)[0]
        return self

    def predict_log_jet(self,genus):
        g=np.asarray(genus,dtype=float).reshape(-1)
        if self.coefficients is None:raise ValueError('Model not fitted')
        return (self.basis(g)@self.coefficients)*g[:,None]

    def predict_ratios(self,genus):return exp_jet(self.predict_log_jet(genus))

    def to_json(self):
        return {'algorithm':'inverse-genus log-series multioutput least squares',
                'degree':self.degree,'training_window':[self.min_train_genus,self.max_train_genus],
                'coefficients':self.coefficients.tolist(),
                'output':'a[g,j]/a[g,0], approximate floating point, j=0..4',
                'normalization':'a[g,0] is known independently from the Bernoulli constant map'}


def score(model,data,start,end):
    gg=data['g'];mask=(gg>=start)&(gg<=end)
    true=exp_jet(data['cum'][mask]*gg[mask,None])[:,1:]
    estimated=model.predict_ratios(gg[mask])[:,1:]
    # This evaluates coefficient ratios, NOT exact rational equality.
    error=np.abs(estimated/true-1)
    b=np.array([[(-1)**j*math.comb(int(g),j) for j in range(1,5)] for g in gg[mask]],dtype=float)
    base=np.abs(b/true-1)
    def m(x):return {'mean_relative_error':float(x.mean()),'median_relative_error':float(np.median(x)),
                     'p90_relative_error':float(np.quantile(x,.90)),'maximum_relative_error':float(x.max()),
                     'by_j_median':[float(np.median(x[:,j])) for j in range(4)]}
    return {'genus_range':[start,end],'n_genera':int(sum(mask)),'model':m(error),'binomial_baseline':m(base)}


def train_and_evaluate(data_file,out_dir,train_start=40,train_end=53,degree=2):
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    data=load_exact_dataset(data_file)
    model=InverseGenusRegressor(degree,train_start,train_end).fit(data)
    summary={'training':model.to_json(),
      'held_out': [score(model,data,a,b) for a,b in [(54,100),(101,250),(251,1000),(54,1000)]],
      'limits': ['Reference exact low-X dataset is derived from BCOV + exact enumerative vanishing conditions',
                 'No g>53 coefficient labels enter fitting',
                 'Numerical extrapolation is not exact rational reconstruction',
                 'Window and order selected exploratorily; evaluate rolling splits before claiming out-of-sample universality']}
    (out_dir/'model.json').write_text(json.dumps(model.to_json(),indent=2))
    (out_dir/'benchmark.json').write_text(json.dumps(summary,indent=2))
    import csv
    with (out_dir/'heldout_predictions.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['genus','j','predicted_ratio','true_ratio','relative_error'])
        for g in [54,60,80,100,200,500,1000]:
            idx=np.where(data['g']==g)[0][0]
            true=exp_jet(data['cum'][idx:idx+1]*g)[0]
            pred=model.predict_ratios([g])[0]
            for j in range(1,5):wr.writerow([g,j,format(pred[j],'.17g'),format(true[j],'.17g'),format(abs(pred[j]/true[j]-1),'.17g')])
    return summary

def score_physics_anchored(model,data,max_genus=1000):
    """Conditional audit: feed the *known* exact a1 as input; predict a2..a4.

    This benchmark supplies the first cumulant from the all-genus Bernoulli
    theorem. The cached value is verified independently at selected genera;
    it is NOT a blind a1 prediction. For new g compute this via leading_only.
    """
    gg=data['g'];mask=(gg>model.max_train_genus)&(gg<=max_genus)
    g=gg[mask];true_c=data['cum'][mask]*g[:,None]
    pred_c=model.predict_log_jet(g)
    pred_c[:,0]=true_c[:,0]  # exact known analytic input, not a model output
    true=exp_jet(true_c)[:,2:]; predicted=exp_jet(pred_c)[:,2:]
    err=np.abs(predicted/true-1)
    return {'predicted_columns':[2,3,4], 'supplied_exact_column':1,
            'genus_range':[int(g[0]),int(g[-1])],
            'median_relative_error':float(np.median(err)),
            'mean_relative_error':float(np.mean(err)),
            'p90_relative_error':float(np.quantile(err,.9)),
            'max_relative_error':float(np.max(err)),
            'by_column_median':[float(np.median(err[:,j])) for j in range(3)],
            'input_a1_policy':'Supplied known all-genus Bernoulli coefficient. Not counted as predicted.'}

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--data',required=True);parser.add_argument('--output',default='results/hybrid')
    parser.add_argument('--train-start',type=int,default=40);parser.add_argument('--train-end',type=int,default=53)
    args=parser.parse_args()
    report=train_and_evaluate(args.data,args.output,args.train_start,args.train_end)
    data=load_exact_dataset(args.data)
    model=InverseGenusRegressor(2,args.train_start,args.train_end).fit(data)
    report['physics_anchored']=score_physics_anchored(model,data)
    from pathlib import Path
    (Path(args.output)/'benchmark.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'raw':report['held_out'],'physics_anchored':report['physics_anchored']},indent=2))


