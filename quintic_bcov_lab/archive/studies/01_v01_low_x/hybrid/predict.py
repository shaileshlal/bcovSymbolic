"""Produce new normalized-X forecasts (not exact rationals).

Example: python -m hybrid.predict --genus 1001 --train-end 1000
Trains only on already-computed lower genera and supplies all-genus a0/a1.
"""
import argparse,sys,json,math
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from .model import InverseGenusRegressor,exp_jet
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from leading_only import bernoulli_abs,constant_coefficients,linear_coefficient
sys.set_int_max_str_digits(0)

def forecast(genus,train_end,data_file,training_width=14):
    if genus<=train_end:raise ValueError('Forecast genus must be > training end')
    with np.load(data_file) as data:
        g=data['g'];c=data['cum']
    data={'g':g,'cum':c}
    model=InverseGenusRegressor(2,train_end-training_width+1,train_end).fit(data)
    cumulants=model.predict_log_jet([genus]);
    bb=bernoulli_abs(genus);c0=constant_coefficients(bb)[genus]
    c1=linear_coefficient(genus,bb,constant_coefficients(bb))
    cumulants[0,0]=float(c1/c0)
    values=exp_jet(cumulants)[0]
    return {'genus':genus, 'fit_genera':[model.min_train_genus,train_end],
            'a0_exact':str(c0),'a1_exact':str(c1),
            'a2_a0_predicted':float(values[2]),'a3_a0_predicted':float(values[3]),
            'a4_a0_predicted':float(values[4]),
            'cumulants_predicted':[float(v) for v in cumulants[0]],
            'warning':'a2,a3,a4 are numerical approximations, not exact rationals; no novel coefficient has been checked by BCOV.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--genus',type=int,default=1001)
    p.add_argument('--train-end',type=int,default=1000)
    p.add_argument('--cached-data',default=str(Path(__file__).resolve().parent/'results'/'float_dataset.npz'))
    p.add_argument('--out',default=str(Path(__file__).resolve().parent/'results'/'forecast_g1001.json'))
    a=p.parse_args()
    result=forecast(a.genus,a.train_end,a.cached_data)
    Path(a.out).write_text(json.dumps(result,indent=2))
    print({k:v for k,v in result.items() if k not in ('a0_exact','a1_exact')})
