"""Independently recover all a_(g,n), n>=g, by a binomial triangular transform.

No unknown pure-X coefficient enters the prediction. Uses the known anomaly
body, local conifold frame, and universal gap. Compare only after prediction.
"""
import argparse
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json
import sys


def verify(checkpoint,engine,output):
    sys.path.insert(0,str(Path(engine).resolve()))
    from bcov_yy import Poly
    from exact_series import evaluate_poly
    from period_frames import frame
    from physical_boundaries import conifold_coefficient
    reports=[]
    for a in json.loads(Path(checkpoint).read_text())['amplitudes']:
        g=a['genus'];d=3*g-3
        values,k,tau=frame('conifold',6*g+16)
        prescribed=conifold_coefficient(g)*tau**(2-2*g)/k**(g-1)
        body=evaluate_poly(Poly.from_data(a['anomaly']),values)
        R=prescribed-body
        predicted={}
        for j in range(d,g-1,-1):
            predicted[j]=(-1)**j*R.coeff(-j)-sum(
                (comb(n,j)*predicted[n] for n in range(j+1,d+1)),Q(0))
        actual=[Q(x) for x in a['ambiguity']]
        assert all(v==actual[n] for n,v in predicted.items()),g
        reports.append({'genus':g,'coefficients_predicted':len(predicted),
                        'predictions':{n:str(v) for n,v in sorted(predicted.items())},'all_match':True})
        print(f'g={g}: {len(predicted)} conifold-binomial predictions verified',flush=True)
    Path(output).write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',default='physical/amplitudes.json')
    p.add_argument('--engine',default='engine')
    p.add_argument('--output',default='results/triangular_checks.json')
    a=p.parse_args();verify(a.checkpoint,a.engine,a.output)
