"""Run from bundle root: python -m hybrid.tests"""
import json,os
from pathlib import Path
from fractions import Fraction as Q
from .digits import encode_fraction,decode_fraction
from .lattice import integral_coset,quintic_response
from .boundary import boundary_system
from .rank_forecast import unknown_degrees
from .model import load_exact_dataset,InverseGenusRegressor,score
from .anomaly_learning import benchmark
from engine.bcov_yy import Poly

def main():
    root=Path(__file__).resolve().parents[1]
    refs=json.loads((root/'reference_amplitudes.json').read_text())['amplitudes']
    for v in [Q(0),Q(-1953125,36288),Q(5)**1000/Q(3)**250]:
        assert decode_fraction(encode_fraction(v))==v
    # Independent exact boundary checks. g4 includes corrected YY coefficient vector.
    for amp in refs[:3]:
        g=amp['genus'];a=boundary_system(g,Poly.from_data(amp['anomaly']),mum_degrees=[0])
        true=tuple(map(Q,amp['ambiguity']))
        assert a.verify(true)
        assert len(a.kernel)==(1 if g==4 else 0)
        if g<=3:assert a.particular==true
    g7=refs[5];part=boundary_system(7,Poly.from_data(g7['anomaly']),mum_degrees=[0])
    assert len(part.kernel)==2 and part.verify(tuple(map(Q,g7['ambiguity'])))
    audit=json.loads((root/'results/lattice_audit.json').read_text())
    for d in [4,5,6]:
        rec=next(x for x in audit['rows'] if x['g']==3 and x['degree']==d)
        coset=integral_coset(rec['intercept'],quintic_response(3,d),rec['cap'])
        z,n=coset.closest_index(0)
        assert z==Q(rec['chosen_coefficient']) and n==int(rec['chosen_GV'])
    assert unknown_degrees(53,((51,20),))==[]
    assert unknown_degrees(54)==[21]
    import numpy as np
    with np.load(root/'hybrid/results/float_dataset.npz') as cache:
        dataset={'g':cache['g'],'cum':cache['cum']}
    fit=InverseGenusRegressor().fit(dataset)
    holdout=score(fit,dataset,54,1000)
    assert holdout['model']['median_relative_error']<.001
    ant=benchmark(root/'reference_amplitudes.json')
    assert ant['model_sign_accuracy']>.75
    print('PASS: token roundtrip, rank, exact physics, lattice, low-X holdout, anti-holo benchmark')
    print('model median',holdout['model']['median_relative_error'],'anti sign',ant['model_sign_accuracy'])

if __name__=='__main__':main()
