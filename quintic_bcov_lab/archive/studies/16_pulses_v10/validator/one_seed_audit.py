#!/usr/bin/env python3
"""Post-validation sensitivity: one withheld integer fixes four additional sites.
Not part of the no-oracle benchmark. Reference is explicitly consulted here.
"""
import sys,json
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'predictor'))
from recursive_defects import read,phase_fit
from deep_envelope import W
from pulse_profile import K

def main():
    ref={(int(z['degree']),int(z['shell'])):int(z['pt']) for z in read(ROOT/'validator/reference_Q20_Q53.csv') if z['pt']}
    period,phase=phase_fit(read(ROOT/'data/shallow_calibration_Q1_Q20.csv'))
    params=json.loads((ROOT/'results/fitted_pulses.json').read_text())
    profiles={(z['residue'],z['front_distance']):(F(z['P']),F(z['Q'])) for z in params['profiles']}
    # Seed: true reference PT(43,10) (NOT a prediction).
    d,s,t,k=43,3,10,8
    P,Q=profiles[s,K(s,t)-k]
    scale=F(ref[d,t],phase[s,k%period]*W(s,t,k))
    B=scale-(-1)**t*(P+Q/F(W(s,t,k)))
    checked=[]
    for d in (47,48,52,53):
        s=d%5;t=10;k=d//5;r=K(s,t)-k
        val=B if r<=0 else B+(-1)**t*(profiles[s,r][0]+profiles[s,r][1]/F(W(s,t,k)))
        out=val*phase[s,k%period]*W(s,t,k)
        assert out.denominator==1
        assert out.numerator==ref[d,t]
        checked.append({'degree':d,'shell':t,'prediction':str(out.numerator),'exact':True})
    record={'seed_degree':43,'seed_shell':10,'seed_PT':str(ref[43,10]),
            'inferred_B_23_shell10':str(B),'new_predicted_exact':4,
            'new_predictions':checked,'remaining_unknown_without_second_seed':{'degree':42,'shell':10}}
    (ROOT/'validator/one_seed_sensitivity.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
