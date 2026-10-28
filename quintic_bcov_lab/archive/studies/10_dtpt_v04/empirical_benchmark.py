"""Pre-frontier PT-calibrated arithmetic predictor, strict degree-ordered holdouts.

One PT leading value in each of the five degree mod-5 classes (Q=6,..,10)
is used to infer constant rank-zero DT amplitudes Ω_mu. Subsequent labels
through Q=30 are reference-only. Exact wall eligibility is reported.
"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import csv, json
from fractions import Fraction as F
from pt_from_dt_solver import pt_mC,parity_sign,applicable,wall_candidates,omega_rankzero,qhat,simple_vacuum_pt

HERE=Path(__file__).parent
BASE=HERE/'data'/'residue_family_benchmark.csv'
SOURCE=[]
with BASE.open() as f:
 for row in csv.DictReader(f):
  SOURCE.append((int(row['degree']),int(row['reference_integer'])))
# Extremal degree 5k inputs. Exclude Q=5 which is not in the single-wall regime.
SOURCE += [(10,-50),(15,-100),(20,175),(25,275),(30,-400)]
SOURCE=sorted(set(SOURCE))
assert len(SOURCE)==25, len(SOURCE)
LABELS=dict(SOURCE)

DELTA={0:0,1:1,2:2,3:2,4:1}

def site(Q):
    s=Q%5;delta=DELTA[s]; m=pt_mC(Q)+delta
    chi=m+Q-5; sign=parity_sign(chi)*chi
    return {'degree':Q,'residue':s,'genus':1-m,'m':m,'delta':delta,
            'chi':chi,'kinematic_multiplier':sign,
            'rankzero_mu':min(s,5-s),'rankzero_qhat':str(qhat(Q,m)),
            'wall_inequality_pass':applicable(Q,m),
            'nonvacuum_candidate_count':len(wall_candidates(Q,m))}

# Fit only Q<=10: one prototype per residue, including the Q=10 extremal seed.
TRAIN={Q:y for Q,y in SOURCE if Q<=10}
assert set(TRAIN)==set(range(6,11))
learned={}
for Q,y in TRAIN.items():
    st=site(Q)
    assert y%st['kinematic_multiplier']==0
    omega=y//st['kinematic_multiplier']
    mu=st['rankzero_mu']
    if mu in learned: assert learned[mu]==omega
    else: learned[mu]=omega
assert learned=={1:8625,2:-1218500,0:5}

ROWS=[]
for Q,ref in SOURCE:
    st=site(Q)
    pred=st['kinematic_multiplier']*learned[st['rankzero_mu']]
    assert pred==ref,(Q,pred,ref)
    assert pred==simple_vacuum_pt(Q,st['m'])
    row={**st,'split':('calibration' if Q<=10 else ('forward_pre_54' if Q<=20 else 'frontier_holdout')),
         'reference_PT_leading':ref,'prediction_PT':pred,'error':pred-ref,
         'independent_modular_prediction':simple_vacuum_pt(Q,st['m']),
         'nonvacuum_summands_possible':len(wall_candidates(Q,st['m']))>0}
    ROWS.append(row)

# Out-of-table predicted degrees. A proper theorem certificate requires both
# applicability and no possible nonvacuum terms, not just a matching fit.
FORECAST=[]
for Q in range(31,41):
    st=site(Q)
    pred=st['kinematic_multiplier']*learned[st['rankzero_mu']]
    FORECAST.append({**st,'predicted_PT_leading':pred,
        'theorem_single_wall_certified':st['wall_inequality_pass'] and not st['nonvacuum_candidate_count']})

(HERE/'results').mkdir(exist_ok=True)
with (HERE/'results'/'leading_PT_benchmark.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(ROWS[0]));w.writeheader();w.writerows(ROWS)
with (HERE/'results'/'forward_Q31_Q40.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(FORECAST[0]));w.writeheader();w.writerows(FORECAST)

splitcounts={split:{'n':sum(r['split']==split for r in ROWS),
                   'exact':sum(r['split']==split and r['prediction_PT']==r['reference_PT_leading'] for r in ROWS)}
            for split in set(r['split'] for r in ROWS)}
result={'method':'Single-calibration-per-residue, learned polar rank-zero DT coefficient * chi*parity',
        'Q6_to_Q10_calibration':TRAIN,'learned_rankzero_coefficients':learned,
        'scores':splitcounts,
        'range_without_any_leading_PT_labels':'11..30',
        'Q21_prediction':next(r for r in ROWS if r['degree']==21),
        'Q31_to_Q40_unchecked_predictions':FORECAST}
(HERE/'results'/'benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
print('Learned rank0 DT coefficients from five early PT values:',learned)
print('Heldout:',splitcounts)
for r in ROWS:
    if r['degree']>20:print('Q',r['degree'],'g',r['genus'],'expected',r['reference_PT_leading'],'pred',r['prediction_PT'],
       'wall-applicable',r['wall_inequality_pass'],'nonvacuum candidates',r['nonvacuum_candidate_count'])
print('Forecast:')
for r in FORECAST:
    print(r['degree'],r['genus'],r['predicted_PT_leading'],r['theorem_single_wall_certified'])
