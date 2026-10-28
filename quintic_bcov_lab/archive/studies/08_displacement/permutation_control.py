"""Exact permutation controls for descriptive top-family pattern discovery.

Not a statistical test of a new physics hypothesis: the extremal dataset is
algorithmically derived from a known theorem and the hypothesis class was
chosen after seeing the prior quintic work.
"""
from itertools import combinations, permutations
from pathlib import Path
import csv, json
R=Path(__file__).resolve().parent
train=[a for a in csv.DictReader((R/'data'/'extremal_displacements.csv').open()) if a['split']=='train']
sgn=[1 if int(a['n'])>0 else -1 for a in train]
mag=[abs(int(a['n'])) for a in train]
assert len(sgn)==8 and len(mag)==8
balanced_total=0;period4_hits=0
for plus_positions in combinations(range(8),4):
    labels=[1 if i in plus_positions else -1 for i in range(8)]
    balanced_total+=1
    if all(labels[i]==labels[i+4] for i in range(4)):period4_hits+=1
mag_total=0;poly2_hits=0
for tup in permutations(mag):
    mag_total+=1
    if all(tup[i+3]-3*tup[i+2]+3*tup[i+1]-tup[i]==0 for i in range(5)):
        poly2_hits+=1
out={'conditioned_sign_pattern':{'balanced_sign_sequences':balanced_total,
      'period4_repetitions':period4_hits,'fraction':period4_hits/balanced_total},
    'quadratic_pattern':{'permutations_of_observed_magnitudes':mag_total,
      'exact_quadratic_sequences':poly2_hits,'fraction':poly2_hits/mag_total},
    'interpretation':'Illustrative combinatorial null. Not a blind discovery significance assessment; labels come from an analytic theorem and hypothesis was prefigured.'}
(R/'results'/'permutation_control.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
