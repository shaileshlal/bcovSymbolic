"""This is the ONLY script that reads the comparison-only Q>=21 labels."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PRED={int(r['Q']):int(r['predicted_PT']) for r in csv.DictReader((ROOT/'results/isolated_predictions_Q21_Q40.csv').open())}
REFERENCE={int(r['Q']):int(r['PT_integer']) for r in csv.DictReader((ROOT/'data/comparison_only_Q21_Q30.csv').open())}
assert len(REFERENCE)==10
exact={q:PRED[q]==v for q,v in REFERENCE.items()}
assert all(exact.values())
res={'withheld_exact_matches':sum(exact.values()),'withheld_total':len(exact),
     'Q21_prediction_without_oracle':PRED[21],
     'Q31_to_Q40_unchecked':{str(q):PRED[q] for q in range(31,41)}}
(ROOT/'results/isolated_validation.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
