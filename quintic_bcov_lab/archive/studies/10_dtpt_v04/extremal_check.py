"""Independent geometric checks of two out-of-table d=5k predictions."""
from math import comb
from pathlib import Path
import csv,json

def geometric_5k(k):
 n=comb(k+3,3)-comb(k-2,3)
 return (-1 if (n+3)%2 else 1)*5*n
p=Path(__file__).resolve().parent
pred={int(r['Q']):int(r['predicted_PT']) for r in csv.DictReader((p/'results/isolated_predictions_Q21_Q40.csv').open())}
checks=[]
for k in [7,8]:
 q=5*k
 geometric=geometric_5k(k)
 assert pred[q]==geometric
 checks.append({'Q':q,'genus':1+(5*k*(k+1))//2,'predicted':pred[q], 'geometric_formula':geometric,'independent_exact_match':True})
(p/'results/extremal_forecast_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(checks)
