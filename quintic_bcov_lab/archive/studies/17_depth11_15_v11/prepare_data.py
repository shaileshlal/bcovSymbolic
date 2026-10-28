"""Prepare exact integer calibration labels separately from the reference directory.
This utility is for release construction only. Predictor uses data/*.csv alone.
"""
from pathlib import Path
import csv
root=Path(__file__).parent
refs=[]
for r in csv.DictReader(open(root/'validator/reference_Q20_Q53.csv')):
 if r['status']=='certified_within_conditional_framework' and 20<=int(r['degree'])<=40:
  refs.append({'degree':r['degree'],'shell':r['shell'],'pt':r['pt']})
for r in csv.DictReader(open(root/'validator/reference_deep_11_15.csv')):
 if r['status']=='certified_conditional' and 20<=int(r['degree'])<=40:
  refs.append({'degree':r['degree'],'shell':r['shell'],'pt':r['pt']})
refs.sort(key=lambda r:(int(r['degree']),int(r['shell'])))
out=root/'data/observed_calibration_Q20_Q40_shell6_15.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=('degree','shell','pt'));w.writeheader();w.writerows(refs)
out=root/'data/withheld_sites_Q41_Q53_shell11_15.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=('degree','shell'));w.writeheader();w.writerows(({'degree':d,'shell':t} for d in range(41,54) for t in range(11,16)))
print('train',len(refs), 'train new shell',sum(int(r['shell'])>=11 for r in refs))
