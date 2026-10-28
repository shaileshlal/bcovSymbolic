"""Prepare physically separate files with targets removed from predictor inputs."""
from pathlib import Path
import csv
from empirical_benchmark import SOURCE
ROOT=Path(__file__).resolve().parent

def save(name,fieldnames,rows):
 p=ROOT/'data'/name
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fieldnames);w.writeheader();w.writerows(rows)

save('calibration_Q6_Q20.csv',['Q','PT_integer'],
     [{'Q':q,'PT_integer':v} for q,v in SOURCE if q<=20])
save('withheld_sites_Q21_Q40.csv',['Q'],[{'Q':q} for q in range(21,41)])
save('comparison_only_Q21_Q30.csv',['Q','PT_integer'],
     [{'Q':q,'PT_integer':v} for q,v in SOURCE if 21<=q<=30])
print('Wrote strict train/metadata/reference separation')
