#!/usr/bin/env python3
"""Conditional shell-composability experiment (NOT an independent holdout).

Uses first two actual NEW shell pulses to fit a linear factor in the shell
variable y and checks deeper shells. This deliberately uses target-block labels;
never count it as an unconditional forecasting result.
"""
from fractions import Fraction
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run():
    data=list(csv.DictReader((ROOT/'results/defect_shells.csv').open()))
    D={}
    for r in data:
        s,t=int(r['residue']),int(r['shell'])
        D[s,t]=(int(r['defect_at_first_update']),int(r['defect_at_next_update']))
    rows=[]
    for s in range(5):
        older={t:v[0] for (ss,t),v in D.items() if ss==s}
        newer={t:v[1] for (ss,t),v in D.items() if ss==s}
        old_support=[t for t,v in older.items() if v]
        new_support=[t for t,v in newer.items() if v]
        if not old_support or not new_support:continue
        h=min(old_support);hn=min(new_support);shift=hn-h
        a=Fraction(newer[hn],older[h]);b=None
        if hn+1 in newer:
            b=Fraction(newer[hn+1]-a*older.get(h+1,0),older[h])
        for t in sorted(newer):
            if b is None or t<hn+2:continue
            guess=a*older.get(t-shift,0)+b*older.get(t-shift-1,0)
            truth=newer[t]
            if not truth:continue
            rows.append(dict(residue=s,shell=t,shell_shift=shift,linear_coefficient=str(a),extra_y_coefficient=str(b),
                             actual_defect=str(truth),predicted_defect=str(guess),exact=(guess==truth),
                             relative_error=float(abs(guess-truth)/abs(truth))))
    assert len(rows)==2 and all(r['residue']==1 for r in rows)
    with (ROOT/'results/conditional_shell_transfer.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
    print('Conditional shell-factor transfer, not blind prediction:',rows)
if __name__=='__main__':run()
