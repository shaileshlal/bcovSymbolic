#!/usr/bin/env python3
"""ISOLATED VALIDATOR: rank-zero DT reconstruction and wall crossing.

Never import this from predictor/*.py. Assumptions: the quintic rank-zero
DT modular ansatz and the stability/wall-crossing conditions of the cited
AFKPS work; these are *not* independent unconditional enumerative proofs.
"""
import csv,math
from pathlib import Path
from pt_from_dt_solver import pt_from_rankzero_dt,pt_mC,wall_candidates,applicable

ROOT=Path(__file__).resolve().parent.parent
DEPTH={0:6,1:5,2:4,3:4,4:5}
DELTA={0:0,1:1,2:2,3:2,4:1}

def reference(deg,t):
    m=pt_mC(deg)+DELTA[deg%5]+t
    if not applicable(deg,m):raise ValueError('Wall equation outside applicability',deg,t,m)
    value,base,terms=pt_from_rankzero_dt(deg,m)
    return {'degree':deg, 'shell':t, 'm':m, 'pt':str(value), 'vacuum_pt':str(base),
            'wall_candidate_count':len(wall_candidates(deg,m)),
            'nonzero_wall_count':sum(1 for r in terms if r['contribution']),
            'reference_type':'conditional_rankzero_DT_wallcrossing'}

def read(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def write(rows,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def extremal_multiples_five():
    # Liu--Ruan extremal formula, independent of the modular DT *ansatz*.
    vals=[]
    for d in (35,40,45,50):
        k=d//5
        x=math.comb(k+3,3)-math.comb(k-2,3)
        val=(-1)**(x+3)*5*x
        row=reference(d,0)
        assert int(row['pt'])==val,(d,val,row['pt'])
        vals.append({'degree':d,'leading_PT':str(val),'extremal_GV':str(val),'check':'pass'})
    return vals

def main():
    # Reproduce independently published PT reference (no training fitted here).
    archive=read(ROOT/'predictor/data/observed_Q6_Q30.csv')
    checks=0
    for row in archive:
        d=int(row['degree']);t=int(row['shell'])
        if d<21:continue
        got=reference(d,t)
        assert int(got['pt'])==int(row['pt']), (d,t,got['pt'],row['pt'])
        checks+=1
    assert checks==48,checks
    future=[reference(d,t) for d in range(31,54) for t in range(DEPTH[d%5])]
    assert len(future)==109
    assert all(r['nonzero_wall_count']==0 for r in future)
    out=ROOT/'results/reference_Q31_Q53.csv'
    write(future,out)
    extremal=extremal_multiples_five()
    write(extremal,ROOT/'results/extremal_geometry_Q35_Q50.csv')
    print(f'Validated reference on {checks} published PT values Q21..30; {len(future)} new sites Q31..53; 4 geometric checks')
    print('Wall candidates with nonzero potential:',[(r['degree'],r['shell'],r['wall_candidate_count']) for r in future if r['wall_candidate_count']])

if __name__=='__main__':main()
