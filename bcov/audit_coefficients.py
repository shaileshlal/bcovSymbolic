"""Exact low-X audit of Klemm-supplied ff_2..ff_15 in HKQ and original YY frames.
Only expands at v=0 and requested X degree; no full polynomial expansion.
"""
from __future__ import annotations
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import gzip,json,csv,argparse

def bernoulli_abs(N):
    b={1:Q(1,6)}
    for g in range(2,N+1):
        b[g]=sum((Q(comb(2*g,2*h))*b[h]*b[g-h] for h in range(1,g)),Q(0))/Q(2*g+1)
    return b

def constant(g,b):
    from math import factorial
    return Q(100)*Q(625)**(g-1)*b[g]*b[g-1]/Q(2*g*(2*g-2)*factorial(2*g-2))

def linear(g,b,c):
    from math import factorial
    if g==2:return Q(-28795,144)
    S=sum((Q((2*h-2)*(2*g-2*h-2))*c[h]*c[g-h] for h in range(2,g-1)),Q(0))
    I=Q(23,25)*Q(-625)**(g-1)*b[g]/Q(2*g*factorial(2*g-2))
    return -Q(577,625)*(g-1)*c[g]+Q(2*g-4)*(62-Q(24*g,5))*c[g-1]-Q(12,5)*S-I

def pure_x_hkq(path, degree=4):
    # P^H_g = sum coeff*(1-X)^r X^n*(v1-3/5)^j*(v2+2/25)^k*(v3-X*(v2+4/125))^ell.
    out={}
    with gzip.open(path,'rt') as f:
        for line in f:
            obj=json.loads(line);g=obj['genus'];co=[Q(0)]*(degree+1)
            for r,n,j,k,ell,num,den in obj['terms']:
                v=n+ell
                if v>degree:continue
                w=Q(int(num),int(den))*Q(-3,5)**j*Q(2,25)**k*Q(-4,125)**ell
                for s in range(min(r,degree-v)+1):
                    co[v+s]+=w*((-1)**s)*comb(r,s)
            out[g]=co
    return out

def audit(source,out):
    p=pure_x_hkq(source);N=max(p);b=bernoulli_abs(N);c={g:constant(g,b) for g in range(2,N+1)}
    rows=[]
    for g,coeff_h in p.items():
        k=Q(5)**(3*g-3)
        coeff_y=[k*x for x in coeff_h]
        c0=c[g]
        c1=linear(g,b,c)
        rows.append(dict(genus=g, a0_exact=str(coeff_y[0]),a0_expected=str(c0),a0_ok=coeff_y[0]==c0,
                         a1_exact=str(coeff_y[1]),a1_expected=str(c1),a1_ok=coeff_y[1]==c1,
                         r1=str(coeff_y[1]/coeff_y[0]),r2=str(coeff_y[2]/coeff_y[0]),
                         r3=str(coeff_y[3]/coeff_y[0]),r4=str(coeff_y[4]/coeff_y[0])))
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    for r in rows:print('g=',r['genus'],'a0=',r['a0_ok'],'a1=',r['a1_ok'],'r1=',r['r1'][:92])
    assert all(r['a0_ok'] and r['a1_ok'] for r in rows), 'recurrence disagreement! inspect audit CSV'
    print('PASS: all constant and linear Bernoulli identities through genus',N)
    return rows

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source',default=str(Path(__file__).resolve().parent/'data'/'bcov_hkq_atomic_genus_2_15.jsonl.gz'));ap.add_argument('--output',default=str(Path(__file__).resolve().parent/'low_x_reaudit_genus_2_15.csv'));o=ap.parse_args();audit(o.source,o.output)
