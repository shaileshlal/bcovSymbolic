"""Exact denominator-lattice invariance tests for pure-X polynomial under X->1-X."""
from __future__ import annotations
from fractions import Fraction as Q
from math import comb, lcm
from pathlib import Path
import gzip,json,csv,argparse

def coeff_H(obj):
    g=obj['genus'];D=3*g-3;res=[Q(0)]*(D+1)
    for r,n,j,k,ell,num,den in obj['terms']:
        co=Q(int(num),int(den))*Q(-3,5)**j*Q(2,25)**k*Q(-4,125)**ell
        for a in range(r+1):
            t=n+ell+a
            if t<=D:res[t]+=co*((-1)**a)*comb(r,a)
    return res

def invert_X(f,g):
    # a_Z(X)=(-1)^(g-1) a_H(1-X)
    D=len(f)-1;return [(-1)**(g-1+s)*sum((f[k]*comb(k,s) for k in range(s,D+1)),Q(0)) for s in range(D+1)]

def lcd(v):
    out=1
    for x in v:out=lcm(out,x.denominator)
    return out

def other_than_five(x):
    while x%5==0 and x>0:x//=5
    return x

def main(source,output,gmax):
    rows=[]
    with gzip.open(source,'rt') as fh:
        for line in fh:
            obj=json.loads(line);g=obj['genus']
            if g>gmax:break
            h=coeff_H(obj);z=invert_X(h,g);yy=[x*5**(3*g-3) for x in h]
            dh,dz,dy=lcd(h),lcd(z),lcd(yy)
            assert dh==dz
            assert other_than_five(dh)==other_than_five(dy)
            rows.append({'genus':g,'degree':3*g-3,'lcd_HKQ':str(dh),'lcd_small_z':str(dz),'lcd_original_YY':str(dy),
                         'lcd_non5':str(other_than_five(dh)), 'denominator_lattice_check':True})
            print('g',g,'HKQ/z LCD equal?',dh==dz,'YY agrees away from 5?',other_than_five(dh)==other_than_five(dy),flush=True)
    with open(output,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    print('PASS',len(rows),'genus tables')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',default=str(Path(__file__).resolve().parent/'data'/'bcov_hkq_atomic_genus_2_15.jsonl.gz'));p.add_argument('--output',default=str(Path(__file__).resolve().parent/'denominator_invariants.csv'));p.add_argument('--genus',type=int,default=8);o=p.parse_args();main(o.source,o.output,o.genus)
