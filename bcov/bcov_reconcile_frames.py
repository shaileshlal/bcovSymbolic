"""Coordinate-corrected HKQ/YY conversion of uploaded BCOV atomic files.

Input is a gzip JSONL file written by bcov_sparse_converter.py, containing
P_g in the z-frame as sums of coeff * Xz^r*(Xz-1)^n * az^j*bz^k*cz^l.

Output HKQ representation uses
    alpha=v1-3/5, beta=v2+2/25, gamma=v3-X*(v2+4/125),
    P_g^HKQ=sum coeff_H * (1-X)^r * X^n * alpha^j * beta^k * gamma^l.

The ratio of coefficients is EXACTLY (-1)^(g-1+n+j+l).

Original YY uses the fifth-root coordinate; its connection generators
satisfy v_i^YY=5^i v_i^HKQ, X^YY=X^HKQ, and
P_g^YY = 5^(3*g-3)*P_g^HKQ(v_i^YY/5^i,X).
"""
from pathlib import Path
from fractions import Fraction
import argparse, gzip, json

def convert_file(path,output):
    stats=[]
    with gzip.open(path,'rt') as fi,gzip.open(output,'wt') as fo:
        for line in fi:
            obj=json.loads(line)
            g=obj['genus']
            rows=[]
            for r,n,j,k,l,num,den in obj['terms']:
                val=Fraction(int(num),int(den))
                val*=(-1)**(g-1+n+j+l)
                rows.append([r,n,j,k,l,str(val.numerator),str(val.denominator)])
            fo.write(json.dumps({'genus':g,'basis':['1-X','X','v1-3/5','v2+2/25','v3-X*(v2+4/125)'],'terms':rows},separators=(',',':'))+'\n')
            stats.append((g,len(rows)))
    return stats

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',nargs='?',default='bcov_yy_atomic_genus_2_15.jsonl.gz')
    parser.add_argument('--output',default='bcov_hkq_atomic_genus_2_15.jsonl.gz')
    args=parser.parse_args()
    print(convert_file(args.input,args.output))
