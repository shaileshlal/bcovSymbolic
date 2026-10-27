"""Exact, sparse parser and propagator-to-YY atomic converter for ff_2..ff_15.

Requires only the Python standard library. Does not expand (X-1)^n or shifted
YY variables a=v1+3/5, b=v2+2/25, c=v3+(1-X)(v2+4/125).

Input symbols: z, Cf[z], S2[z], tS1[z], tS0[z].
The *interpretation* Cf[z]=C_zzz and the quoted propagator choices MUST be
validated against independently normalized YY free energies. This script
is therefore a convention-dependent transformation, NOT a validated
holomorphic-ambiguity extraction.
"""
from __future__ import annotations
import ast, gzip, json, re, sys, glob, argparse
from fractions import Fraction
from collections import defaultdict, Counter
from pathlib import Path

sys.setrecursionlimit(30000)
ZER=(0,0,0,0,0)  # exponents of (z,Cf,S2,tS1,tS0)
NAMES=('z','Cf','S2','tS1','tS0')
MON={n: tuple(int(i==j) for i in range(5)) for j,n in enumerate(NAMES)}

def clean(p):
    return {m:v for m,v in p.items() if v}

def add(a,b,sign=1):
    out=a.copy()
    for m,c in b.items():
        out[m] = out.get(m,Fraction(0)) + sign*c
        if not out[m]: del out[m]
    return out

def mul(a,b):
    if len(a)>len(b): a,b=b,a
    if len(a)==1:
        (m,c),=a.items()
        return clean({tuple(x+y for x,y in zip(m,n)):c*d for n,d in b.items()})
    out=defaultdict(Fraction)
    for m,c in a.items():
        for n,d in b.items():
            out[tuple(x+y for x,y in zip(m,n))]+=c*d
    return clean(out)

def power(a,n):
    if not isinstance(n,int) or n<0: raise ValueError(f'bad power: {n}')
    if n==0: return {ZER:Fraction(1)}
    if len(a)==1:
        (m,c),=a.items()
        return {tuple(x*n for x in m):c**n}
    out={ZER:Fraction(1)}
    while n:
        if n&1: out=mul(out,a)
        n >>=1
        if n:a=mul(a,a)
    return out

def parse_ast(t):
    t=re.sub(r'\b(Cf|S2|tS1|tS0)\[z\]',r'\1',t).replace('^','**')
    tree=ast.parse(t,mode='eval')
    def ev(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int):
            return {ZER:Fraction(n.value)}
        if isinstance(n,ast.Name) and n.id in MON:
            return {MON[n.id]:Fraction(1)}
        if isinstance(n,ast.UnaryOp) and isinstance(n.op, ast.USub):
            return {m:-v for m,v in ev(n.operand).items()}
        if isinstance(n,ast.UnaryOp) and isinstance(n.op, ast.UAdd):return ev(n.operand)
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Pow):
                if not isinstance(n.right,ast.Constant) or not isinstance(n.right.value,int): raise ValueError('power not integer')
                return power(ev(n.left),n.right.value)
            if isinstance(n.op,ast.Add):return add(ev(n.left),ev(n.right))
            if isinstance(n.op,ast.Sub):return add(ev(n.left),ev(n.right),-1)
            if isinstance(n.op,ast.Mult):return mul(ev(n.left),ev(n.right))
            if isinstance(n.op,ast.Div):
                lhs,rhs=ev(n.left),ev(n.right)
                if len(rhs)!=1 or ZER not in rhs: raise ValueError('denominator not constant')
                c=rhs[ZER]
                return {m:v/c for m,v in lhs.items()}
        raise ValueError(f'unsupported expression: {ast.dump(n)[:150]}')
    return ev(tree.body)

def transform(term,g,kappa=5,mu=Fraction(1,3125)):
    # maps monomial z^p Cf^i S2^j tS1^k tS0^l to
    # coeff * X^r (X-1)^n * a^j b^k c^l, where
    # a=v1+3/5; b=v2+2/25; c=v3+(1-X)(v2+4/125)
    (p,i,j,k,l),coeff=term
    n=p-3*i+2*j+k
    r=g-1+i-j-k-l-n
    coeff=coeff * (-1)**j * Fraction(1,2)**l * Fraction(kappa)**(g-1+i-j-k-l) * (mu**n)
    return (r,n,j,k,l),coeff

def process_genus(path,g,kappa=5,mu=Fraction(1,3125)):
    poly=parse_ast(Path(path).read_text())
    converted=defaultdict(Fraction)
    for t in poly.items():
        e,v=transform(t,g,kappa=kappa,mu=mu)
        converted[e]+=v
    converted=clean(converted)
    z_res=Counter()
    minr=0
    for (r,n,j,k,l),c in converted.items():
        z_res[n]+=1
        minr=min(r,minr)
    return poly,converted, {'genus':g,'raw_monomials':len(poly),'atomic_terms':len(converted),'min_z_exponent':min(z_res),'max_z_exponent':max(z_res),'min_X_exponent':minr,'max_shift_degree':max(j+k+l for r,n,j,k,l in converted),'size_bytes':Path(path).stat().st_size}

def write_atomic(outfile,genus_data):
    # each row: [r,n,j,k,l,num,den]. This exactly encodes the atomic representation.
    with gzip.open(outfile,'wt') as f:
        for g,(_,poly,meta) in genus_data.items():
            rows=[[r,n,j,k,l,str(v.numerator),str(v.denominator)]
                  for (r,n,j,k,l),v in sorted(poly.items())]
            f.write(json.dumps({'genus':g,'terms':rows},separators=(',',':'))+'\n')

if __name__ == '__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir',default='.',help='directory containing ff_2.txt etc')
    ap.add_argument('--output',default='bcov_yy_atomic_genus_2_15.jsonl.gz')
    ap.add_argument('--manifest',default='bcov_yy_manifest.json')
    ap.add_argument('--kappa',type=int,default=5)
    ap.add_argument('--mu',type=Fraction,default=Fraction(1,3125))
    args=ap.parse_args()
    paths=sorted(glob.glob(str(Path(args.input_dir)/'ff_*.txt')),
        key=lambda p:int(re.search(r'ff_(\d+)',p).group(1)))
    if not paths: ap.error('no ff_N.txt inputs found')
    data={}
    for path in paths:
        g=int(re.search(r'ff_(\d+)',path).group(1))
        original,converted,meta=process_genus(path,g,args.kappa,args.mu)
        print(meta,flush=True)
        data[g]=(original,converted,meta)
    write_atomic(args.output,data)
    Path(args.manifest).write_text(json.dumps([meta for _,_,meta in data.values()],indent=2))
    print('Wrote',args.output,flush=True)
