"""Exact residual quintic YY ambiguity: no floats in coefficients.

Coordinate x = psi**(-5) = 3125*z_PF; X=x/(x-1).
This module constructs *reference-derived* canonical affine coordinates,
NOT an independently boundary-solved affine representative at all genera.
A homogeneous kernel factorization is proved by orbifold/conifold/MUM
boundary analysis and checked by direct boundary linear algebra for g<=9.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json, csv, math


def dimensions(g: int) -> tuple[int,int,int]:
    if g < 2: raise ValueError('g >= 2')
    K=3*g-3
    m=(K+4)//5
    D=(2*g-2)//5
    assert m+D==g-1
    return K,m,D


def kernel_vector(g: int, r: int) -> list[F]:
    """X**(r+1)*(1-X)**m as a coefficient vector."""
    K,m,D=dimensions(g)
    if r < 0 or r >= D: raise ValueError('residual index outside kernel')
    v=[F(0) for _ in range(K+1)]
    for k in range(m+1): v[1+r+k]=F((-1)**k*comb(m,k))
    return v


def add_vectors(vectors, weights, length):
    z=[F(0)]*length
    for v,w in zip(vectors,weights):
        for j,x in enumerate(v): z[j]+=w*x
    return z


def reference_split(g: int, physical: list[F]) -> dict:
    """Canonical boundary gauge: anchor middle slots m+1,...,g-1 are zero.

    Each c_r is computed from the physical reference. This extraction DOES NOT
    predict c_r. For an independent prediction, replace the reference-derived
    anchor by an independently computed exact boundary particular solution.
    """
    K,m,D=dimensions(g)
    if len(physical)!=K+1: raise ValueError('ambiguity length mismatch')
    basis=[kernel_vector(g,r) for r in range(D)]
    c=[F(0)]*D
    # Upper triangular: the coefficient of X**(m+1+r) in basis r is (-1)**m.
    for r in range(D-1,-1,-1):
        idx=m+1+r
        higher=sum((c[t]*basis[t][idx] for t in range(r+1,D)), F(0))
        c[r]=F((-1)**m)*(physical[idx]-higher)
    delta=add_vectors(basis,c,K+1)
    anchor=[p-d for p,d in zip(physical,delta)]
    assert all(anchor[i]==0 for i in range(m+1,g))
    assert all(anchor[i]+delta[i]==physical[i] for i in range(K+1))
    assert delta[0]==0 and all(delta[i]==0 for i in range(g,K+1))
    # Reduced section: f=(x-1)**(g-1)*a(x/(x-1)).
    # delta f= sum b_r x**(r+1)*(1-x)**(D-r-1), r=0..D-1.
    b=[(-1)**(g-2-r)*c[r] for r in range(D)]
    xQ=[F(0)]*(D+1)
    for r,br in enumerate(b):
        for j in range(D-r):
            xQ[1+r+j] += br * F((-1)**j * comb(D-1-r,j))
    assert xQ[0]==0
    # Alternative transformation: difference section from X polynomial.
    sec=[F(0)]*(g)
    for j,dj in enumerate(delta[:g]):
        for k in range(g-j):
            sec[j+k]+=dj*F((-1)**(g-1-j-k)*comb(g-1-j,k))
    assert sec[:D+1]==xQ and all(s==0 for s in sec[D+1:])
    return {'genus':g, 'K':K, 'm':m,'D':D, 'pivot_indices':list(range(m+1,g)),
            'physical':physical,'anchor':anchor,'residual':delta,
            'c':c,'bernstein_x':b,'xQ':xQ,'basis':basis}


def _signed_log10_abs(x: F) -> float:
    if not x:return 0.
    def log_big(a: int) -> float:
        s=str(abs(a));return len(s)-1+math.log10(int(s[:15])/10**(len(s[:15])-1))
    return log_big(x.numerator)-log_big(x.denominator)


def mixed_features(anomaly, xdegree: int) -> dict:
    """Summaries of the exact nonholomorphic source at the selected X slice."""
    matches=[]
    for exps,n,d in anomaly:
        if exps[4]==xdegree and any(exps[1:4]):
            matches.append((exps,F(int(n),int(d))))
    mag=[abs(v) for _,v in matches]
    return {'mixed_terms_at_Xj':len(matches),
            'mixed_v1_terms_at_Xj':sum(e[1]>0 for e,_ in matches),
            'mixed_v2_terms_at_Xj':sum(e[2]>0 for e,_ in matches),
            'mixed_v3_terms_at_Xj':sum(e[3]>0 for e,_ in matches),
            'mixed_log10_max_at_Xj':max((_signed_log10_abs(t) for t in mag),default=None)}


def load_reference(path: Path):
    data=json.loads(path.read_text())
    if data.get('convention') is None: raise ValueError('Missing YY gauge convention')
    for obj in data['amplitudes']:
        g=int(obj['genus'])
        yield g,[F(x) for x in obj['ambiguity']],obj['anomaly']


def serialize(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,tuple):return [serialize(x) for x in v]
    if isinstance(v,list):return [serialize(x) for x in v]
    if isinstance(v,dict):return {k:serialize(x) for k,x in v.items() if k!='basis'}
    return v


def build(reference: Path, outdir: Path):
    outdir.mkdir(parents=True,exist_ok=True)
    dataset=[];samples=[];validation=[]
    for g,physical,anomaly in load_reference(reference):
        item=reference_split(g,physical)
        dataset.append(serialize(item))
        for r in range(item['D']):
            j=item['m']+1+r
            row={'genus':g,'residual_index':r,'X_pivot_degree':j,
                 'K':item['K'],'m':item['m'],'D':item['D'],
                 'relative_pivot_degree':j/item['K'],
                 'a_pivot_phys':str(physical[j]),
                 'c_X_kernel':str(item['c'][r]),
                 'b_x_bernstein':str(item['bernstein_x'][r]),
                 'label_sign':0 if not item['c'][r] else (1 if item['c'][r]>0 else -1),
                 'label_signed_log10':(0.0 if not item['c'][r] else (1 if item['c'][r]>0 else -1)*(_signed_log10_abs(item['c'][r]) if _signed_log10_abs(item['c'][r])>18 else math.log10(1+10**_signed_log10_abs(item['c'][r])))),
                 **mixed_features(anomaly,j)}
            samples.append(row)
        validation.append({'genus':g,'original_coefficients':3*g-2,
                           'residual_dimension':item['D'],
                           'round_trip_exact':True,'xQ_degree':max([i for i,v in enumerate(item['xQ']) if v],default=-1)})
    (outdir/'intermediate_exact.json').write_text(json.dumps({'schema':'quintic-yy-residual-v1',
        'gauge':'X=x/(x-1), x=psi^-5; f=(x-1)^(g-1)*a(X)',
        'provenance':'corrected reference amplitudes g=2..12; anchors REFERENCE-DERIVED, not independently computed',
        'coordinates':'c: X(1-X)^m R(X), b: x^(r+1)(1-x)^(D-r-1)',
        'generations':dataset},indent=2)+'\n')
    if samples:
        with (outdir/'intermediate_training.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(samples[0].keys()))
            writer.writeheader();writer.writerows(samples)
    (outdir/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    return dataset,samples,validation


def bcov02_manifest(root: Path) -> dict:
    """Non-invasive bridge for the separate v0.2 structural-learner corpus.

    V0.2 tasks are optional model features. We never use a task as a
    numerical boundary/ambiguity label without a verified matching schema.
    """
    root=Path(root)
    result={'present':root.exists(),'root':str(root),'tasks':0,'registries':0,'examples':[]}
    if not root.exists():return result
    for p in root.rglob('tasks.jsonl'):
        result['tasks']+=sum(1 for line in p.open() if line.strip())
        result['examples'].append(str(p.relative_to(root)))
    for p in root.rglob('structure_registry.json'):
        result['registries']+=1
        result['examples'].append(str(p.relative_to(root)))
    return result


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',default='data/reference_amplitudes.json')
    ap.add_argument('--output',default='data')
    a=ap.parse_args()
    data,samples,aud=build(Path(a.source),Path(a.output))
    print(f'Exact split: {len(data)} genera; {len(samples)} residual targets; {len(aud)} round-trip checks')
