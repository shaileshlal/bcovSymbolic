"""Exact scan of pure X coefficients and all monomial denominator valuations.

Standard library only. Input is a BCOV schema-1 checkpoint. No recurrence is
inferred by floating point fitting, and no denominator conjecture is treated
as proved. Explicit top-edge formulas are verified against each available genus.

Usage: python scan_coefficients.py physical/amplitudes.json --output results
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial, lcm
from pathlib import Path
import argparse
import json


@lru_cache(maxsize=None)
def bernoulli(n):
    if n==0: return Q(1)
    return -sum((comb(n+1,k)*bernoulli(k) for k in range(n)),Q(0))/(n+1)


@lru_cache(maxsize=None)
def bernoulli_magnitude_from_lower(g):
    if g==1: return Q(1,6)
    if g<1: raise ValueError('Positive genus index required')
    return sum((comb(2*g,2*h)*bernoulli_magnitude_from_lower(h)*
                bernoulli_magnitude_from_lower(g-h) for h in range(1,g)),Q(0))/(2*g+1)


def leading(g):
    return Q(5**(3*g-3))*abs(bernoulli(2*g))/(2*g*(2*g-2))


def vertex_edge(h):
    # Includes the genus-one one-point seed: k_1=5/12.
    return Q(5**(3*h-2))*abs(bernoulli(2*h))/(2*h)


def mixed_edge(g):
    """[v1 X^(3g-4)] P_g, from the top part of the BCOV source."""
    return -Q(1,2)*(5*(2*g-3)*vertex_edge(g-1)+sum(
        (vertex_edge(h)*vertex_edge(g-h) for h in range(1,g)),Q(0)))


def next_to_leading(g):
    return -Q(8,5)*(g-1)*leading(g)-2*mixed_edge(g)


def constant(g):
    # Equivalent to (-625)^(g-1) times the MUM constant map, chi=-200.
    return Q(100*625**(g-1))*abs(bernoulli(2*g)*bernoulli(2*g-2))/(
        2*g*(2*g-2)*factorial(2*g-2))


def primes(n):
    composite=set()
    out=[]
    for p in range(2,n+1):
        if p not in composite:
            out.append(p)
            composite.update(range(p*p,n+1,p))
    return out


def factors(n):
    """Exact trial factorization, used only for small-prime denominator products."""
    if n<1: raise ValueError('Positive integer required')
    out={}
    p=2
    while p*p<=n:
        while n%p==0:
            out[p]=out.get(p,0)+1
            n//=p
        p+=1
    if n>1: out[n]=out.get(n,0)+1
    return out


def small_numerator_factors(n,bound=97):
    n=abs(n)
    if not n:return {'small_prime_factors':{},'unfactored_cofactor':'0'}
    out={}
    for p in primes(bound):
        while n%p==0:
            out[p]=out.get(p,0)+1
            n//=p
    return {'small_prime_factors':out,'unfactored_cofactor':str(n)}


def val(n,p):
    if not n:return None
    n=abs(n);k=0
    while n%p==0:k+=1;n//=p
    return k


def rref_nullspace(rows,width):
    """Exact homogeneous linear algebra for exploratory rational-ratio fits."""
    a=[[Q(c) for c in row] for row in rows]
    pivots=[];r=0
    for col in range(width):
        pivot=next((j for j in range(r,len(a)) if a[j][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        lead=a[r][col];a[r]=[c/lead for c in a[r]]
        for j in range(len(a)):
            if j!=r and a[j][col]:
                factor=a[j][col]
                a[j]=[c-factor*d for c,d in zip(a[j],a[r])]
        pivots.append(col);r+=1
        if r==len(a):break
    basis=[]
    for free in (j for j in range(width) if j not in pivots):
        v=[Q(0)]*width;v[free]=Q(1)
        for i,col in enumerate(pivots):v[col]=-a[i][free]
        basis.append(v)
    return basis


def peval(coeffs,x):
    return sum((c*Q(x)**j for j,c in enumerate(coeffs)),Q(0))


def fit_ratios(sequence,train_transitions=5):
    """Test a_(k+1)/a_k=P(k)/Q(k), degrees <=2, on held-out transitions.

    A fit requires a one-dimensional nullspace and nonzero denominator at all
    tested points. The maximal 2/2 fit has five free parameters and exactly
    five training ratios; its interpolation alone is not evidence.
    """
    points=[(i,b/a) for (i,a),(j,b) in zip(sequence,sequence[1:]) if j==i+1 and a]
    train=points[:train_transitions];test=points[train_transitions:]
    if len(train)<train_transitions or not test:
        return {'status':'insufficient_data'}
    fits=[]
    for m in range(3):
        for n in range(3):
            rows=[[Q(x)**j for j in range(m+1)]+[-ratio*Q(x)**j for j in range(n+1)] for x,ratio in train]
            null=rref_nullspace(rows,m+n+2)
            if len(null)!=1:continue
            p,q=null[0][:m+1],null[0][m+1:]
            if not any(q) or any(peval(q,x)==0 for x,_ in points):continue
            failures=[{'index':x,'predicted':str(peval(p,x)/peval(q,x)),'actual':str(y)}
                      for x,y in test if peval(p,x)/peval(q,x)!=y]
            fits.append({'numerator_degree':m,'denominator_degree':n,
                         'p':list(map(str,p)),'q':list(map(str,q)),
                         'heldout_count':len(test),'heldout_failures':failures,
                         'passes_heldout':not failures})
    return {'status':'scanned','train_indices':[x for x,_ in train],
            'heldout_indices':[x for x,_ in test],'fitted_models':fits,
            'any_passes_heldout':any(x['passes_heldout'] for x in fits)}


def scan(path,output):
    data=json.loads(Path(path).read_text())
    if data['schema']!=1 or not data['convention'].startswith('YY2004-quintic-psi-'):
        raise ValueError('Unsupported coefficient convention')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    arrays={};stats: list[dict]=[];all_coefficients=[];checks=[]
    for amplitude in data['amplitudes']:
        g=amplitude['genus'];d=3*g-3
        if amplitude['ambiguity'] is None:raise ValueError(f'Unresolved genus {g}')
        c=[Q(x) for x in amplitude['ambiguity']];c+= [Q(0)]*(d+1-len(c))
        arrays[g]=c
        all_terms={tuple(m):Q(int(a),int(b)) for m,a,b in amplitude['anomaly']}
        all_terms.update({(0,0,0,0,n):v for n,v in enumerate(c) if v})
        den_pure=lcm(*(v.denominator for v in c))
        den_all=lcm(*(v.denominator for v in all_terms.values()))
        fp=factors(den_pure);fa=factors(den_all)
        expected={p:2*g//(p-1) for p in primes(2*g+1) if p>=7}
        checks.append({'genus':g,'leading':leading(g)==c[d],
                       'next_to_leading':next_to_leading(g)==c[d-1],
                       'constant':constant(g)==c[0],
                       'mixed_v1_edge':mixed_edge(g)==all_terms.get((0,1,0,0,d-1),0),
                       'bernoulli_denominator_sieve':factors(bernoulli(2*g).denominator)==
                          {p:1 for p in primes(2*g+1) if 2*g%(p-1)==0}})
        stats.append({'genus':g,'terms':len(all_terms),'pure_coefficients':len(c),
                      'pure_denominator_lcm':str(den_pure),'all_denominator_lcm':str(den_all),
                      'pure_denominator_factors':fp,'all_denominator_factors':fa,
                      'empirical_sieve_p_ge_7':expected,
                      'pure_sieve_matches':{p:v for p,v in fp.items() if p>=7}==expected,
                      'all_sieve_matches':{p:v for p,v in fa.items() if p>=7}==expected,
                      'all_5_integral':fa.get(5,0)==0,
                      'minimum_5_adic_valuation':min(val(v.numerator,5)-val(v.denominator,5) for v in all_terms.values()),
                      'constant_map_denominator_factors':factors(c[0].denominator),
                      'all_lcm_over_constant_map_denominator':str(Q(den_all,c[0].denominator)),
                      'constant_map_already_saturates_p_ge_7':{p:e for p,e in factors(c[0].denominator).items() if p>=7}==expected,
                      'all_no_prime_above_2g_plus_1':max(fa,default=1)<=2*g+1})
        for n,v in enumerate(c):
            all_coefficients.append({'genus':g,'n':n,'top_offset':d-n,'coefficient':str(v),
                'denominator_factors':factors(v.denominator),
                'v5':None if not v else val(v.numerator,5)-val(v.denominator,5),
                'normalized_by_5_to_d':str(v/5**d),**small_numerator_factors(v.numerator)})
    if not all(all(v for k,v in row.items() if k!='genus') for row in checks):
        raise AssertionError('An exact edge identity failed')
    models={}
    for offset in range(4):
        seq=[(g,c[3*g-3-offset]/5**(3*g-3)) for g,c in sorted(arrays.items()) if 3*g-3>=offset]
        models[f'genus_direction_top_offset_{offset}']=fit_ratios(seq)
    for g,c in arrays.items():
        if len(c)>6:models[f'within_genus_{g}']=fit_ratios(list(enumerate(c)))
    predictions=[{'genus':g,'n_top':3*g-3,'leading':str(leading(g)),
                  'next_to_leading':str(next_to_leading(g)),
                  'constant':str(constant(g)),
                  'mixed_v1_edge':str(mixed_edge(g)),
                  'full_amplitude_computed':g in arrays} for g in range(2,54)]
    assert all(bernoulli_magnitude_from_lower(g)==abs(bernoulli(2*g)) for g in range(1,54))
    summary={'bernoulli_binomial_convolution_verified_through':53,
             'convention':data['convention'],'computed_genera':sorted(arrays),
             'pure_coefficients_scanned':len(all_coefficients),
             'all_monomials_scanned':sum(s['terms'] for s in stats),
             'exact_identity_checks':checks,'denominator_statistics':stats,
             'denominator_pattern_status':'Empirical only for full coefficient arrays; no all-genus proof supplied.',
             'fresh_holdout_genera':[g for g in arrays if g>=9],
             'low_complexity_ratio_models':models}
    def save(name,obj):(output/name).write_text(json.dumps(obj,indent=2)+'\n')
    save('summary.json',summary)
    save('pure_x_coefficients.json',all_coefficients)
    save('exact_edge_predictions_through_53.json',predictions)
    print(json.dumps({'genera':sorted(arrays),'pure_coefficients':len(all_coefficients),
        'all_monomials':summary['all_monomials_scanned'],
        'edge_checks_passed':True,
        'full_denominator_sieve_matches':all(s['all_sieve_matches'] for s in stats),
        'any_validated_small_ratio_model':any(v.get('any_passes_heldout',False) for v in models.values())},indent=2))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint',nargs='?',default='physical/amplitudes.json')
    p.add_argument('--output',default='results')
    args=p.parse_args();scan(args.checkpoint,args.output)
