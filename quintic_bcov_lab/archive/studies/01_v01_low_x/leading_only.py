"""Exact scalar computations. No Airy descendant or dressed equations.

Standard Python only. Tangent-number integer arithmetic computes Bernoulli
numbers independently of the archived Bernoulli table.
"""
from fractions import Fraction as Q
from math import factorial, gcd, log10, lcm
from pathlib import Path
import argparse, gzip, hashlib, json, sys, time

def bernoulli_abs(limit):
    # Tangent-number triangular recurrence, T_(2n-1).
    t=[0]*(limit+1)
    t[1]=1
    for k in range(2,limit+1):t[k]=(k-1)*t[k-1]
    for k in range(2,limit+1):
        for j in range(k,limit+1):t[j]=(j-k)*t[j-1]+(j-k+2)*t[j]
    return [Q(0)]+[Q(2*g*t[g],2**(2*g)*(2**(2*g)-1)) for g in range(1,limit+1)]

def constant_coefficients(bb):
    return [Q(0),Q(0)]+[Q(100*625**(g-1),2*g*(2*g-2)*factorial(2*g-2))*bb[g]*bb[g-1] for g in range(2,len(bb))]

def linear_coefficient(g,bb,c):
    if g==2:return Q(-28795,144)
    terms=[]
    for h in range(2,g//2+1):
        terms.append(((1 if 2*h==g else 2)*(2*h-2)*(2*g-2*h-2)*c[h].numerator*c[g-h].numerator,
                      c[h].denominator*c[g-h].denominator))
    den=lcm(*(d for _,d in terms));s=Q(sum(n*(den//d) for n,d in terms),den)
    return (-Q(577,625)*(g-1)*c[g]+(2*g-4)*(62-Q(24*g,5))*c[g-1]
            -Q(12,5)*s-Q(23,25)*(-625)**(g-1)*bb[g]/(2*g*factorial(2*g-2)))

def leading_airy(limit):
    # U(s)=sum (-1)^n (6n)!/[(3n)!(2n)!] s^n; log U /288^n.
    u=[1];ell=[0];alpha=[Q(0),Q(0)]
    for n in range(1,limit):
        numerator=u[-1]
        for k in range(6*n-5,6*n+1):numerator*=k
        denominator=(3*n-2)*(3*n-1)*(3*n)*(2*n-1)*(2*n)
        assert numerator%denominator==0
        u.append(-numerator//denominator)
        ell.append(n*u[n]-sum(ell[j]*u[n-j] for j in range(1,n)))
        alpha.append(Q(ell[n],n*288**n))
    assert alpha[2:5]==[Q(-5,24),Q(5,16),Q(-1105,1152)]
    return alpha

def logq(q):return log10(abs(q.numerator))-log10(q.denominator)
def description(q):
    p,d=str(q.numerator),str(q.denominator)
    return dict(numerator_digits=len(p.lstrip('-')),denominator_digits=len(d),
                numerator_sha256=hashlib.sha256(p.encode()).hexdigest(),
                denominator_sha256=hashlib.sha256(d.encode()).hexdigest(),
                log10_absolute=logq(q))

def main(limit,out):
    sys.set_int_max_str_digits(0);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    bb=bernoulli_abs(limit);c=constant_coefficients(bb)
    print('Bernoulli constants ready',time.perf_counter()-start,flush=True)
    alpha=leading_airy(limit)
    print('Leading Airy ready',time.perf_counter()-start,flush=True)
    selected=sorted(set([2,3,4,10,25,53,100,200,250,500,750,limit]))
    selected=[g for g in selected if g<=limit]
    exact={};stats={}
    for g in selected:
        a1=linear_coefficient(g,bb,c)
        exact[str(g)]={k:str(v) for k,v in dict(a0=c[g],a1=a1,alpha=alpha[g]).items()}
        stats[str(g)]={k:description(Q(v)) for k,v in exact[str(g)].items()}
        stats[str(g)]['log10_abs_a0_over_alpha']=logq(c[g]/alpha[g])
        print('selected',g,stats[str(g)]['a1']['numerator_digits'],flush=True)
    assert exact['3']['a1']=='17301875/9072'
    assert exact['4']['a1']=='107421875/24192'
    stats['scope']='a0 and leading alpha at every genus; a1 at the explicitly listed genera'
    stats['max_genus']=limit;stats['seconds']=time.perf_counter()-start
    stats['lattice']={'matrix':[[240,0],[-504,6048]],'smith':[gcd(240,gcd(504,6048)),240*6048//24],
        'leading_scale_low_X_counterexample_g3_a0':str(288**2*c[3]),
        'a0_over_alpha_g2':str(c[2]/alpha[2]),'a0_over_alpha_g3':str(c[3]/alpha[3])}
    (out/'scalar_audit.json').write_text(json.dumps(stats,indent=2))
    with gzip.open(out/'scalar_exact.json.gz','wt') as f:json.dump(exact,f)
    with gzip.open(out/'all_constant_and_leading.json.gz','wt') as f:
        json.dump({str(g):{'a0':str(c[g]),'alpha':str(alpha[g])} for g in range(2,limit+1)},f)
    print('PASS',stats['seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--genus',type=int,default=1000)
    p.add_argument('--out',type=Path,default=Path(__file__).parent/'results')
    args=p.parse_args();main(args.genus,args.out)
