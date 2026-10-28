"""Direct BCOV in Q[u,v1,v2,v3,X]/(v1,v2,v3,X)^5.

No leading or descendant Airy coefficient enters this calculation.
Seeds: the five low-X coefficients at g=2,3,4 (only four at g=2).
Higher genera: constant maps, degree <=4 GV data from the seeds, and
the known degree <=4 vanishing for g>=5. This is stronger than integrality.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,lcm
import argparse, gzip, json, sys, time
sys.path.insert(0,str(Path(__file__).parent/'engine'))
import bcov_yy as alg
from exact_series import Jet,evaluate_poly
from period_frames import frame
from physical_boundaries import seed_gv,record_gv
from leading_only import bernoulli_abs,constant_coefficients,linear_coefficient

def integer_source(g,vertices,previous_two):
    """Accumulate genus splittings over one common denominator, then reduce.

    This only changes exact rational arithmetic scheduling; no series dressing
    or asymptotic normalization is performed.
    """
    pi,pd=previous_two.integer_form();products=[]
    for h in range(1,g//2+1):
        ai,ad=vertices[h].integer_form();bi,bd=vertices[g-h].integer_form()
        products.append((ai,bi,ad*bd,1 if 2*h==g else 2))
    den=lcm(pd,*(d for _,_,d,_ in products))
    total={m:c*(den//pd) for m,c in pi.items()}
    for ai,bi,d,multiplicity in products:
        partial={}
        for a,ca in ai.items():
            degree_a=sum(a[1:])
            for b,cb in bi.items():
                if degree_a+sum(b[1:])>3:continue
                m=tuple(x+y for x,y in zip(a,b))
                partial[m]=partial.get(m,0)+ca*cb
        factor=multiplicity*(den//d)
        for m,c in partial.items():total[m]=total.get(m,0)+factor*c
    return alg.Poly({m:Q(-c,2*den) for m,c in total.items() if c})

def run(limit,out,reference,comparison=None,native=False):
    sys.set_int_max_str_digits(0);start=time.perf_counter();out.mkdir(parents=True,exist_ok=True)
    alg.J_DEGREE_LIMIT=4
    if native:
        from native_source import NativeSource
        source_function=NativeSource()
    else:source_function=integer_source
    def project(p,d):return alg.Poly({m:c for m,c in p.terms.items() if sum(m[1:])<=d})
    ring=alg.YYRing()
    assert all(sum(m[1:])>=1 for p in ring.table[1:] for m in p.terms)
    full={a['genus']:a for a in json.loads(reference.read_text())['amplitudes']}
    old={r['genus']:r for r in json.loads(comparison.read_text())['rows']} if comparison else {}
    bb=bernoulli_abs(limit);constants=constant_coefficients(bb)
    raw_values,raw_k,raw_q=frame('mum',12)
    values=[v.truncate(5) for v in raw_values];k=raw_k.truncate(5);q=raw_q.truncate(5)
    qp=[q**d for d in range(5)];gv=seed_gv(4)
    vertices={1:Q(125,12)-alg.v1/2-Q(28,3)*alg.u+Q(5,12)*alg.X}
    previous_two=project(ring.covariant(vertices[1],1,1),3)
    all_coeff={};rows=[];seed_gv_data={};last_terms={}
    for g in range(2,limit+1):
        t=time.perf_counter();source=source_function(g,vertices,previous_two)
        assert source.degree('u')<=2,(g,'u degree')
        r0,r1,r2=[source.coefficient('u',j) for j in range(3)]
        gradients=(r0,project(-r1+5*alg.X*r2,3),-r2)
        anomaly=alg.integrate_gradient(*gradients)
        if g<=4:
            coeff=list(map(Q,full[g]['ambiguity'][:5]))
        else:
            # sine kernels for the only nonzero low-degree GV genera h=0,2,3.
            fact=factorial(2*g-2)
            kernels={0:bb[g]/(2*g*fact),2:Q(2*(-1)**g,fact),
                     3:Q((-1)**(g-1)*(2**(2*g-1)-8),fact)}
            target=Jet({0:constants[g]/(-625)**(g-1)},5)
            for d in range(1,5):
                Nd=Q(0)
                for r in range(1,d+1):
                    if d%r:continue
                    Nd+=r**(2*g-3)*sum((kernels[h]*gv[h,d//r] for h in kernels),Q(0))
                target+=Nd*qp[d]
            residual=target/k**(g-1)-evaluate_poly(anomaly,values)
            coeff=[]
            for j in range(5):
                c=(-1)**j*residual.coeff(j);coeff.append(c);residual-=c*values[4]**j
            assert all(residual.coeff(j)==0 for j in range(5)),(g,'MUM')
        polynomial=anomaly+sum((c*alg.X**j for j,c in enumerate(coeff)),alg.Poly())
        resid=project(polynomial.partial('v1')-alg.u*polynomial.partial('v2')
                     -alg.u*(alg.u+5*alg.X)*polynomial.partial('v3')-source,3)
        assert not resid,(g,'anomaly residual')
        if g in full:
            assert anomaly==alg.Poly.from_data(full[g]['anomaly']),(g,'full anomaly')
            assert coeff==list(map(Q,full[g]['ambiguity'][:5])),(g,'full ambiguity')
        if g in old:assert list(map(str,coeff))==old[g]['coefficients'],(g,'prior direct runner')
        assert coeff[0]==constants[g],(g,'constant identity')
        assert coeff[1]==linear_coefficient(g,bb,constants),(g,'linear identity')
        if g<=4:
            record_gv(g,polynomial,4,gv)
            seed_gv_data[str(g)]={str(d):str(gv[g,d]) for d in range(1,5)}
            if g==4:assert all(gv[4,d]==0 for d in range(1,5))
        vertices[g]=project(ring.covariant(project(polynomial,3),g,0),3)
        previous_two=project(ring.covariant(vertices[g],g,1),3)
        all_coeff[str(g)]=list(map(str,coeff))
        rows.append({'genus':g,'coefficients':list(map(str,coeff)),
                     'terms':len(polynomial),'seed':g<=4,'seconds':time.perf_counter()-t})
        last_terms={','.join(map(str,m)):str(c) for m,c in polynomial.terms.items()}
        if g<=5 or g%10==0 or g==limit:
            result={'method':'Direct BCOV J-adic quotient, degree 4; no Airy input',
                    'integer_backend':'GMP native source' if native else 'Python integers',
                    'max_genus':g,'seeds':'g=2,3,4 low-X coefficients',
                    'boundary':'constant maps + seeded degree<=4 GV + g>=5 degree<=4 vanishing',
                    'seed_GV':seed_gv_data,'elapsed_seconds':time.perf_counter()-start,
                    'rows':rows,'last_polynomial_terms':last_terms}
            with gzip.open(out/'direct_J4_exact.json.gz','wt') as f:json.dump(result,f)
            (out/'direct_status.json').write_text(json.dumps({k:v for k,v in result.items() if k not in ['rows','last_polynomial_terms']},indent=2))
            print('g',g,'terms',len(polynomial),'seconds',round(rows[-1]['seconds'],2),'total',round(time.perf_counter()-start,2),flush=True)
        if (out/'STOP').exists():break
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--genus',type=int,default=200)
    p.add_argument('--out',type=Path,default=Path(__file__).parent/'results')
    p.add_argument('--reference',type=Path,default=Path(__file__).parent/'reference_amplitudes.json')
    p.add_argument('--comparison',type=Path)
    p.add_argument('--native',action='store_true',help='Use optional integer-only GMP accelerator')
    a=p.parse_args();run(a.genus,a.out,a.reference,a.comparison,a.native)
