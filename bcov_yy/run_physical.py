"""Checkpointed physical recursion, target genus 53 by default.

Use --plan first, then benchmark a small --genus. Full genus 53 is a target,
not a runtime/memory promise for this expanded-polynomial Python backend.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter
from fractions import Fraction as Q
import bcov_yy as algebra
from bcov_yy import BCOV
from physical_boundaries import seed_gv, solve_boundary, record_gv, required_gv
from monodromy import validate


def plan(max_genus=53):
    if not 2<=max_genus<=53: raise ValueError("Supported target genera are 2..53")
    result=[]
    for g in range(2,max_genus+1):
        d=(2*g-2)//5
        result.append({'genus':g,'unknowns':3*g-2,'conifold':2*g-2,
                       'orbifold':(3*g+1)//5,'mum':d+1,'gv_max_degree':d,
                       'jet_precision':6*g+16,
                       'special_gv':{'degree':20,'value':175} if g==51 else None})
    return result


def atomic_json(path,data):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2)+'\n')
    tmp.replace(path)


def run(max_genus,output,mu='mu',max_terms=1000000,max_products=50000000,degree=None):
    if not 2<=max_genus<=53: raise ValueError('Genus must lie in 2..53')
    validate()
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    algebra.set_resource_limits(max_terms,max_products)
    degree=max((2*max_genus-2)//5,degree or 0)
    cp=output/'amplitudes.json'
    state=output/'run.json'
    # Loading validates all previous anomaly polynomials, not just a checksum.
    s=BCOV.load(cp) if cp.exists() else BCOV()
    gv=seed_gv(degree)
    diagnostics=[]
    for g,a in sorted(s.amplitudes.items()):
        if a.ambiguity is None: break
        if s.residual(g): raise ValueError(f'Checkpoint residual at genus {g}')
        physical,_=solve_boundary(g,a.anomaly,gv)
        stored=tuple(a.ambiguity)+(Q(0),)*(len(physical)-len(a.ambiguity))
        if physical!=stored:
            raise ValueError(f'Checkpoint fails physical boundaries at genus {g}')
        record_gv(g,a.polynomial(),degree,gv)
        if any(gv[g,d].denominator!=1 for d in range(1,degree+1)):
            raise ValueError(f'Checkpoint has nonintegral GV data at genus {g}')
    status={'schema':1,'target_genus':max_genus,'mu':mu,'coordinate':'x=z/mu',
            'physical_q':'x*exp(h/f)/3125','gv_degree':degree,
            'backend':'Python Fraction; expanded sparse polynomials',
            'max_terms':max_terms,'max_products':max_products,
            'completed_genus':max((g for g,a in s.amplitudes.items() if a.ambiguity is not None),default=1)}
    atomic_json(output/'plan.json',plan(max_genus))
    try:
        for g in range(2,max_genus+1):
            if g in s.amplitudes and s.amplitudes[g].ambiguity is not None: continue
            start=perf_counter()
            a=s.integrate(g)
            vector,diag=solve_boundary(g,a.anomaly,gv)
            s.set_ambiguity(g,vector)
            p=s.amplitudes[g].polynomial()
            if s.residual(g): raise ValueError(f'BCOV residual at genus {g}')
            record_gv(g,p,degree,gv)
            for d in range(1,degree+1):
                if gv[g,d].denominator!=1:
                    raise ValueError(f'Nonintegral GV invariant at ({g},{d})')
                known=required_gv(g,d)
                if known is not None and gv[g,d]!=known:
                    raise ValueError(f'Enumerative residual at ({g},{d})')
            diag.update({'terms':len(p),'seconds':perf_counter()-start,
                         'max_coefficient_bits':max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())
                                                    for c in p.terms.values())})
            diagnostics.append(diag)
            # Commit only fully validated physical amplitudes.
            s.save(cp)
            atomic_json(output/f'boundary_g{g:02d}.json',diag)
            atomic_json(output/'gv.json',[{'g':h,'d':d,'value':str(v)} for (h,d),v in sorted(gv.items())])
            status.update({'completed_genus':g,'status':'running'})
            atomic_json(state,status)
            print(f"g={g} terms={len(p)} boundary_rank={diag['rank']} time={diag['seconds']:.3f}s",flush=True)
            s.compact_recipes()
        status['status']='complete'
    except (Exception,KeyboardInterrupt) as error:
        status.update({'status':'stopped','error':str(error),'exception':type(error).__name__})
        atomic_json(state,status)
        raise
    atomic_json(state,status)
    return s,gv


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genus',type=int,default=53)
    parser.add_argument('--output',default='physical_run')
    parser.add_argument('--mu',default='mu',help='Nonzero coordinate scale; symbolic label, not a numerical substitution')
    parser.add_argument('--max-terms',type=int,default=1000000,help='Per-operation term limit; 0 disables')
    parser.add_argument('--max-products',type=int,default=50000000,help='Polynomial pair-product limit; 0 disables')
    parser.add_argument('--gv-degree',type=int,help='Additional GV degrees for verification')
    parser.add_argument('--plan',action='store_true')
    args=parser.parse_args()
    if args.plan:
        print(json.dumps(plan(args.genus),indent=2))
    else:
        run(args.genus,args.output,args.mu,args.max_terms,args.max_products,args.gv_degree)
