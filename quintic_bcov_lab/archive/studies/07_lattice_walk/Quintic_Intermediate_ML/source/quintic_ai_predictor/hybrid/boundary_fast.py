"""Exact physical boundary rows for the mirror-quintic YY ambiguity.

Each row is Q-linear in a_g[0:3g-2]. The B-model conifold, orbifold and
MUM constraints are independent of ML. Unknown GV values must remain free;
never substitute n_{g,d}=0 unless geometric vanishing is established.
"""
from __future__ import annotations
from fractions import Fraction as Q
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'engine'))
from physical_boundaries import (Jet,frame,conifold_coefficient,constant_map,
                                  lower_gw,ExactLinearSystem)
from exact_series import evaluate_poly

@dataclass
class AffineBoundary:
    genus:int
    particular:tuple[Q,...]
    kernel:tuple[tuple[Q,...],...]
    free_indices:tuple[int,...]
    constraint_rows:tuple[tuple[tuple[Q,...],Q,str],...]
    @property
    def rank(self):return len(self.particular)-len(self.free_indices)
    def evaluate(self,free_values=()):
        if len(free_values)!=len(self.free_indices):raise ValueError('one value required per free slot')
        ans=list(self.particular)
        for val,v in zip(free_values,self.kernel):
            for i in range(len(ans)):ans[i]+=Q(val)*v[i]
        return tuple(ans)
    def verify(self,x):
        if len(x)!=len(self.particular):return False
        return all(sum((a*b for a,b in zip(row,x)),Q(0))==rhs for row,rhs,_ in self.constraint_rows)
    def project_by_free_coordinates(self,candidate):
        """Set free coefficients from a *rational* guess; solve all pivots exactly.

        This is an algebraic projection, NOT the Euclidean nearest projection.
        Decimal candidate coordinates may be converted via Q(str(...)), but
        such conversion gives exact satisfaction, not physical validity.
        """
        return self.evaluate([Q(candidate[j]) for j in self.free_indices])


def boundary_system(g,anomaly,gv=None,include=('orbifold','conifold','mum'),mum_degrees=None,precision=None):
    """Return an exact affine ambiguity family. By default enforce all known MUM degrees.

    For an explicitly partial MUM test use mum_degrees=[0] to impose constant
    maps only. If mum is omitted no MUM input is consulted.
    """
    if g<2:raise ValueError('g >=2 required')
    include=set(include)
    num=3*g-2
    n=max(3*g+8,precision or 0)
    D=(2*g-2)//5
    indices=list(range(D+1)) if mum_degrees is None else list(mum_degrees)
    if any(d<0 or d>D for d in indices):raise ValueError('MUM degrees must be <= floor((2g-2)/5)')
    rows=[]
    for frame_name in ('orbifold','conifold','mum'):
        if frame_name not in include:continue
        vals,k,flat=frame(frame_name,n)
        K=k**(g-1);base=evaluate_poly(anomaly,vals)*K
        cols=[];xp=Jet({0:1},n)
        for j in range(num):
            cols.append(xp*K);xp=xp*vals[4]
        if frame_name=='orbifold':
            ii=range((3*g-3+4)//5)
            target=lambda i:Q(0)
        elif frame_name=='conifold':
            gap=2*g-2
            ii=range(-gap,0)
            target_series=flat**(-gap)*conifold_coefficient(g)
            target=lambda i:target_series.coeff(i)
        else:
            _,_,q=frame('mum',n)
            ii=indices
            target_series=Jet({0:constant_map(g)},n)
            for d in indices:
                if d>0:
                    if gv is None:raise ValueError('Missing lower-genus GV dictionary for nonconstant MUM constraints')
                    target_series+=lower_gw(g,d,gv)*q**d
            target=lambda i:target_series.coeff(i)
        for i in ii:
            rows.append((tuple(c.coeff(i) for c in cols),Q(target(i)-base.coeff(i)),f'{frame_name}:{i}'))
    echelon=ExactLinearSystem(num)
    for row,rhs,label in rows:echelon.add(row,rhs,label)
    pivots=sorted(echelon.rows)
    free=tuple(i for i in range(num) if i not in echelon.rows)
    def solve_with(free_values,rhs=True):
        values=[Q(0)]*num
        for j,v in zip(free,free_values):values[j]=Q(v)
        for p in reversed(pivots):
            row=echelon.rows[p]
            values[p]=(row[-1] if rhs else Q(0))-sum((row[j]*values[j] for j in range(p+1,num)),Q(0))
        return tuple(values)
    p=solve_with([0]*len(free))
    k=[]
    for i in range(len(free)):
        x=[0]*len(free);x[i]=1
        k.append(solve_with(x,False))
    result=AffineBoundary(g,p,tuple(k),free,tuple(rows))
    assert result.verify(p)
    assert all(all(sum((r[j]*v[j] for j in range(num)),Q(0))==0 for r,_,_ in rows) for v in k)
    return result
