"""Physical quintic boundary fixing in the engine's zero-at-v=0 convention.

Conifold gap + orbifold regularity + MUM constant maps + GV/Castelnuovo
constraints. No imported ambiguity tables and no zero-ambiguity fallback.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial, comb
from exact_series import Jet, evaluate_poly
from period_frames import frame, mum_periods
from bcov_yy import X


@lru_cache(maxsize=None)
def bernoulli(n):
    if n==0: return Q(1)
    return -sum((comb(n+1,k)*bernoulli(k) for k in range(n)),Q(0))/(n+1)


def constant_map(g, chi=-200):
    if g<2: raise ValueError('Stable genera only')
    return Q((-1)**g*chi,2)*abs(bernoulli(2*g)*bernoulli(2*g-2))/(
        2*g*(2*g-2)*factorial(2*g-2))


def conifold_coefficient(g):
    # Rational frame nu=delta+7/10 delta^2+...; integral flat coordinate
    # differs by sqrt(5). The monodromy pairing contributes 5**(1-g).
    return Q((-1)**(g-1),5**(g-1))*bernoulli(2*g)/(2*g*(2*g-2))


def required_gv(g,d):
    """Return proven boundary value, or None (never presume unknowns vanish)."""
    if d<=0: raise ValueError('Positive curve degree required')
    if 10*g > d*d+5*d+10: return Q(0)
    if (g,d)==(51,20): return Q(175)
    return None


@lru_cache(maxsize=None)
def gv_kernel(g,h):
    if h>g: return Q(0)
    k=g-h
    S=Jet({j:Q((-1)**j,2**(2*j)*factorial(2*j+1)) for j in range(k+1)},k+1)
    return (S**(2*h-2)).coeff(k)


def divisors(d):
    return (r for r in range(1,d+1) if d%r==0)


def gw_from_gv(g,d,gv):
    return sum((Q(r)**(2*g-3)*gv_kernel(g,h)*gv[h,d//r]
                for r in divisors(d) for h in range(g+1)),Q(0))


def lower_gw(g,d,gv):
    # Under Castelnuovo at these degrees, all same-genus proper divisors
    # vanish as well. Keep that assertion explicit rather than dropping them.
    total=Q(0)
    for r in divisors(d):
        same=required_gv(g,d//r)
        if same is None:
            raise ValueError(f'Missing enumerative constraint n[{g},{d//r}]')
        total+=Q(r)**(2*g-3)*(same+sum(
            (gv_kernel(g,h)*gv[h,d//r] for h in range(g)),Q(0)))
    return total


def q_coefficients(series,degree):
    """Lagrange inversion without expanding z(q): q=x exp(h/f)/3125."""
    f,h=mum_periods(degree+3)
    m=h/f
    derivative=series.derivative()
    out=[series.coeff(0)]
    for d in range(1,degree+1):
        out.append(Q(3125**d,d)*(derivative*((-d*m).exp())).coeff(d-1))
    return out


def seed_gv(degree=20):
    n=degree+8
    values,factor,q=frame('mum',n)
    f,h=mum_periods(n)
    T=1+(h/f).euler()
    x=Jet({1:1},n)
    yukawa=5*((1-x)*f*f*T**3)**-1
    K=q_coefficients(yukawa,degree)
    u,v1,_,_,XX=values
    p=Q(125,12)-v1/2-Q(28,3)*u+Q(5,12)*XX
    derivative=p/(-5*T)
    F1=q_coefficients(derivative,degree)
    gv={}
    for d in range(1,degree+1):
        N0=K[d]/d**3
        gv[0,d]=N0-sum((Q(1,r**3)*gv[0,d//r] for r in divisors(d) if r>1),Q(0))
        N1=F1[d]/d
        known=sum((Q(1,r)*(gv[0,d//r]/12+(gv[1,d//r] if r>1 else 0))
                   for r in divisors(d)),Q(0))
        gv[1,d]=N1-known
    return gv


def record_gv(g,polynomial,degree,gv):
    values,k,_=frame('mum',degree+8)
    F=evaluate_poly(polynomial,values)*k**(g-1)
    coeff=q_coefficients(F,degree)
    for d in range(1,degree+1):
        known=sum((Q(r)**(2*g-3)*gv_kernel(g,h)*gv[h,d//r]
                   for r in divisors(d) for h in range(g+1)
                   if h<g or r>1),Q(0))
        gv[g,d]=coeff[d]-known
    return coeff


class ExactLinearSystem:
    """Incremental exact elimination, with labels on contradictions/rank defects."""
    def __init__(self,size):
        self.size=size
        self.rows={}
        self.labels=[]

    def add(self,row,rhs,label):
        row=[Q(v) for v in row]+[Q(rhs)]
        if len(row)!=self.size+1: raise ValueError('Wrong equation width')
        for pivot,old in sorted(self.rows.items()):
            c=row[pivot]
            if c:
                row=[a-c*b for a,b in zip(row,old)]
        pivot=next((i for i in range(self.size) if row[i]),None)
        if pivot is None:
            if row[-1]: raise ValueError(f'Inconsistent boundary equation: {label}: {row[-1]}')
            return
        c=row[pivot]
        self.rows[pivot]=[a/c for a in row]
        self.labels.append(label)

    def solve(self):
        if len(self.rows)!=self.size:
            raise ValueError(f'Boundary rank {len(self.rows)}/{self.size}; no physical continuation')
        values=[Q(0)]*self.size
        for pivot,row in sorted(self.rows.items(),reverse=True):
            values[pivot]=row[-1]-sum((row[i]*values[i] for i in range(pivot+1,self.size)),Q(0))
        return tuple(values)


def solve_boundary(g,anomaly,gv,precision=None):
    """Solve the actual pure-X vector, rather than imposing a foreign YY basis.

    precision can only increase the conservative jet bound. All needed
    coefficients remain certified even during negative-valuation products.
    """
    if not 2<=g<=53: raise ValueError('This enumerative dataset supports genera 2..53')
    n=max(6*g+16,precision or 0)
    count=3*g-2
    system=ExactLinearSystem(count)
    diagnostics={}
    def rows(kind,indices,target):
        values,k,flat=frame(kind,n)
        K=k**(g-1)
        base=evaluate_poly(anomaly,values)*K
        # Horner-style shared X powers.
        columns=[]
        xp=Jet({0:1},n)
        for j in range(count):
            columns.append(xp*K)
            xp=xp*values[4]
        for i in indices:
            system.add([c.coeff(i) for c in columns],target(i,flat)-base.coeff(i),f'{kind}[{i}]')
        diagnostics[kind]={'rank_after':len(system.rows),'equations':len(indices)}
    orb_count=(3*g-3+4)//5
    rows('orbifold',list(range(orb_count)),lambda i,t:Q(0))
    gap=2*g-2
    # Multiplication by flat**(-gap) is a Laurent jet. Its subleading delta
    # poles are required, even though the gap is stated in the flat coordinate.
    gap_series=frame('conifold',n)[2]**(-gap)*conifold_coefficient(g)
    rows('conifold',list(range(-gap,0)),lambda i,t:gap_series.coeff(i))
    D=(2*g-2)//5
    _,_,q=frame('mum',n)
    target=Jet({0:constant_map(g)},n)
    for d in range(1,D+1): target=target+lower_gw(g,d,gv)*q**d
    rows('mum',list(range(D+1)),lambda i,t:target.coeff(i))
    solution=system.solve()
    diagnostics.update({'genus':g,'unknowns':count,'rank':len(system.rows),
                        'curve_degree':D,'jet_precision':n,
                        'constant_map':str(constant_map(g)),
                        'conifold_coefficient':str(conifold_coefficient(g))})
    return solution,diagnostics
