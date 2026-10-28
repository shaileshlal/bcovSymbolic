"""Rational local PF frames in x=z/mu, delta=1-x, and y=psi**5.

The PF parameters stay quintic (1/5,2/5,3/5,4/5); mu is a coordinate
scale, never a change of geometry. No analytic-continuation constants enter.
"""
from fractions import Fraction as Q
from functools import lru_cache
from exact_series import Jet

A = tuple(Q(j,5) for j in range(1,5))


def product(xs):
    p=Q(1)
    for x in xs: p*=x
    return p


def mum_periods(n):
    f, h = [Q(1)], [Q(0)]
    logarithmic = Q(0)
    for k in range(1,n):
        f.append(f[-1]*product(k-1+a for a in A)/k**4)
        logarithmic += sum((1/(k-1+a) for a in A), Q(0))-Q(4,k)
        h.append(f[-1]*logarithmic)
    return Jet(dict(enumerate(f)),n), Jet(dict(enumerate(h)),n)


def orbifold_period(exponent,n):
    # omega_r = psi**r sum c_k y**k, y=psi**5.
    c=[Q(1)]
    for k in range(1,n):
        r=exponent+5*k
        c.append(c[-1]*Q((r-5)**4, product(r-j for j in range(1,5))))
    return Jet(dict(enumerate(c)),n)


def _dual_mul(a,b):
    return a[0]*b[0], a[0]*b[1]+a[1]*b[0]


@lru_cache(maxsize=8192)
def conifold_action(r):
    """Coefficients of L(delta**r), and their r derivatives.

    L=theta_x**4-x prod(theta_x+a); theta_x=-(1-delta)d_delta.
    Differentiating in r yields the non-log part of L(delta**r log delta).
    Keys are shifts relative to r. Dual arithmetic keeps everything rational.
    """
    def theta(terms):
        out={}
        for shift,v in terms.items():
            b=_dual_mul(v,(r+shift,Q(1)))
            for j,sign in [(shift,1),(shift-1,-1)]:
                a=out.get(j,(Q(0),Q(0)))
                out[j]=(a[0]+sign*b[0],a[1]+sign*b[1])
        return out
    def add(a,b,scale=1):
        c=dict(a)
        for j,v in b.items():
            w=c.get(j,(Q(0),Q(0)))
            c[j]=(w[0]+scale*v[0],w[1]+scale*v[1])
        return c
    t={0:(Q(1),Q(0))}
    for _ in range(4): t=theta(t)
    p={0:(Q(1),Q(0))}
    for a in A: p=add(theta(p),p,a)
    # x p = (1-delta) p
    return add(add(t,p,-1),{j+1:v for j,v in p.items()})


def conifold_periods(n):
    def extend(c):
        for k in range(3,n):
            power=k-3
            known=sum((v*conifold_action(Q(j)).get(power-j,(0,0))[0]
                       for j,v in enumerate(c)), Q(0))
            lead=conifold_action(Q(k))[-3][0]
            c.append(-known/lead)
        return Jet(dict(enumerate(c)),n)
    # Nonvanishing analytic A-period; initial values fix a convenient rational frame.
    f=extend([Q(1),Q(0),Q(0)])
    # Monodromy image: L(nu log delta) must have no delta**(-1) term.
    b=-conifold_action(Q(1)).get(-2,(0,0))[1]/conifold_action(Q(2))[-3][1]
    nu=extend([Q(0),Q(1),b])
    return f,nu


def generators(omega, theta, theta_t, X):
    u=theta(omega)/omega
    v1=theta(theta_t)/theta_t+2*u
    B2=theta(theta(omega))/omega
    B3=theta(theta(theta(omega)))/omega
    v2=B2-u*v1
    v3=B3+u*v2-5*u*X*(v1-2)
    return (u,v1,v2,v3,X)


@lru_cache(maxsize=12)
def frame(kind,n):
    t=Jet({1:1},n)
    if kind=='mum':
        f,h=mum_periods(n)
        theta=lambda a:-5*a.euler()
        mirror=h/f
        theta_t=-5*(1+mirror.euler())
        X=t/(t-1)
        values=generators(f,theta,theta_t,X)
        factor=f*f*(t-1)/625
        return values,factor,t*mirror.exp()/3125
    if kind=='conifold':
        f,nu=conifold_periods(n)
        theta=lambda a:5*(1-t)*a.derivative()
        flat=nu/f
        values=generators(f,theta,theta(flat),1- t**-1)
        factor=-t*f*f/625
        return values,factor,flat
    if kind=='orbifold':
        f=orbifold_period(1,n)
        h=orbifold_period(2,n)
        # Strip psi powers analytically; no fractional powers in stored jets.
        theta=lambda a:5*a.euler()
        u=1+theta(f)/f
        T=1+theta(h/f)/(h/f)
        v1=T+theta(T)/T+2*u
        B2=theta(u)+u*u
        B3=theta(B2)+u*B2
        v2=B2-u*v1
        X=(1-t)**-1
        v3=B3+u*v2-5*u*X*(v1-2)
        # Full frame factor = psi**(-3) times factor; impose missing psi power separately.
        return (u,v1,v2,v3,X),f*f*(1-t)/625,None
    raise ValueError('Unknown frame')


def coefficients_in_z(jet):
    """Atomic coefficients of jet(z/mu): (z exponent, rational, mu exponent).

    Returns finite certified Laurent data; mu remains an unevaluated nonzero
    coordinate scale. The omitted remainder is O((z/mu)**jet.precision).
    """
    return [(k,c,-k) for k,c in sorted(jet.c.items())]
