#!/usr/bin/env python3
"""Exact quintic Yamaguchi--Yau P2/P3 basis transforms.

Yamaguchi--Yau (hep-th/0406078), equations (3.33), (4.1), (4.8), (A.1).
Uses their X=1/(1-psi**5) and D=psi*d/dpsi. No floats.

1. Uncentered Bell/Wronskian generators J1,J2,J3.
2. Centered H1,H2,H3 with preserved X-only projection.
3. Pure X ambiguity slices, conifold P=1-X reexpansion.
4. Fourier-dual basis via exact quintic mirror q series.

Run: python yy_p2_p3_exact.py
Requires: sympy; Python 3.10+.
"""
import sympy as s
from fractions import Fraction as Q
from math import factorial
v1,v2,v3,X = s.symbols('v1 v2 v3 X')
j1,j2,j3,Y = s.symbols('J1 J2 J3 Y')
H1,H2,H3 = s.symbols('H1 H2 H3')
v1, v2, v3, X = s.symbols('v1 v2 v3 X')

Dmap = {v1: -v1 ** 2 - 2 * v2 - 10 * X + 5 * v1 * X, v2: -v1 * v2 + v3, v3: v2 ** 2 - 24 * X - 25 * v2 * X - 5 * v1 * v2 * X + 10 * v3 * X, X: 5 * X * (X - 1)}

def ring_deriv(f):
    return s.expand(sum((s.diff(f, v) * d for v, d in Dmap.items())))

J1 = v1

J2 = s.expand(J1 * J1 + ring_deriv(J1))

J3 = s.expand(J1 * J2 + ring_deriv(J2))

j1, j2, j3, Y = s.symbols('J1 J2 J3 Y')

v2inv = s.solve(J2 - j2, v2)[0].subs(v1, j1)

v3inv = s.solve((J3 - j3).subs({v1: j1, v2: v2inv}), v3)[0]

P2 = s.Rational(3125, 144) - s.Rational(15625, 288) * v1 + s.Rational(125, 24) * v1 ** 2 - s.Rational(5, 24) * v1 ** 3 - s.Rational(3125, 36) * v2 + s.Rational(25, 6) * v1 * v2 + s.Rational(350, 9) * v3 - s.Rational(28795, 144) * X - s.Rational(835, 144) * v1 * X + s.Rational(5, 6) * v1 ** 2 * X - s.Rational(2375, 12) * v2 * X + s.Rational(205, 144) * X ** 2 - s.Rational(325, 288) * v1 * X ** 2 + s.Rational(25, 48) * X ** 3

P2J = s.Poly(s.expand(P2.subs({v1: j1, v2: v2inv, v3: v3inv}, simultaneous=True)), j1, j2, j3, X).as_expr()

c0 = 781250 - 2734375 * v1 ** 3 + 787500 * v1 ** 4 - 94500 * v1 ** 5 + 4536 * v1 ** 6 + 6562500 * v2 - 16721250 * v1 * v2 - 2625000 * v1 ** 2 * v2 + 819000 * v1 ** 3 * v2 - 54432 * v1 ** 4 * v2 - 18112500 * v2 ** 2 - 1772400 * v1 * v2 ** 2 + 295344 * v1 ** 2 * v2 ** 2 - 936320 * v2 ** 3 - 4935000 * v3 + 12337500 * v1 * v3 - 1184400 * v1 ** 2 * v3 + 47376 * v1 ** 3 * v3 + 19740000 * v2 * v3 - 947520 * v1 * v2 * v3 - 4421760 * v3 ** 2

c1 = 27683000 - 72635850 * v1 + 12252135 * v1 ** 2 - 3366615 * v1 ** 3 + 604044 * v1 ** 4 - 41580 * v1 ** 5 - 81544680 * v2 - 54284034 * v1 * v2 + 3202584 * v1 ** 2 * v2 + 93240 * v1 ** 3 * v2 - 99165864 * v2 ** 2 + 3824016 * v1 * v2 ** 2 + 45473064 * v3 + 1318632 * v1 * v3 - 189504 * v1 ** 2 * v3 + 45007200 * v2 * v3

c2 = -112828006 - 12527550 * v1 + 5722185 * v1 ** 2 - 1658685 * v1 ** 3 + 176400 * v1 ** 4 - 233375520 * v2 - 3865134 * v1 * v2 + 104160 * v1 ** 2 * v2 - 113818740 * v2 ** 2 - 323736 * v3 + 256620 * v1 * v3

c3 = 3339968 - 4795350 * v1 + 2353785 * v1 ** 2 - 444325 * v1 ** 3 + 819840 * v2 - 266910 * v1 * v2 - 118440 * v3

c4 = 1696500 - 1683150 * v1 + 686175 * v1 ** 2 + 119700 * v2

c5 = 477000 - 598500 * v1

c6 = 225000

P3 = s.Rational(5, 72576) * sum((c * X ** k for k, c in enumerate([c0, c1, c2, c3, c4, c5, c6])))

P3J = s.Poly(s.expand(P3.subs({v1: j1, v2: v2inv, v3: v3inv}, simultaneous=True)), j1, j2, j3, X).as_expr()

H1, H2, H3 = s.symbols('H1 H2 H3')

h2 = J2 + 10 * X

h3 = J3 - 50 * X + 100 * X ** 2

invs = {v1: H1, v2: (5 * X * H1 - H2) / 2, v3: (-H3 + 5 * X * H2 + (25 * X ** 2 - 35 * X) * H1) / 2}

P3H = s.Poly(s.expand(P3.subs(invs, simultaneous=True)), H1, H2, H3, X).as_expr()

# Exact truncated Taylor-series operations, with xi=z/(1-5^5 z).
def zero(n):
    return [Q(0)] * n

def add(a, b):
    return [x + y for x, y in zip(a, b)]

def scale(a, c):
    return [x * c for x in a]

def mul(a, b):
    n = len(a)
    c = zero(n)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:n - i]):
                if y:
                    c[i + j] += x * y
    return c

def pw(a, k):
    n = len(a)
    o = zero(n)
    o[0] = 1
    for _ in range(k):
        o = mul(o, a)
    return o

def inv(a):
    n = len(a)
    o = zero(n)
    o[0] = 1 / a[0]
    for i in range(1, n):
        o[i] = -sum((a[j] * o[i - j] for j in range(1, i + 1))) / a[0]
    return o

def exp(a):
    assert a[0] == 0
    n = len(a)
    o = zero(n)
    o[0] = 1
    for i in range(1, n):
        o[i] = sum((Q(j) * a[j] * o[i - j] for j in range(1, i + 1))) / i
    return o

def compose(a, b):
    n = len(a)
    ret = zero(n)
    pow_b = zero(n)
    pow_b[0] = 1
    for c in a:
        ret = add(ret, scale(pow_b, c))
        pow_b = mul(pow_b, b)
    return ret

def mirror_data(m):
    n = m + 1
    z = zero(n)
    for i in range(1, n):
        z[i] = Q((-3125) ** (i - 1))
    t = zero(n)
    h = zero(n)
    for i in range(n):
        a = Q(factorial(5 * i), factorial(i) ** 5)
        t[i] = a
        if i:
            h[i] = 5 * a * sum((Q(1, j) for j in range(i + 1, 5 * i + 1)), Q(0))
    w = compose(t, z)
    e = compose(h, z)
    logr = mul(e, inv(w))
    q = mul(z, exp(logr))
    zz = scale(z, -3125)
    zz[0] += 1
    Ebase = mul(mul(w, w), zz)
    return (q, Ebase)

def dual_polynomials(g, m, q, Ebase):
    n = m + 1
    E = pw(Ebase, g - 1)
    Ei = inv(E)
    Psi = [mul(pw(q, k), Ei) for k in range(n)]
    assert all((Psi[k][i] == (1 if i == k else 0) for k in range(n) for i in range(k + 1)))
    return Psi

def coordinates_in_dual(poly_coeff, Psi):
    m = len(poly_coeff) - 1
    out = []
    rem = list(poly_coeff)
    for k in range(m + 1):
        c = rem[k]
        out.append(c)
        rem = add(rem, scale(Psi[k], -c))
    assert all((x == 0 for x in rem))
    return out

def fmt(a):
    return str(a) if a.denominator != 1 else str(a.numerator)

def deriv(a):
    n = len(a)
    out = zero(n)
    for i in range(1, n):
        out[i - 1] = i * a[i]
    return out

def Dpsi(a):
    n = len(a)
    xi = zero(n)
    xi[1] = 1
    coeff = mul(scale(xi, -5), add([Q(1)] + zero(n - 1), scale(xi, 3125)))
    return mul(coeff, deriv(a))

def revert(q):
    n = len(q)
    x = zero(n)
    x[1] = 1
    for m in range(2, n):
        x[m] = -compose(q, x)[m]
    assert compose(q, x) == [Q(0), Q(1)] + zero(n - 2)
    return x

def hol_yy(g, m):
    n = m + 1
    q, Ebase = mirror_data(m)
    z = zero(n)
    for i in range(1, n):
        z[i] = Q((-3125) ** (i - 1))
    w0raw = [Q(factorial(5 * i), factorial(i) ** 5) for i in range(n)]
    w = compose(w0raw, z)
    h = q[1:] + [Q(0)]
    xi = zero(n)
    xi[1] = 1
    dt = mul(add([Q(1)] + zero(n - 1), scale(xi, 3125)), add([Q(1)] + zero(n - 1), mul(xi, mul(deriv(h), inv(h)))))
    b = mul(Dpsi(w), inv(w))
    v1 = add(scale(b, 2), mul(Dpsi(dt), inv(dt)))
    v2 = add(mul(Dpsi(Dpsi(w)), inv(w)), scale(mul(v1, b), -1))
    X = scale(xi, -3125)
    v3 = add(mul(Dpsi(Dpsi(Dpsi(w))), inv(w)), scale(mul(b, add(add(scale(v2, -1), scale(mul(X, v1), 5)), scale(X, -10))), -1))
    E = pw(Ebase, g - 1)
    N0 = Q((-1) ** (g - 1), 625 ** (g - 1))
    return (q, [v1, v2, v3, X], scale(E, N0))

def eval_poly(pol, gens):
    n = len(gens[0])
    V = s.symbols('v1 v2 v3 X')
    r = zero(n)
    for exps, coeff in s.Poly(pol, *V).terms():
        t = zero(n)
        t[0] = Q(coeff.p, coeff.q)
        for val, e in zip(gens, exps):
            t = mul(t, pw(val, e))
        r = add(r, t)
    return r

def main():
    assert s.expand(P2J.subs({j1:v1,j2:J2,j3:J3}, simultaneous=True)-P2)==0
    assert s.expand(P3J.subs({j1:v1,j2:J2,j3:J3}, simultaneous=True)-P3)==0
    assert s.expand(P3H.subs({H1:v1,H2:h2,H3:h3},simultaneous=True)-P3)==0
    P2H=s.expand(P2.subs(invs,simultaneous=True))
    assert s.expand(P2H.subs({H1:v1,H2:h2,H3:h3},simultaneous=True)-P2)==0
    print('J2:',J2,'\nJ3:',J3)
    print('H2:',h2,'\nH3:',h3)
    for g,P,PJ,PH,m in ((2,P2,P2J,P2H,3),(3,P3,P3J,P3H,6)):
        projection=P.subs({v1:0,v2:0,v3:0})
        projection_j=PJ.subs({j1:0,j2:0,j3:0})
        assert s.expand(PH.subs({H1:0,H2:0,H3:0})-projection)==0
        print(f'\nGenus {g}:')
        print('monomials original, J, centered H:', *[len(s.Poly(x,*sy).terms()) for x,sy in ((P,(v1,v2,v3,X)),(PJ,(j1,j2,j3,X)),(PH,(H1,H2,H3,X)))])
        print('original X-only projection:',s.factor(projection))
        print('uncentered J-induced projection shift:',s.factor(projection_j-projection))
        print('X-only projection at P=1-X:',s.expand(projection.subs(X,1-Y)))
        print('H3 squared coefficient:',s.factor(s.Poly(PH,H3).coeff_monomial(H3**2)))
        q,Ebase=mirror_data(m)
        Psi=dual_polynomials(g,m,q,Ebase)
        xi=s.symbols('xi')
        fxi=s.Poly(s.expand(projection.subs(X,-3125*xi)),xi)
        fxi_coeff=[Q(fxi.nth(j).p,fxi.nth(j).q) for j in range(m+1)]
        b=coordinates_in_dual(fxi_coeff,Psi)
        print('Fourier-dual coordinates b[g,k]:',list(map(str,b)))
        print('Dual Psi[g,k], first two:', [list(map(str,p)) for p in Psi[:2]])
        q_full,gens,N=hol_yy(g,5)
        f_q=compose(mul(N,eval_poly(P,gens)),revert(q_full))
        print('Full A-model F_g coefficients (reliable through q^4):',list(map(str,f_q[:5])))
    # Known quintic degree-one and degree-two GV sanity checks
    f2q=compose(mul(hol_yy(2,5)[2],eval_poly(P2,hol_yy(2,5)[1])),revert(hol_yy(2,5)[0]))
    assert f2q[0]==Q(-5,144) and f2q[1]==Q(575,48) and f2q[2]==Q(5125,2)
    f3q=compose(mul(hol_yy(3,5)[2],eval_poly(P3,hol_yy(3,5)[1])),revert(hol_yy(3,5)[0]))
    assert f3q[0]==Q(5,36288) and f3q[1]==Q(2875,6048)
    print('\nAll exact identities and enumerative sanity checks PASSED.')

if __name__=='__main__': main()
