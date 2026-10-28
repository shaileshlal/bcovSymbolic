"""Exact genus-54 BCOV ambiguity lattice response and MUM q-expansion.

This is a homogeneous/linear BCOV response, not a recomputation of the
complete physical genus-54 anomaly. No held-out GV/PT table is imported.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
from sympy import symbols,expand,simplify
import sys,json
sys.path.insert(0,str(Path(__file__).parent/'code'/'engine'))
from period_frames import mum_periods
from physical_boundaries import q_coefficients

Q=21;g=54;n=319125
m=(3*g-3+4)//5
D=(2*g-2)//5
assert m==32 and D==21
X,x=symbols('X x')
# Delta a(X) = integer*5^107 X^21(1-X)^32.
scale_5=4*g-4-5*D
assert scale_5==107
# The reduced section f=(x-1)^(g-1) a(X) is n*5^107 x^21.
response=(x-1)**(g-1)*(X**D*(1-X)**m).subs(X,x/(x-1))
assert simplify(response-x**D)==0
# Zero constant term, conifold gap and orbifold are automatic in the homogeneous subspace.
# MUM response: Delta F_g=n*q^D * f(x)^(2g-2)*exp(-D*h(x)/f(x))
# which is q^D times an ordinary power series in the large-volume mirror map.
F,H=mum_periods(8)
R=(F**(2*g-2))*((-D*H/F).exp())
vals=q_coefficients(R,4)
assert vals[0]==1
assert vals[1]==240*(g-1)-770*D == -3450
result={
 'g':g,'degree':D,'n_predicted':n,
 'YY_homogeneous_difference':f'{n}*5^{scale_5}*X^{D}*(1-X)^{m}',
 'reduced_section_difference':f'{n}*5^{scale_5}*x^{D}',
 'physical_MUM_q_response':'Delta F_g = n*q^D * omega0(z(q))^(2g-2) * (z(q)/q)^D',
 'q_shifts_per_unit_charge':list(map(str,vals)),
 'q_shifts_for_predicted_319125':list(map(str,[v*n for v in vals])),
 'limitation':'This certifies the ambiguity response/section and MUM behavior, NOT full physical BCOV genus54; the charge is not fixed by the anomaly alone.'
}
(Path(__file__).parent/'results').mkdir(exist_ok=True)
(Path(__file__).parent/'results'/'bcov_g54_lattice_response.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
