from fractions import Fraction as Q
import json
from pathlib import Path
from bcov_yy import Poly,reference_genus2, YYRing,u,v1,v2,v3,X
from quintic_frames import *
from math import comb

EXPECTED_HKQ=(Q(25,144)-Q(625,288)*v1+Q(25,24)*v1**2-Q(5,24)*v1**3
              -Q(625,36)*v2+Q(25,6)*v1*v2+Q(350,9)*v3
              -Q(5759,3600)*X-Q(167,720)*v1*X+Q(1,6)*v1**2*X
              -Q(475,12)*v2*X+Q(41,3600)*X**2-Q(13,288)*v1*X**2+X**3/240)
assert yy_to_hkq(reference_genus2(),2)==EXPECTED_HKQ
assert hkq_to_yy(EXPECTED_HKQ,2)==reference_genus2()
assert z_to_yy(yy_to_z(reference_genus2(),2),2)==reference_genus2()
for g in (2,3,6):
 p=(Q(4,27)+u*v1+v2*v3/29+X**(3*g-3)+v1*X**2)
 assert hkq_to_yy(yy_to_hkq(p,g),g)==p
 assert z_to_hkq(hkq_to_z(p,g),g)==p
 # Derivation is coordinate-covariant: theta_H X=X(X-1), theta_z Xz=Xz(Xz-1)
ring=YYRing()
assert theta_hkq(X,ring)==X*(X-1)
assert theta_z(X,ring)==X*(X-1)
# The YY Euler operator acting on v1 maps through theta_H scaling.
assert theta_hkq(v1,ring)== yy_to_hkq(ring.theta(v1),2)*5 # YY theta v1Y divided by 5*5
obj=FramedPolynomial(reference_genus2(),2,'yy')
assert obj.to('hkq').poly==EXPECTED_HKQ
assert obj.to('small_z').to('yy')==obj
try:
    obj+obj.to('hkq')
except ValueError:pass
else:raise AssertionError('Frame mix accepted')
assert json.loads((Path(__file__).resolve().parent/'CONVENTIONS.json').read_text())['gu_gauge']==GU_GAUGE
print('PASS: HKQ g2, exact frame changes, typed frame guard, Euler differential ring')
