"""Exact coordinate/normalization maps for one-parameter quintic BCOV.

YY: original Yamaguchi--Yau fifth-root coordinate s, psi=s^5.
HKQ: psi=1/(3125*z), z is the small large-radius coordinate.
The symbol X in HKQ/YY is 1/(1-psi); in small-z X_z=1/(1-3125*z).
P_g=C^(g-1) F_g. Output Poly uses the variables of the declared target frame.

No periods, q-series, or hidden normalizations are applied here.
"""
from __future__ import annotations
from fractions import Fraction as Q
from dataclasses import dataclass
from typing import Literal
from bcov_yy import Poly, YYRing, u,v1,v2,v3,X


def _genus(g:int):
    if type(g) is not int or g<2:raise ValueError('Normalized polynomial conversion requires integer genus >= 2')


def yy_to_hkq(p:Poly,g:int)->Poly:
    _genus(g)
    # uY=5uH, v1Y=5v1H, v2Y=25v2H, v3Y=125v3H, X unchanged
    return Poly.coerce(p).substitute({'u':5*u,'v1':5*v1,'v2':25*v2,'v3':125*v3})/5**(3*g-3)


def hkq_to_yy(p:Poly,g:int)->Poly:
    _genus(g)
    return (5**(3*g-3)*Poly.coerce(p).substitute({'u':u/5,'v1':v1/5,'v2':v2/25,'v3':v3/125}))


def hkq_to_z(p:Poly,g:int)->Poly:
    _genus(g)
    return (-1)**(g-1)*Poly.coerce(p).substitute({'u':-u,'v1':-v1,'v2':v2,'v3':-v3,'X':1-X})


def z_to_hkq(p:Poly,g:int)->Poly:
    return hkq_to_z(p,g)  # involution including the global P_g parity


def yy_to_z(p:Poly,g:int)->Poly:
    return hkq_to_z(yy_to_hkq(p,g),g)


def z_to_yy(p:Poly,g:int)->Poly:
    return hkq_to_yy(z_to_hkq(p,g),g)


def theta_hkq(p:Poly,ring:YYRing|None=None)->Poly:
    """HKQ Euler derivation, obtained by conjugation from original YY theta.
    This does not rescale amplitudes by genus; it transforms scalar connection generators.
    """
    ring=ring or YYRing()
    yy=Poly.coerce(p).substitute({'u':u/5,'v1':v1/5,'v2':v2/25,'v3':v3/125})
    return ring.theta(yy).substitute({'u':5*u,'v1':5*v1,'v2':25*v2,'v3':125*v3})/5


def theta_z(p:Poly,ring:YYRing|None=None)->Poly:
    """Small-z Euler derivation via HKQ coordinate inversion (theta_z=-theta_H)."""
    hkq=Poly.coerce(p).substitute({'u':-u,'v1':-v1,'v2':v2,'v3':-v3,'X':1-X})
    dh=theta_hkq(hkq,ring)
    return -dh.substitute({'u':-u,'v1':-v1,'v2':v2,'v3':-v3,'X':1-X})

# Gu/Klemm propagator gauge. These values are a choice, not universal tensors.
GU_GAUGE={
  'q^z_zz':'-8/(5*z)', 'q_z^zz':'z/25', 'q_z^z':'-6/625',
  'q_zz':'2/(25*z^2)', 'q_z':'2/(3125*z)'
}
# In the small-z frame: Szz=-z^2/(5*Xz)*(v1+3/5),
# tSz= z/(5*Xz)*(v2+2/25),
# tS= 1/(10*Xz)*(v3+(1-Xz)*(v2+4/125)).


Frame = Literal["yy", "hkq", "small_z"]

@dataclass(frozen=True)
class FramedPolynomial:
    """Type-safe normalized physical P_g for one fixed coordinate convention.

    The same symbolic name X denotes DIFFERENT functions in HKQ and small-z.
    Conversion is explicit and the object refuses silent cross-frame addition.
    """
    poly: Poly
    genus: int
    frame: Frame

    def __post_init__(self):
        _genus(self.genus)
        if self.frame not in ('yy','hkq','small_z'):
            raise ValueError(f'Unknown BCOV polynomial frame: {self.frame!r}')
        if not isinstance(self.poly, Poly):
            raise TypeError('Polynomial data must be an exact Poly')

    def to(self,target:Frame)->'FramedPolynomial':
        if target not in ('yy','hkq','small_z'):
            raise ValueError(f'Unknown BCOV polynomial frame: {target!r}')
        if target==self.frame:return self
        conv={
            ('yy','hkq'):yy_to_hkq, ('hkq','yy'):hkq_to_yy,
            ('hkq','small_z'):hkq_to_z, ('small_z','hkq'):z_to_hkq,
            ('yy','small_z'):yy_to_z, ('small_z','yy'):z_to_yy
        }
        return FramedPolynomial(conv[(self.frame,target)](self.poly,self.genus),self.genus,target)

    def __add__(self,other:'FramedPolynomial')->'FramedPolynomial':
        if not isinstance(other,FramedPolynomial):return NotImplemented
        if (self.genus,self.frame)!=(other.genus,other.frame):
            raise ValueError('Cannot add polynomials with different genus or frame')
        return FramedPolynomial(self.poly+other.poly,self.genus,self.frame)
