"""Shared exact BCOV/YY backend, or a user-supplied compatible BCOV package.

The bundled engine is pinned to the exact YY code used by the historical
quintic studies; this is NOT an assertion that an unavailable local bcov0.2
installation has been loaded. An alternative backend is opt-in and validated.
"""
from __future__ import annotations
import importlib
import os
import sys
from functools import lru_cache
from pathlib import Path

EXPECTED_CONVENTION = 'YY2004-quintic-psi-X=1/(1-psi^5)-v3-unscaled-v1'

@lru_cache(maxsize=1)
def module():
    name = os.environ.get('BCOV_BACKEND_MODULE', '').strip()
    if name:
        result = importlib.import_module(name)
    else:
        # Legacy modules import their internal companions by unqualified names.
        engine_path = Path(__file__).resolve().parent / 'vendor' / 'engine'
        if str(engine_path) not in sys.path:
            sys.path.insert(0, str(engine_path))
        result = importlib.import_module('bcov_yy')
    for required in ('Poly', 'BCOV', 'YYRing', 'reference_genus2'):
        if not hasattr(result, required):
            raise TypeError(f'BCOV backend {result.__name__} lacks {required}')
    if result.BCOV.convention != EXPECTED_CONVENTION:
        raise ValueError('Unsupported BCOV convention: '+str(result.BCOV.convention))
    return result


def api():
    m = module()
    return {name:getattr(m,name) for name in ('BCOV','Poly','YYRing','reference_genus2','YY_GENUS2_AMBIGUITY')}


def validate_small():
    """Independent genus-two anomaly and published reference regression."""
    m=module()
    b=m.BCOV()
    a=b.integrate(2)
    b.set_ambiguity(2,m.YY_GENUS2_AMBIGUITY)
    assert b.amplitudes[2].polynomial() == m.reference_genus2()
    assert b.residual(2) == 0
    return {'backend':m.__name__, 'convention':m.BCOV.convention,
            'genus_two_terms':len(b.amplitudes[2].polynomial()),'bcov_residual_zero':True}
