"""Adapter to bcovSymbolic's exact YY engine with explicit convention checks.

Historical subprocesses keep their archived engine. Only the unified analysis
API uses this adapter; BCOV_BACKEND_MODULE can explicitly select a compatible
alternative for comparison.
"""
from __future__ import annotations

import importlib
import os
from functools import lru_cache
from types import ModuleType
from typing import Any

EXPECTED_CONVENTION = 'YY2004-quintic-psi-X=1/(1-psi^5)-v3-unscaled-v1'


@lru_cache(maxsize=1)
def module() -> ModuleType:
    """Load a compatible engine without modifying sys.path or legacy imports."""
    name = os.environ.get('BCOV_BACKEND_MODULE', '').strip() or 'bcov.bcov_yy'
    result = importlib.import_module(name)
    for required in ('Poly', 'BCOV', 'YYRing', 'reference_genus2', 'YY_GENUS2_AMBIGUITY'):
        if not hasattr(result, required):
            raise TypeError(f'BCOV backend {name} lacks {required}')
    if result.BCOV.convention != EXPECTED_CONVENTION:
        raise ValueError('Unsupported BCOV convention: '+str(result.BCOV.convention))
    return result


def api() -> dict[str, Any]:
    """Return the required engine symbols for existing analysis callers."""
    m = module()
    return {name: getattr(m, name) for name in (
        'BCOV', 'Poly', 'YYRing', 'reference_genus2', 'YY_GENUS2_AMBIGUITY',
    )}


def validate_small() -> dict[str, str | int | bool]:
    """Check the reference genus-two polynomial and anomaly equation exactly."""
    m = module()
    solver = m.BCOV()
    solver.integrate(2)
    solver.set_ambiguity(2, m.YY_GENUS2_AMBIGUITY)
    if solver.amplitudes[2].polynomial() != m.reference_genus2() or solver.residual(2) != 0:
        raise ValueError('Backend failed the exact genus-two regression')
    return {'backend': m.__name__, 'convention': m.BCOV.convention,
            'genus_two_terms': len(solver.amplitudes[2].polynomial()),
            'bcov_residual_zero': True}
