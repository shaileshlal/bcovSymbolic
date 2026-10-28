"""Exact quintic BCOV toolkit; all coefficients use rational arithmetic."""
from .bcov_yy import (
    Amplitude, BCOV, Poly, Recipe, RecipeArena, YYRing,
    YY_GENUS2_AMBIGUITY, X, integrate_gradient, reference_genus2, u, v1, v2, v3,
)
from .quintic_frames import Frame, FramedPolynomial

__version__ = "0.3.0"
__all__ = [
    "Amplitude", "BCOV", "Poly", "Recipe", "RecipeArena", "YYRing",
    "YY_GENUS2_AMBIGUITY", "X", "integrate_gradient", "reference_genus2",
    "u", "v1", "v2", "v3", "Frame", "FramedPolynomial",
]
