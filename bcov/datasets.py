"""Read packaged source data without depending on the working directory.

Atomic rows are HKQ amplitudes in the factored basis documented in the
convention ledger. Unknown GV entries are never interpreted as zero here.
"""
from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from fractions import Fraction
import gzip
from importlib.resources import files
import json
from pathlib import Path

from .bcov_yy import Poly, X, v1, v2, v3
from .quintic_frames import Frame, FramedPolynomial

HKQ_BASIS = ("1-X", "X", "v1-3/5", "v2+2/25", "v3-X*(v2+4/125)")


@dataclass(frozen=True)
class AtomicTerm:
    """Coefficient times (1-X)^r X^n alpha^j beta^k gamma^ell in HKQ."""

    r: int
    n: int
    j: int
    k: int
    ell: int
    coefficient: Fraction

    def __post_init__(self) -> None:
        """Reject negative or noninteger exponents and inexact coefficients."""
        if any(type(x) is not int or x < 0 for x in (self.r, self.n, self.j, self.k, self.ell)):
            raise ValueError("Atomic exponents must be nonnegative integers")
        if not isinstance(self.coefficient, (int, Fraction)):
            raise TypeError("Atomic coefficients must be exact rationals")
        object.__setattr__(self, "coefficient", Fraction(self.coefficient))


@dataclass(frozen=True)
class AtomicAmplitude:
    """One source amplitude, with its genus and immutable HKQ atomic terms."""

    genus: int
    terms: tuple[AtomicTerm, ...]

    def __post_init__(self) -> None:
        """Validate the genus carried by the source amplitude."""
        if type(self.genus) is not int or self.genus < 2:
            raise ValueError("Atomic amplitudes require an integer genus >= 2")

    def expand(self, frame: Frame = "hkq") -> FramedPolynomial:
        """Expand the full polynomial; intended for small genera only."""
        result = Poly()
        alpha, beta = v1 - Fraction(3, 5), v2 + Fraction(2, 25)
        gamma = v3 - X * (v2 + Fraction(4, 125))
        for t in self.terms:
            result += t.coefficient * (1-X)**t.r * X**t.n * alpha**t.j * beta**t.k * gamma**t.ell
        return FramedPolynomial(result, self.genus, "hkq").to(frame)


def iter_amplitudes(source: str | Path | None = None) -> Iterator[AtomicAmplitude]:
    """Stream source genera; default to the compressed dataset inside the wheel."""
    resource = files("bcov").joinpath("data/bcov_hkq_atomic_genus_2_15.jsonl.gz")
    with (Path(source).open("rb") if source is not None else resource.open("rb")) as raw:
        with gzip.open(raw, "rt", encoding="utf8") as stream:
            seen: set[int] = set()
            for line in stream:
                obj = json.loads(line)
                if tuple(obj.get("basis", ())) != HKQ_BASIS:
                    raise ValueError("Atomic data must declare the expected HKQ basis")
                amplitude = AtomicAmplitude(obj["genus"], tuple(
                    AtomicTerm(r, n, j, k, ell, Fraction(int(num), int(den)))
                    for r, n, j, k, ell, num, den in obj["terms"]
                ))
                if amplitude.genus in seen:
                    raise ValueError(f"Duplicate atomic genus {amplitude.genus}")
                seen.add(amplitude.genus)
                yield amplitude


def load_amplitude(genus: int, source: str | Path | None = None) -> AtomicAmplitude:
    """Load a supplied genus, failing explicitly if no source amplitude exists."""
    if type(genus) is not int or genus < 2:
        raise ValueError("Choose an integer genus >= 2")
    for amplitude in iter_amplitudes(source):
        if amplitude.genus == genus:
            return amplitude
    raise KeyError(f"No supplied atomic amplitude at genus {genus}")
