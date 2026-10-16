"""Exact quintic BCOV recursion in the original Yamaguchi--Yau convention.

Python >= 3.10; standard library only. Run ``python bcov_yy.py --genus 3``.
The polynomial algebra is commutative. theta is a derivation, not a variable.
See README.md for normalization, equations, boundary data, and limitations.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import lcm
import json
import os
from pathlib import Path
import tempfile
from time import perf_counter
from types import MappingProxyType
from typing import Iterable, Mapping

NAMES = ("u", "v1", "v2", "v3", "X")
WEIGHTS = (1, 1, 2, 3, 1)
ZERO_MONOMIAL = (0, 0, 0, 0, 0)
Monomial = tuple[int, int, int, int, int]

# Opt-in resource limits used by the physical driver. Zero means unlimited.
MAX_TERMS = 0
MAX_PRODUCTS = 0


class ResourceLimitError(RuntimeError):
    pass


def set_resource_limits(max_terms=0, max_products=0):
    global MAX_TERMS, MAX_PRODUCTS
    if max_terms < 0 or max_products < 0:
        raise ValueError("Resource limits must be nonnegative")
    MAX_TERMS, MAX_PRODUCTS = max_terms, max_products



def rational(value: int | Fraction) -> Fraction:
    """Forbid implicit conversion of floating point input into exact data."""
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    raise TypeError("Exact coefficients must be int or Fraction; floats are forbidden")


def variable_index(variable: str | int) -> int:
    if isinstance(variable, str):
        return NAMES.index(variable)
    if type(variable) is int and 0 <= variable < len(NAMES):
        return variable
    raise ValueError(f"Unknown variable: {variable!r}")


class Poly:
    """Canonical sparse QQ[u,v1,v2,v3,X], immutable through the public API.

    Only occupied exponent tuples are stored; coefficients are reduced Fractions.
    No automatic pretty printing, symbolic simplification, or degree truncation.
    """
    __slots__ = ("_terms", "_hash", "_integer_form")

    def __init__(self, terms: Mapping[Monomial, int | Fraction]
                 | Mapping[tuple[int, ...], int | Fraction] | None = None):
        data = {}
        for powers, coefficient in (terms or {}).items():
            if len(powers) != 5 or any(type(p) is not int or p < 0 for p in powers):
                raise ValueError("Monomials require five nonnegative integer exponents")
            coefficient = rational(coefficient)
            if coefficient:
                data[tuple(powers)] = coefficient
        if MAX_TERMS and len(data) > MAX_TERMS:
            raise ResourceLimitError(f"Polynomial terms {len(data)} exceed budget {MAX_TERMS}")
        self._terms = MappingProxyType(data)
        self._hash = None
        self._integer_form = None

    @classmethod
    def constant(cls, value: int | Fraction) -> Poly:
        return cls({ZERO_MONOMIAL: rational(value)})

    @classmethod
    def variable(cls, variable: str | int) -> Poly:
        powers = [0] * 5
        powers[variable_index(variable)] = 1
        return cls({tuple(powers): 1})

    @classmethod
    def coerce(cls, value: Poly | int | Fraction) -> Poly:
        return value if isinstance(value, Poly) else cls.constant(value)

    @property
    def terms(self) -> Mapping[Monomial, Fraction]:
        return self._terms

    def __bool__(self):
        return bool(self._terms)

    def __len__(self):
        return len(self._terms)

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(frozenset(self._terms.items()))
        return self._hash

    def __eq__(self, other):
        if isinstance(other, Poly):
            return self._terms == other._terms
        if type(other) is int or isinstance(other, Fraction):
            return self == Poly.constant(other)
        return NotImplemented

    def __add__(self, other):
        other = Poly.coerce(other)
        data = dict(self._terms)
        for m, c in other._terms.items():
            data[m] = data.get(m, 0) + c
            if not data[m]:
                del data[m]
        return Poly(data)

    __radd__ = __add__

    def __neg__(self):
        return self.scale(-1)

    def __sub__(self, other):
        return self + (-Poly.coerce(other))

    def __rsub__(self, other):
        return Poly.coerce(other) - self

    def scale(self, coefficient: int | Fraction) -> Poly:
        c = rational(coefficient)
        if not c or not self:
            return Poly()
        if c == 1:
            return self
        return Poly({m: c * a for m, a in self._terms.items()})

    def integer_form(self):
        """Cached common denominator; inner products then use integer arithmetic."""
        if self._integer_form is None:
            denominator = lcm(*(c.denominator for c in self._terms.values()))
            integers = {m:c.numerator*(denominator//c.denominator)
                        for m,c in self._terms.items()}
            self._integer_form = (MappingProxyType(integers), denominator)
        return self._integer_form

    def __mul__(self, other):
        if not isinstance(other, Poly):
            return self.scale(other)
        if not self or not other:
            return Poly()
        if len(self) == 1 and ZERO_MONOMIAL in self._terms:
            return other.scale(self._terms[ZERO_MONOMIAL])
        if len(other) == 1 and ZERO_MONOMIAL in other._terms:
            return self.scale(other._terms[ZERO_MONOMIAL])
        if MAX_PRODUCTS and len(self)*len(other) > MAX_PRODUCTS:
            raise ResourceLimitError(f"Pair products {len(self)*len(other)} exceed budget {MAX_PRODUCTS}")
        left, ld = self.integer_form()
        right, rd = other.integer_form()
        data = {}
        for a, ca in left.items():
            for b, cb in right.items():
                m = tuple(i + j for i, j in zip(a, b))
                data[m] = data.get(m, 0) + ca * cb
            if MAX_TERMS and len(data) > MAX_TERMS:
                raise ResourceLimitError(f"Intermediate terms exceed budget {MAX_TERMS}")
        denominator = ld*rd
        return Poly({m:Fraction(c,denominator) for m,c in data.items() if c})

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = rational(scalar)
        if not scalar:
            raise ZeroDivisionError("Division by zero")
        return self.scale(1 / scalar)

    def __pow__(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("Only nonnegative integer powers are polynomial")
        result, base = Poly.constant(1), self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent >>= 1
            if exponent:
                base = base * base
        return result

    def partial(self, variable: str | int) -> Poly:
        i = variable_index(variable)
        data = {}
        for m, c in self._terms.items():
            if m[i]:
                p = list(m)
                p[i] -= 1
                data[tuple(p)] = c * m[i]
        return Poly(data)

    def integrate(self, variable: str | int) -> Poly:
        """Primitive with zero integration constant (not a boundary condition)."""
        i = variable_index(variable)
        data = {}
        for m, c in self._terms.items():
            p = list(m)
            p[i] += 1
            data[tuple(p)] = c / p[i]
        return Poly(data)

    def coefficient(self, variable: str | int, power: int) -> Poly:
        """Coefficient polynomial with that variable removed from its exponents."""
        if type(power) is not int or power < 0:
            raise ValueError("Coefficient power must be a nonnegative integer")
        i = variable_index(variable)
        data = {}
        for m, c in self._terms.items():
            if m[i] == power:
                p = list(m)
                p[i] = 0
                data[tuple(p)] = c
        return Poly(data)

    def degree(self, variable: str | int | None = None) -> int:
        """Weighted degree when variable=None; degree(0 polynomial) is -1."""
        if variable is None:
            return max((sum(a * b for a, b in zip(m, WEIGHTS)) for m in self._terms), default=-1)
        i = variable_index(variable)
        return max((m[i] for m in self._terms), default=-1)

    def homogeneous_part(self, weight: int) -> Poly:
        return Poly({m: c for m, c in self._terms.items()
                     if sum(a * b for a, b in zip(m, WEIGHTS)) == weight})

    def substitute(self, replacements: Mapping[str | int, Poly | int | Fraction]) -> Poly:
        """Simultaneous substitution; replacements do not act on one another."""
        replacements = {variable_index(k): Poly.coerce(v) for k, v in replacements.items()}
        cache = {}
        result = Poly()
        for m, c in self._terms.items():
            term = Poly.constant(c)
            for i, exponent in enumerate(m):
                if exponent:
                    key = (i, exponent)
                    if key not in cache:
                        cache[key] = replacements.get(i, GENERATORS[i]) ** exponent
                    term = term * cache[key]
            result = result + term
        return result

    def evaluate(self, values: Mapping[str, int | Fraction]) -> Fraction:
        required = {NAMES[i] for m in self._terms for i, e in enumerate(m) if e}
        if not required.issubset(values):
            raise ValueError(f"Missing values for {sorted(required - values.keys())}")
        vals = {variable_index(k): rational(v) for k, v in values.items()}
        return sum((c * _product(vals[i] ** e for i, e in enumerate(m) if e)
                    for m, c in self._terms.items()), Fraction(0))

    def to_data(self) -> list:
        return [[list(m), str(c.numerator), str(c.denominator)]
                for m, c in sorted(self._terms.items())]

    @classmethod
    def from_data(cls, rows: list) -> Poly:
        data = {}
        for m, numerator, denominator in rows:
            key = tuple(m)
            if key in data:
                raise ValueError("Duplicate monomial in checkpoint")
            data[key] = Fraction(int(numerator), int(denominator))
        return cls(data)

    def __repr__(self):
        return f"Poly(terms={len(self)}, weighted_degree={self.degree()})"

    def __str__(self):
        pieces = []
        order = lambda item: (sum(a*b for a,b in zip(item[0], WEIGHTS)), item[0])
        for m, c in sorted(self._terms.items(), key=order):
            factors = [name if power == 1 else f"{name}^{power}"
                       for name, power in zip(NAMES, m) if power]
            body = "*".join(factors)
            magnitude = abs(c)
            term = body if body and magnitude == 1 else str(magnitude) + ("*" + body if body else "")
            pieces.append((" - " if c < 0 else " + ") + term)
        if not pieces:
            return "0"
        result = "".join(pieces)
        return ("-" + result[3:]) if result.startswith(" - ") else result[3:]


def _product(values: Iterable[Fraction]) -> Fraction:
    result = Fraction(1)
    for value in values:
        result *= value
    return result


GENERATORS = tuple(Poly.variable(i) for i in range(5))
u, v1, v2, v3, X = GENERATORS


class YYRing:
    """Quintic theta=psi*d/dpsi and normalized covariant differentiation.

    The table is derived from the AB relations in YY (3.16),(3.22),(3.29),
    so the complete u dependence is retained before anomaly coefficient extraction.
    """
    def __init__(self, cache_size: int = 8192):
        A = v1 - 1 - 2*u
        B2 = v2 + u*v1
        B3 = v3 - u*v2 + 5*u*X*(v1-2)
        du = B2 - u**2
        dX = 5*X*(X-1)
        A2 = -4*B2 - 2*A*u - 2*u + 2*u**2 - 2*A + 10*X*u + 5*X*A - 5*X - 1
        dv1 = A2 - A**2 + 2*du
        dv2 = B3 - u*B2 - du*v1 - u*dv1
        dB3 = 10*X*B3 - 35*X*B2 + 50*X*u - 24*X - u*B3
        dv3 = dB3 + du*v2 + u*dv2 - 5*(du*X*(v1-2) + u*dX*(v1-2) + u*X*dv1)
        self.table = (du, dv1, dv2, dv3, dX)
        # Instance-local, bounded cache; no process-global polynomial retention.
        self._transition = lru_cache(maxsize=cache_size)(self._monomial_derivative)

    def _monomial_derivative(self, m: Monomial) -> Poly:
        data = {}
        for i, exponent in enumerate(m):
            if not exponent:
                continue
            base = list(m)
            base[i] -= 1
            for shift, coefficient in self.table[i].terms.items():
                out = tuple(a + b for a, b in zip(base, shift))
                data[out] = data.get(out, 0) + exponent * coefficient
        return Poly(data)

    def theta(self, p: Poly | int | Fraction) -> Poly:
        data = {}
        for m, coefficient in Poly.coerce(p).terms.items():
            for out, value in self._transition(m).terms.items():
                data[out] = data.get(out, 0) + coefficient * value
        return Poly(data)

    def theta_z(self, p: Poly | int | Fraction) -> Poly:
        """z*d/dz for z=(5*psi)^(-5); this is not d/dz."""
        return self.theta(p) / -5

    def covariant(self, p: Poly, genus: int, insertions: int) -> Poly:
        """P_g^(n) -> P_g^(n+1); genus and tensor rank are external metadata."""
        if type(genus) is not int or genus < 0 or type(insertions) is not int or insertions < 0:
            raise ValueError("Genus and insertion count must be nonnegative integers")
        connection = -insertions*(v1-2*u) + (2*genus-2)*(u-Fraction(5,2)*X)
        return self.theta(p) + connection*p

    def clear_cache(self):
        self._transition.cache_clear()


class Recipe:
    """An interned arithmetic DAG handle; expansion is explicit and memoized."""
    __slots__ = ("arena", "node")

    def __init__(self, arena: RecipeArena, node: int):
        self.arena, self.node = arena, node

    def __add__(self, other):
        return self.arena.binary("add", self, other)

    __radd__ = __add__

    def __mul__(self, other):
        return self.arena.binary("mul", self, other)

    __rmul__ = __mul__

    def __neg__(self):
        return self * -1

    def expand(self) -> Poly:
        return self.arena.expand(self)


class RecipeArena:
    """Keep repeated amplitude leaves, products and derivative vertices atomic.

    A recipe DAG delays expansion; it does not eliminate the expanded polynomial
    required by the current gradient solver. clear_values() releases expansions.
    """
    def __init__(self, ring: YYRing):
        self.ring = ring
        self.nodes: list[tuple] = []
        self._intern: dict[tuple, int] = {}
        self._values: dict[int, Poly] = {}

    def _make(self, spec: tuple) -> Recipe:
        if spec not in self._intern:
            self._intern[spec] = len(self.nodes)
            self.nodes.append(spec)
        return Recipe(self, self._intern[spec])

    def leaf(self, p: Poly | int | Fraction) -> Recipe:
        return self._make(("leaf", Poly.coerce(p)))

    def coerce(self, p: Recipe | Poly | int | Fraction) -> Recipe:
        if isinstance(p, Recipe):
            if p.arena is not self:
                raise ValueError("Recipes from different arenas cannot be mixed")
            return p
        return self.leaf(p)

    def binary(self, op: str, a, b) -> Recipe:
        if op not in ("add", "mul"):
            raise ValueError("Unknown recipe operation")
        a, b = self.coerce(a), self.coerce(b)
        left, right = sorted((a.node, b.node))
        return self._make((op, left, right))

    def covariant(self, p, genus: int, insertions: int) -> Recipe:
        return self._make(("cov", self.coerce(p).node, genus, insertions))

    def expand(self, recipe: Recipe) -> Poly:
        recipe = self.coerce(recipe)
        # Topological order avoids Python recursion limits for long arithmetic DAGs.
        needed, stack = set(), [recipe.node]
        while stack:
            node = stack.pop()
            if node in needed or node in self._values:
                continue
            needed.add(node)
            spec = self.nodes[node]
            if spec[0] in ("add", "mul"):
                stack.extend(spec[1:])
            elif spec[0] == "cov":
                stack.append(spec[1])
        for node in sorted(needed):
            spec = self.nodes[node]
            if spec[0] == "leaf":
                result = spec[1]
            elif spec[0] == "add":
                result = self._values[spec[1]] + self._values[spec[2]]
            elif spec[0] == "mul":
                result = self._values[spec[1]] * self._values[spec[2]]
            else:
                result = self.ring.covariant(self._values[spec[1]], spec[2], spec[3])
            self._values[node] = result
        return self._values[recipe.node]

    def clear_values(self):
        self._values.clear()

    def stats(self) -> dict:
        return {"nodes": len(self.nodes), "expanded_nodes": len(self._values),
                "operations": dict(Counter(spec[0] for spec in self.nodes))}


class Amplitude:
    """Anomaly representative plus a separately stored rational ambiguity vector.

    Vector slot j multiplies X^j. None means unresolved, not physically zero.
    The anomaly representative vanishes at v1=v2=v3=0 by integration convention.
    """
    __slots__ = ("genus", "anomaly", "ambiguity")

    def __init__(self, genus: int, anomaly: Poly, ambiguity=None):
        if type(genus) is not int or genus < 2:
            raise ValueError("An amplitude requires integer genus >= 2")
        anomaly = Poly.coerce(anomaly)
        if anomaly.degree("u") > 0 or anomaly.degree() > 3*genus-3:
            raise ValueError("Amplitude violates u independence or YY degree bound")
        if anomaly.substitute({"v1":0,"v2":0,"v3":0}):
            raise ValueError("Store pure-X terms in the separate ambiguity vector")
        vector = None if ambiguity is None else tuple(rational(c) for c in ambiguity)
        if vector is not None and len(vector) > 3*genus-2:
            raise ValueError("Ambiguity exceeds degree 3g-3")
        self.genus, self.anomaly, self.ambiguity = genus, anomaly, vector

    def with_ambiguity(self, coefficients: Iterable[int | Fraction]) -> Amplitude:
        return Amplitude(self.genus, self.anomaly, coefficients)

    def polynomial(self) -> Poly:
        if self.ambiguity is None:
            raise ValueError(f"Genus {self.genus} ambiguity is unresolved; choose explicit coefficients")
        result = self.anomaly
        for power, coefficient in enumerate(self.ambiguity):
            if coefficient:
                result = result + coefficient*X**power
        return result

    def to_data(self) -> dict:
        return {"genus": self.genus, "anomaly": self.anomaly.to_data(),
                "ambiguity": None if self.ambiguity is None else [str(c) for c in self.ambiguity]}

    @classmethod
    def from_data(cls, data: dict) -> Amplitude:
        vector = data["ambiguity"]
        return cls(data["genus"], Poly.from_data(data["anomaly"]),
                   None if vector is None else tuple(Fraction(c) for c in vector))


def integrate_gradient(d1: Poly, d2: Poly, d3: Poly) -> Poly:
    """Solve the complete three-variable anomaly, keeping the X constant free."""
    gradient = (d1, d2, d3)
    variables = ("v1", "v2", "v3")
    if any(p.degree("u") > 0 for p in gradient):
        raise ValueError("Anomaly gradient must be independent of u")
    for i in range(3):
        for j in range(i):
            if gradient[i].partial(variables[j]) != gradient[j].partial(variables[i]):
                raise ValueError(f"Nonintegrable anomaly gradient: {variables[i]}, {variables[j]}")
    result = Poly()
    for variable, target in reversed(list(zip(variables, gradient))):
        residual = target - result.partial(variable)
        result = result + residual.integrate(variable)
    if any(result.partial(variable) != target for variable, target in zip(variables, gradient)):
        raise ValueError("Gradient integration failed")
    return result


class BCOV:
    """Exact YY recursion, with explicit ambiguity choices between genera."""
    convention = "YY2004-quintic-psi-X=1/(1-psi^5)-v3-unscaled-v1"

    def __init__(self):
        self.ring = YYRing()
        self.arena = RecipeArena(self.ring)
        self.amplitudes: dict[int, Amplitude] = {}
        self._vertices: dict[tuple[int,int], Recipe] = {}
        self._sources: dict[int, Recipe] = {}

    def vertex(self, genus: int, insertions: int) -> Recipe:
        if type(genus) is not int or type(insertions) is not int:
            raise ValueError("Genus and insertions must be integers")
        if genus < 0 or insertions < (3 if genus == 0 else 1 if genus == 1 else 0):
            raise ValueError("Unavailable unstable genus/insertion pair")
        key = (genus, insertions)
        if key in self._vertices:
            return self._vertices[key]
        if genus == 0 and insertions == 3:
            result = self.arena.leaf(1)
        elif genus == 1 and insertions == 1:
            # YY (3.8); genus-one *one-point* seed, not F_1 itself.
            result = self.arena.leaf(Fraction(125,12)-v1/2-Fraction(28,3)*u+Fraction(5,12)*X)
        elif genus >= 2 and insertions == 0:
            if genus not in self.amplitudes:
                raise ValueError(f"Genus {genus} has not been integrated")
            result = self.arena.leaf(self.amplitudes[genus].polynomial())
        else:
            result = self.arena.covariant(self.vertex(genus, insertions-1), genus, insertions-1)
        self._vertices[key] = result
        return result

    def source_recipe(self, genus: int) -> Recipe:
        """R_g=-1/2[P_(g-1)^(2)+sum P_h^(1) P_(g-h)^(1)]."""
        if type(genus) is not int or genus < 2:
            raise ValueError("Recursion starts at genus 2")
        if genus not in self._sources:
            source = self.vertex(genus-1, 2)
            for h in range(1, genus//2+1):
                coefficient = 1 if h == genus-h else 2
                source = source + coefficient*self.vertex(h,1)*self.vertex(genus-h,1)
            self._sources[genus] = Fraction(-1,2)*source
        return self._sources[genus]

    def anomaly_gradient(self, genus: int) -> tuple[Poly, Poly, Poly]:
        source = self.source_recipe(genus).expand()
        if source.degree("u") > 2:
            raise ValueError("Expected a quadratic u source; check lower-genus conventions")
        r0, r1, r2 = (source.coefficient("u", k) for k in range(3))
        # (d_v1 - u*d_v2 - u*(u+5X)*d_v3) P_g = R_g.
        return r0, -r1 + 5*X*r2, -r2

    def integrate(self, genus: int) -> Amplitude:
        if genus in self.amplitudes:
            return self.amplitudes[genus]
        anomaly = integrate_gradient(*self.anomaly_gradient(genus))
        amplitude = Amplitude(genus, anomaly)
        self.amplitudes[genus] = amplitude
        return amplitude

    def set_ambiguity(self, genus: int, coefficients: Iterable[int | Fraction]) -> Amplitude:
        if any(g > genus for g in self.amplitudes):
            raise ValueError("Changing lower boundary data would invalidate computed higher genera; use a new solver")
        amplitude = self.integrate(genus).with_ambiguity(coefficients)
        self.amplitudes[genus] = amplitude
        # Derivatives and sources containing the old choice must be rebuilt.
        self._vertices = {k:v for k,v in self._vertices.items() if k[0] < genus}
        self._sources = {g:v for g,v in self._sources.items() if g <= genus}
        self.arena.clear_values()
        return amplitude

    def residual(self, genus: int, polynomial: Poly | None = None) -> Poly:
        p = self.amplitudes[genus].polynomial() if polynomial is None else polynomial
        return p.partial("v1")-u*p.partial("v2")-u*(u+5*X)*p.partial("v3")-self.source_recipe(genus).expand()

    def compact_recipes(self):
        """Release temporary DAGs/expansions; keep exact amplitude checkpoints.

        Handles previously returned by source_recipe remain valid snapshots if
        the caller still holds them, so release those handles to reclaim memory.
        """
        self.arena = RecipeArena(self.ring)
        self._vertices.clear()
        self._sources.clear()
        self.ring.clear_cache()

    def save(self, path: str | Path):
        """Versioned exact JSON, written by atomic replacement; no pickle."""
        data = {"schema":1, "convention":self.convention, "amplitudes":
                [a.to_data() for _,a in sorted(self.amplitudes.items())]}
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=path.name+".", dir=path.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as stream:
                json.dump(data, stream, indent=2)
                stream.write("\n")
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

    @classmethod
    def load(cls, path: str | Path) -> BCOV:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data["schema"] != 1 or data["convention"] != cls.convention:
            raise ValueError("Unsupported checkpoint schema/convention")
        result = cls()
        expected = 2
        for item in data["amplitudes"]:
            a = Amplitude.from_data(item)
            if a.genus != expected:
                raise ValueError("Checkpoint genera must be consecutive from genus 2")
            # Validate mathematics as well as structure before accepting a checkpoint.
            computed = result.integrate(a.genus)
            if computed.anomaly != a.anomaly:
                raise ValueError(f"Checkpoint fails the anomaly equation at genus {a.genus}")
            result.amplitudes[a.genus] = a
            expected += 1
        return result


# Exact pure-X part of YY (4.1), in this file's original psi convention.
# These are supplied reference boundary data, not inferred from the anomaly.
YY_GENUS2_AMBIGUITY = (Fraction(3125,144), Fraction(-28795,144),
                       Fraction(205,144), Fraction(25,48))


def reference_genus2() -> Poly:
    """Independent published regression target, YY (4.1)."""
    return (Fraction(3125,144)-Fraction(15625,288)*v1+Fraction(125,24)*v1**2
            -Fraction(5,24)*v1**3-Fraction(3125,36)*v2+Fraction(25,6)*v1*v2
            +Fraction(350,9)*v3-Fraction(28795,144)*X-Fraction(835,144)*v1*X
            +Fraction(5,6)*v1**2*X-Fraction(2375,12)*v2*X+Fraction(205,144)*X**2
            -Fraction(325,288)*v1*X**2+Fraction(25,48)*X**3)


def demo(max_genus: int = 3, checkpoint: str | None = None) -> BCOV:
    if max_genus < 2:
        raise ValueError("Choose genus >= 2")
    solver = BCOV()
    print("Exact quintic BCOV recursion; original YY psi convention")
    print("P_1^(1) =", solver.vertex(1,1).expand())
    for genus in range(2, max_genus+1):
        start = perf_counter()
        a = solver.integrate(genus)
        print(f"g={genus}: anomaly terms={len(a.anomaly)}, weight<={a.anomaly.degree()}, "
              f"u-independent={a.anomaly.partial('u') == 0}")
        if genus == 2:
            a = solver.set_ambiguity(genus, YY_GENUS2_AMBIGUITY)
            if a.polynomial() != reference_genus2():
                raise AssertionError("Published genus-two regression failed")
            print("Published genus-two polynomial reproduced exactly.")
            print("P_2 =", a.polynomial())
        else:
            print(f"g={genus}: ambiguity unresolved; {3*genus-2} rational slots available.")
            if genus < max_genus:
                solver.set_ambiguity(genus, ())
                print(f"g={genus}: explicitly choosing zero ambiguity for a formal continuation.")
        p = solver.amplitudes[genus].anomaly
        if solver.residual(genus, p) != 0:
            raise AssertionError("Nonzero BCOV residual")
        print(f"BCOV residual=0, elapsed={perf_counter()-start:.3f}s")
    print("Shared recipe DAG:", solver.arena.stats())
    if checkpoint:
        solver.save(checkpoint)
        print("Exact checkpoint:", checkpoint)
    return solver


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, default=3)
    parser.add_argument("--checkpoint")
    args = parser.parse_args()
    demo(args.genus, args.checkpoint)
