"""Provenance-aware access to supplied exact YY and empirical PT datasets."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, TypedDict, cast

from .backend import EXPECTED_CONVENTION
from .runner import ROOT, resolve_stage
from .sections import from_reference

PROJECT_ROOT = ROOT


class IntermediateLabel(TypedDict):
    """One exact residual coefficient with explicit reference-derived status."""
    genus: int
    residual_index: int
    coefficient: str
    status: str


def read_reference() -> dict[str, Any]:
    """Read the corrected g2–12 reference, verifying gauge and source checksum."""
    if reference_digest() != provenance()['sha256_reference']:
        raise ValueError('Reference amplitude checksum does not match provenance')
    doc = json.loads((ROOT/'shared/reference_amplitudes.json').read_text())
    if doc.get('convention') != EXPECTED_CONVENTION:
        raise ValueError('Reference and BCOV engine use different gauges')
    return cast(dict[str, Any], doc)


def intermediate_labels() -> list[IntermediateLabel]:
    """Derive the 22 affine-section labels; these are not blind predictions."""
    out: list[IntermediateLabel] = []
    for amplitude in read_reference()['amplitudes']:
        genus = amplitude['genus']
        residual = from_reference(genus, amplitude['ambiguity'])
        out.extend({'genus': genus, 'residual_index': j, 'coefficient': str(value),
                    'status': 'reference_derived_not_prediction'}
                   for j, value in enumerate(residual.coefficients))
    return out


def write_intermediate(path: str | Path) -> list[IntermediateLabel]:
    """Export exact rational labels to a caller-chosen output file."""
    rows = intermediate_labels()
    path = Path(path)
    path.parent.mkdir(exist_ok=True, parents=True)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def read_csv(stage: str, relative_path: str) -> list[dict[str, str]]:
    """Read a packaged study CSV, rejecting paths outside that study."""
    base = ROOT/'studies'/resolve_stage(stage)['id']
    path = (base/relative_path).resolve()
    if not path.is_relative_to(base) or not path.is_file():
        raise ValueError('Invalid study dataset path')
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def provenance() -> dict[str, Any]:
    """Return the original source provenance, assumptions, and seed corrections."""
    return cast(dict[str, Any], json.loads((ROOT/'shared/PROVENANCE.json').read_text()))


def reference_digest() -> str:
    """Return the SHA256 of the bundled corrected reference amplitudes."""
    return hashlib.sha256((ROOT/'shared/reference_amplitudes.json').read_bytes()).hexdigest()
