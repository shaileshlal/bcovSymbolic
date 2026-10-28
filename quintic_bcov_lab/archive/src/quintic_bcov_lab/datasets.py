"""Central provenance-aware access to exact YY and empirical PT datasets."""
from __future__ import annotations
import json,csv,hashlib
from pathlib import Path
from fractions import Fraction as Q
from .sections import from_reference

PROJECT_ROOT=Path(__file__).resolve().parents[2]


def read_reference():
    doc=json.loads((PROJECT_ROOT/'shared/reference_amplitudes.json').read_text())
    expected='YY2004-quintic-psi-X=1/(1-psi^5)-v3-unscaled-v1'
    if doc.get('convention')!=expected:
        raise ValueError('Reference and BCOV engine use different gauges')
    return doc


def intermediate_labels():
    out=[]
    for amp in read_reference()['amplitudes']:
        g=amp['genus'];r=from_reference(g,[Q(s) for s in amp['ambiguity']])
        out.extend({'genus':g,'residual_index':j,'coefficient':str(val),
                    'status':'reference_derived_not_prediction'}
                   for j,val in enumerate(r.coefficients))
    return out


def write_intermediate(path):
    rows=intermediate_labels();path=Path(path);path.parent.mkdir(exist_ok=True,parents=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows


def read_csv(stage,relative_path):
    if '..' in Path(relative_path).parts:raise ValueError('invalid relative dataset path')
    p=PROJECT_ROOT/'studies'/stage/relative_path
    with p.open(newline='') as stream:return list(csv.DictReader(stream))


def provenance():return json.loads((PROJECT_ROOT/'shared/PROVENANCE.json').read_text())


def reference_digest():return hashlib.sha256((PROJECT_ROOT/'shared/reference_amplitudes.json').read_bytes()).hexdigest()
