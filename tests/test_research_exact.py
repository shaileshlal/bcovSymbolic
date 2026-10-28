from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import pytest
from quintic_bcov_lab import backend,datasets,runner
from quintic_bcov_lab.sections import dimensions,kernel,from_reference,difference_from_charges
from quintic_bcov_lab.arithmetic import (difference,newton_forecast,valuation,
    period_four_phase,sparse_pulse_reconstruct,envelope,front,terminal_pulse_depth)


def test_bcov_genus_two():
    r=backend.validate_small()
    assert r['bcov_residual_zero'] and r['genus_two_terms']==14


def test_reference_and_kernel_roundtrips():
    labels=[]
    for item in datasets.read_reference()['amplitudes']:
        g=item['genus'];sp=from_reference(g,item['ambiguity'])
        assert sp.validate()
        _,m,D=dimensions(g)
        assert len(sp.coefficients)==D
        for r in range(D):
            assert len(kernel(g,r))==3*g-2
        labels.extend(sp.coefficients)
    assert len(labels)==22


def test_label_provenance():
    labels=datasets.intermediate_labels()
    assert len(labels)==22
    assert all(z['status']=='reference_derived_not_prediction' for z in labels)
    assert len(set((x['genus'],x['residual_index']) for x in labels))==22


def test_exact_lattice_forward_substitution():
    H=[[Q(2),0,0],[Q(7),Q(-5),0],[Q(3),Q(11),Q(13)]]
    delta=[Q(4),Q(-3),Q(1)]
    n=[sum(a*b for a,b in zip(row,delta)) for row in H]
    assert difference_from_charges(9,n,H)==tuple(delta)
    with pytest.raises(ValueError):difference_from_charges(9,n,[[1,2,0],[0,1,0],[0,0,1]])


def test_newton_pulses_exact():
    initial=[Q(0),Q(1),Q(4)]
    p=[Q(0),Q(0),Q(7),Q(0),Q(-2),Q(0)]
    a=sparse_pulse_reconstruct(initial,p)
    assert difference(a,3)==p
    assert newton_forecast([0,1,4,9],order=2)==[Q(16)]
    assert valuation(5**12*7,5)==12
    assert valuation(0,5) is None
    assert [period_four_phase(k) for k in range(8)]==[1,1,-1,-1,1,1,-1,-1]


def test_normalized_front_formula():
    assert envelope(1,10,8)==5*64+(2*1-5)*8+6-20
    assert front(0,10)==9 and front(3,10)==10
    assert terminal_pulse_depth(1,9)==17
    assert terminal_pulse_depth(3,9)==21


def test_registry_and_digests():
    stages=runner.stages()
    assert len(stages)==18
    assert len({s['id'] for s in stages})==len(stages)
    assert runner.resolve_stage('18')['id']=='18_nested_v12'
    assert datasets.reference_digest()==datasets.provenance()['sha256_reference']


def test_bcov_engine_code_identical_across_legacy():
    root=runner.ROOT
    canonical=(root/'src/quintic_bcov_lab/vendor/engine/bcov_yy.py').read_bytes()
    same=[]
    for p in (root/'studies').rglob('bcov_yy.py'):
        assert p.read_bytes()==canonical, str(p)
        same.append(p)
    assert len(same)>=7
