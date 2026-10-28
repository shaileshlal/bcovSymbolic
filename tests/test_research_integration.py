"""Integration boundaries: backend agreement, immutable provenance, and isolation."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

from bcov import BCOV, YY_GENUS2_AMBIGUITY
from bcov.datasets import iter_amplitudes
from bcov.patterns import pure_x
from quintic_bcov_lab import backend, datasets, runner
from quintic_bcov_lab.arithmetic import fraction, valuation


def test_shared_backend_and_original_engine() -> None:
    """Compare integrated and pinned solver polynomials through formal genus 3."""
    assert backend.module().BCOV is BCOV
    path = runner.ROOT/'src/quintic_bcov_lab/vendor/engine/bcov_yy.py'
    spec = importlib.util.spec_from_file_location('_pinned_bcov_regression', path)
    assert spec and spec.loader
    original = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = original
    spec.loader.exec_module(original)
    live, pinned = BCOV(), original.BCOV()
    for solver in (live, pinned):
        solver.integrate(2)
        solver.set_ambiguity(2, YY_GENUS2_AMBIGUITY)
        solver.integrate(3)
    for g in (2, 3):
        assert live.amplitudes[g].anomaly.terms == pinned.amplitudes[g].anomaly.terms


def test_reference_overlap_all_genera() -> None:
    """All 11 reference ambiguity vectors agree with independently stored HKQ data."""
    references = {a['genus']: tuple(Q(x) for x in a['ambiguity'])
                  for a in datasets.read_reference()['amplitudes']}
    checked = 0
    for amplitude in iter_amplitudes():
        if amplitude.genus in references:
            assert pure_x(amplitude) == references[amplitude.genus]
            checked += 1
    assert checked == 11


def test_archive_integrity() -> None:
    """Every imported file, including original reports, remains byte-identical."""
    assert runner.audit_import() == {'matching': 528, 'changed': 0, 'missing': 0}
    assert runner.audit_sources() == {
        'matching': 480, 'changed_after_reproduction': 0, 'missing': 0,
        'shared': 10, 'optional': 2,
    }


def test_workspace_run_and_source_isolation(tmp_path: Path) -> None:
    """Restore reference aliases and run an original test without changing sources."""
    assert runner.prepare(tmp_path) == 10
    assert runner.prepare(tmp_path) == 0
    result = runner.execute('02', workspace=tmp_path)
    assert result['returncode'] == 0, result['stderr_tail']
    assert (tmp_path/'run_logs/02_high_x_rope_test.json').is_file()
    assert runner.audit_import()['changed'] == 0
    with pytest.raises(ValueError):
        runner.prepare(runner.ROOT)
    unrelated = tmp_path/'unrelated'
    unrelated.mkdir()
    (unrelated/'keep.txt').write_text('user data')
    with pytest.raises(ValueError):
        runner.prepare(unrelated)
    assert (unrelated/'keep.txt').read_text() == 'user data'


def test_invalid_requests_fail() -> None:
    """Reject vacuous runs, ambiguous stages, path traversal and inexact labels."""
    with pytest.raises(KeyError):
        runner.resolve_stage('1')
    with pytest.raises(ValueError):
        runner.execute_all(stages_subset=[])
    with pytest.raises(ValueError):
        datasets.read_csv('18', '../17_depth11_15_v11/data/calibration_deep_Q20_Q40.csv')
    with pytest.raises(TypeError):
        fraction(0.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        valuation(8, 4)


def test_failed_execution_is_logged(tmp_path: Path) -> None:
    """A subprocess launch failure must be visible in the exit status and report."""
    with patch('quintic_bcov_lab.runner.subprocess.run', side_effect=OSError('missing interpreter')):
        result = runner.execute('02', workspace=tmp_path)
    assert result['status'] == 'failed' and result['returncode'] == 127
    assert 'missing interpreter' in result['stderr_tail']
