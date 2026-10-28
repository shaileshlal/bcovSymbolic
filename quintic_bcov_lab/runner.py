"""Reproduce archived studies in writable workspaces, preserving source evidence.

Predictors and validators retain their original separation and working paths.
The packaged archive is read-only input; generated logs and results belong to
an explicitly prepared workspace (by default .bcov-research in the current cwd).
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Literal, TypedDict, cast
import zipfile

ROOT = Path(__file__).resolve().parent/'archive'
CATALOG = ROOT/'shared/stage_catalog.json'
Action = Literal['run', 'test']


class Stage(TypedDict):
    """Chronological stage metadata and original argument vectors."""
    id: str
    archive: str
    title: str
    run: list[str]
    test: list[str]
    version: str


class RunReport(TypedDict):
    """Fresh subprocess evidence, distinct from supplied benchmark snapshots."""
    stage: str
    action: str
    status: str
    returncode: int
    elapsed_s: float
    command: list[str]
    stdout_tail: str
    stderr_tail: str
    original_archive: str
    workspace: str
    comparison: str


def stages() -> list[Stage]:
    """Return the original 18-stage catalog in chronological order."""
    return cast(list[Stage], json.loads(CATALOG.read_text()))


def resolve_stage(identifier: str) -> Stage:
    """Resolve an exact ID or unique prefix; reject empty or ambiguous matches."""
    matches = [s for s in stages() if identifier and s['id'].startswith(identifier)]
    if len(matches) != 1:
        raise KeyError(f'Expected a unique stage for {identifier!r}; found {len(matches)}')
    return matches[0]


def workspace_path(workspace: str | Path | None = None) -> Path:
    """Resolve a writable run directory outside the packaged source archive."""
    path = Path(workspace or Path.cwd()/'.bcov-research').resolve()
    if path == ROOT or path.is_relative_to(ROOT) or ROOT.is_relative_to(path):
        raise ValueError('Workspace must be separate from the packaged archive')
    return path


def prepare(workspace: str | Path | None = None) -> int:
    """Copy inputs once and restore deduplicated references in the run workspace.

    Existing managed workspaces retain their generated outputs. An unrelated
    nonempty directory or a workspace from another archive version is rejected.
    """
    dest = workspace_path(workspace)
    marker = dest/'.bcov-research.json'
    digest = json.loads((ROOT.parent/'IMPORT_MANIFEST.json').read_text())['sha256']
    if marker.exists():
        if json.loads(marker.read_text()).get('archive_sha256') != digest:
            raise ValueError('Workspace belongs to a different source archive; use a new directory')
    else:
        if dest.exists() and any(dest.iterdir()):
            raise ValueError('Use an empty directory or an existing managed research workspace')
        dest.mkdir(parents=True, exist_ok=True)
        for name in ('studies', 'shared'):
            shutil.copytree(ROOT/name, dest/name)
        marker.write_text(json.dumps({'archive_sha256': digest}, indent=2)+'\n')
    source = dest/'shared/reference_amplitudes.json'
    added = 0
    manifest = json.loads((ROOT/'shared/source_files.json').read_text())
    for item in manifest:
        if item['status'] != 'deduplicated:shared/reference_amplitudes.json':
            continue
        target = dest/'studies'/item['stage']/item['path']
        if target.is_file():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            target.symlink_to(os.path.relpath(source, target.parent))
        except (OSError, NotImplementedError):
            shutil.copyfile(source, target)
        added += 1
    return added


def install_large_data(source: str | Path, workspace: str | Path | None = None) -> list[str]:
    """Install the separately supplied exact genus-1000 files into a workspace."""
    source = Path(source)
    names = ('direct_J4_exact.json.gz', 'all_constant_and_leading.json.gz')
    if source.is_dir() and any(not (source/name).is_file() for name in names):
        raise ValueError('Both exact low-X files are required')
    prepare(workspace)
    dest = workspace_path(workspace)/'studies/01_v01_low_x/results'
    installed = []
    for name in names:
        target = dest/('gmp' if name.startswith('direct_') else '')/name
        target.parent.mkdir(parents=True, exist_ok=True)
        if source.is_dir():
            shutil.copy2(source/name, target)
        else:
            with zipfile.ZipFile(source) as archive:
                found = [n for n in archive.namelist() if Path(n).name == name]
                if len(found) != 1:
                    raise ValueError(f'{name} absent or ambiguous')
                with archive.open(found[0]) as inp, target.open('wb') as out:
                    shutil.copyfileobj(inp, out)
        installed.append(str(target))
    return installed


def execute(stage_id: str, mode: Action = 'test', timeout: float = 300,
            verbose: bool = False, workspace: str | Path | None = None) -> RunReport:
    """Run one archived command, recording failure, timeout, and complete logs."""
    if mode not in ('test', 'run') or timeout <= 0:
        raise ValueError('Require action test/run and a positive timeout')
    item = resolve_stage(stage_id)
    prepare(workspace)
    work = workspace_path(workspace)
    cwd = work/'studies'/item['id']
    cmd = [sys.executable, *item[mode]]
    env = os.environ.copy()
    # Stage-local engines win over the live unified backend. Suppress bytecode
    # writes so source audits remain meaningful even for installed wheels.
    paths = [str(cwd), str(cwd/'code'), str(cwd/'quintic_ai_predictor')]
    env['PYTHONPATH'] = os.pathsep.join(paths+[env.get('PYTHONPATH', '')])
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    started = time.monotonic()
    if verbose:
        print('RUN', item['id'], mode, flush=True)
    try:
        result = subprocess.run(cmd, cwd=cwd, env=env, text=True,
                                capture_output=True, timeout=timeout, check=False)
        returncode, stdout, stderr = result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired as error:
        returncode = 124
        stdout = _output(error.stdout)
        stderr = 'Timeout exceeded\n'+_output(error.stderr)
    except OSError as error:
        returncode, stdout, stderr = 127, '', str(error)
    report: RunReport = {
        'stage': item['id'], 'action': mode, 'status': 'passed' if returncode == 0 else 'failed',
        'returncode': returncode, 'elapsed_s': round(time.monotonic()-started, 3),
        'command': cmd, 'stdout_tail': stdout[-2500:], 'stderr_tail': stderr[-2500:],
        'original_archive': item['archive'], 'workspace': str(work),
        'comparison': 'Fresh execution of archived study; retain its predictor/validator assumptions',
    }
    logs = work/'run_logs'
    logs.mkdir(exist_ok=True)
    (logs/f'{item["id"]}_{mode}.log').write_text(stdout+'\n--- STDERR ---\n'+stderr)
    (logs/f'{item["id"]}_{mode}.json').write_text(json.dumps(report, indent=2)+'\n')
    if verbose:
        print(report['status'], item['id'], report['elapsed_s'], 'seconds', flush=True)
    return report


def _output(value: str | bytes | None) -> str:
    """Normalize captured timeout output across Python versions."""
    return value.decode(errors='replace') if isinstance(value, bytes) else value or ''


def execute_all(mode: Action = 'test', limit_minutes: float = 10,
                stages_subset: list[str] | None = None, continue_on_error: bool = True,
                verbose: bool = True, workspace: str | Path | None = None) -> list[RunReport]:
    """Run selected studies chronologically; limit_minutes is per-stage timeout."""
    selected = {resolve_stage(s)['id'] for s in stages_subset} if stages_subset is not None else {s['id'] for s in stages()}
    if not selected or limit_minutes <= 0:
        raise ValueError('Select at least one stage and a positive per-stage timeout')
    reports = []
    for item in stages():
        if item['id'] not in selected:
            continue
        report = execute(item['id'], mode, 60*limit_minutes, verbose, workspace)
        reports.append(report)
        if not continue_on_error and report['returncode'] != 0:
            break
    path = workspace_path(workspace)/'run_logs'/f'latest_{mode}_all.json'
    path.write_text(json.dumps(reports, indent=2)+'\n')
    return reports


def audit_sources(workspace: str | Path | None = None) -> dict[str, int]:
    """Audit original stage hashes, or report changes in an existing workspace."""
    root = ROOT if workspace is None else workspace_path(workspace)
    counts = {'matching': 0, 'changed_after_reproduction': 0, 'missing': 0, 'shared': 0, 'optional': 0}
    for record in json.loads((ROOT/'shared/source_files.json').read_text()):
        if record['status'].startswith('deduplicated'):
            counts['shared'] += 1
        elif record['status'].startswith('optional:'):
            counts['optional'] += 1
        elif record['status'] == 'bundled':
            path = root/'studies'/record['stage']/record['path']
            if not path.is_file():
                counts['missing'] += 1
            elif hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
                counts['changed_after_reproduction'] += 1
            else:
                counts['matching'] += 1
    return counts


def audit_import() -> dict[str, int]:
    """Check every imported file against the SHA256 recorded at ZIP extraction."""
    expected = json.loads((ROOT.parent/'IMPORT_MANIFEST.json').read_text())['files']
    counts = {'matching': 0, 'changed': 0, 'missing': 0}
    for relative, digest in expected.items():
        path = ROOT/relative
        key = 'missing' if not path.is_file() else 'matching' if hashlib.sha256(path.read_bytes()).hexdigest() == digest else 'changed'
        counts[key] += 1
    return counts
