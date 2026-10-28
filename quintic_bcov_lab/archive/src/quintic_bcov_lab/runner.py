"""Sequential reproduction of the 18 original studies, with deterministic provenance.

Runs the ORIGINAL archived study scripts in their directory, without moving
withheld reference files into predictor inputs. Individual historical tests
implement their own isolated/no-oracle controls.
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

ROOT=Path(__file__).resolve().parents[2]
CATALOG=ROOT/'shared/stage_catalog.json'


def stages():return json.loads(CATALOG.read_text())


def resolve_stage(identifier):
    identifiers=[x for x in stages() if x['id']==identifier or x['id'].startswith(identifier+'_') or x['id'].startswith(identifier)]
    if len(identifiers)!=1:raise KeyError(f'Expected unique stage, got {len(identifiers)} for {identifier!r}: {[x["id"] for x in identifiers]}')
    return identifiers[0]


def prepare():
    """Re-create deduplicated shared reference aliases omitted from ZIP.

    Symlinks used on Unix, copies on platforms without symlink support.
    No historical output, calibration or hidden target is changed here.
    """
    source=ROOT/'shared/reference_amplitudes.json'
    if not source.is_file():raise FileNotFoundError(source)
    manifest=json.loads((ROOT/'shared/source_files.json').read_text())
    n=0
    for item in manifest:
        if item['status']!='deduplicated:shared/reference_amplitudes.json':continue
        dest=ROOT/'studies'/item['stage']/item['path']
        if dest.is_file():continue
        dest.parent.mkdir(parents=True,exist_ok=True)
        try:dest.symlink_to(os.path.relpath(source,dest.parent))
        except (OSError,NotImplementedError):shutil.copyfile(source,dest)
        n+=1
    return n


def install_large_data(source):
    """Install optional exact genus-1000 coefficients in the old v0.1 layout."""
    import zipfile
    dest=ROOT/'studies/01_v01_low_x/results'
    installed=[]
    names=('direct_J4_exact.json.gz','all_constant_and_leading.json.gz')
    source=Path(source)
    if source.is_dir():
        for name in names:
            p=source/name
            if p.is_file():
                target=dest/('gmp' if name.startswith('direct_') else '')/name
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(p,target);installed.append(str(target))
    else:
        with zipfile.ZipFile(source) as z:
            for name in names:
                found=[n for n in z.namelist() if Path(n).name==name]
                if len(found)!=1:raise ValueError(f'{name} absent or ambiguous')
                target=dest/('gmp' if name.startswith('direct_') else '')/name
                target.parent.mkdir(parents=True,exist_ok=True)
                with z.open(found[0]) as inp,target.open('wb') as out:shutil.copyfileobj(inp,out)
                installed.append(str(target))
    return installed


def execute(stage_id,mode='test',timeout=300,verbose=False):
    item=resolve_stage(stage_id)
    cwd=ROOT/'studies'/item['id']
    prepare()
    cmd=[sys.executable]+list(item[mode])
    env=os.environ.copy()
    paths=[str(ROOT/'src'),str(cwd),str(cwd/'code'),str(cwd/'quintic_ai_predictor')]
    env['PYTHONPATH']=os.pathsep.join(paths+[env.get('PYTHONPATH','')])
    started=time.monotonic()
    if verbose:print('RUN',item['id'],mode,' '.join(cmd),flush=True)
    try:
        ret=subprocess.run(cmd,cwd=cwd,env=env,text=True,capture_output=True,timeout=timeout,check=False)
        returncode=ret.returncode
        stdout=ret.stdout
        stderr=ret.stderr
    except subprocess.TimeoutExpired as e:
        returncode=124;stdout=(e.stdout or b'').decode(errors='replace') if isinstance(e.stdout,bytes) else (e.stdout or '')
        stderr='Timeout exceeded' + ((e.stderr or b'').decode(errors='replace') if isinstance(e.stderr,bytes) else (e.stderr or ''))
    report={'stage':item['id'],'action':mode,'status':'passed' if returncode==0 else 'failed',
            'returncode':returncode,'elapsed_s':round(time.monotonic()-started,3),
            'command':' '.join(cmd),'stdout_tail':stdout[-2500:],'stderr_tail':stderr[-2500:],
            'original_archive':item['archive'],'comparison':'historical archived study; refer to its model/validator split'}
    logdir=ROOT/'run_logs';logdir.mkdir(exist_ok=True)
    (logdir/f'{item["id"]}_{mode}.log').write_text(stdout+'\n--- STDERR ---\n'+stderr)
    if verbose:print(report['status'],item['id'],'elapsed',report['elapsed_s'],report['stdout_tail'][-200:],flush=True)
    return report


def execute_all(mode='test',limit_minutes=20,stages_subset=None,continue_on_error=True,verbose=True):
    report=[]
    selected=stages_subset or [x['id'] for x in stages()]
    for stage_id in selected:
        r=execute(stage_id,mode=mode,timeout=min(600,60*limit_minutes),verbose=verbose)
        report.append(r)
        if not continue_on_error and r['returncode']!=0:break
    outfile=ROOT/'run_logs'/f'latest_{mode}_all.json'
    outfile.parent.mkdir(exist_ok=True)
    outfile.write_text(json.dumps(report,indent=2)+'\n')
    return report


def audit_sources():
    """Recompute SHA256 for unmodified bundled historical inputs/code.

    Historical scripts may rewrite results during new runs, so distinguish
    missing, mutated, and matching files. Does not judge physics correctness.
    """
    manifest=json.loads((ROOT/'shared/source_files.json').read_text())
    counts={'matching':0,'changed_after_reproduction':0,'missing':0,'shared':0,'optional':0}
    for record in manifest:
        if record['status'].startswith('deduplicated'):
            counts['shared']+=1;continue
        if record['status'].startswith('optional:'):
            counts['optional']+=1;continue
        if record['status']!='bundled':continue
        p=ROOT/'studies'/record['stage']/record['path']
        if not p.is_file():counts['missing']+=1
        elif hashlib.sha256(p.read_bytes()).hexdigest()!=record['sha256']:counts['changed_after_reproduction']+=1
        else:counts['matching']+=1
    return counts
