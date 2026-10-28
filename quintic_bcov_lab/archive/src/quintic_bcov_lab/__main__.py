"""python -m quintic_bcov_lab ..."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from . import __version__
from . import backend,datasets,runner


def main(argv=None):
    ap=argparse.ArgumentParser(prog='python -m quintic_bcov_lab',description='Sequential BCOV / quintic arithmetic research suite')
    sub=ap.add_subparsers(dest='command',required=True)
    sub.add_parser('list',help='List every historical study in chronological order')
    p=sub.add_parser('stage',help='Run an original study prediction or its tests')
    p.add_argument('id',help='Stage ID or unique prefix, e.g. 05 or 18')
    p.add_argument('--action',choices=['run','test'],default='test')
    p.add_argument('--timeout',type=int,default=300)
    sub.add_parser('prepare',help='Restore shared exact reference aliases')
    p=sub.add_parser('all',help='Run stages sequentially, respecting predictor/validator separation')
    p.add_argument('--action',choices=['run','test'],default='test')
    p.add_argument('--timeout-minutes',type=int,default=10)
    p.add_argument('--stages',nargs='*',default=None)
    p.add_argument('--stop-on-error',action='store_true')
    sub.add_parser('backend-check',help='Validate exact BCOV engine and reference gauge')
    sub.add_parser('audit',help='SHA provenance audit of original study files')
    p=sub.add_parser('generate-intermediate',help='Export exact reference-derived intermediate labels')
    p.add_argument('--out',default='run_logs/intermediate_exact_labels.csv')
    p=sub.add_parser('install-large-data',help='Install optional exact low-X genus-1000 data')
    p.add_argument('--source',required=True,help='Directory or optional data ZIP')
    args=ap.parse_args(argv)
    if args.command=='list':
        for s in runner.stages():print(s['id'].ljust(27),s['title'])
        return 0
    if args.command=='prepare':print('reference aliases added:',runner.prepare());return 0
    if args.command=='backend-check':
        print(json.dumps(backend.validate_small(),indent=2));return 0
    if args.command=='audit':
        print(json.dumps(runner.audit_sources(),indent=2));return 0
    if args.command=='generate-intermediate':
        path=runner.ROOT/args.out
        print('written',len(datasets.write_intermediate(path)),'exact targets to',path)
        return 0
    if args.command=='install-large-data':
        print('Installed\n'+'\n'.join(runner.install_large_data(args.source)));return 0
    if args.command=='stage':
        r=runner.execute(args.id,mode=args.action,timeout=args.timeout,verbose=True)
        print(json.dumps(r,indent=2));return 0 if r['status']=='passed' else 1
    if args.command=='all':
        reports=runner.execute_all(args.action,args.timeout_minutes,args.stages,not args.stop_on_error)
        passed=sum(x['status']=='passed' for x in reports)
        print('Stage runs:',passed,'passed /',len(reports),'total')
        return 0 if passed==len(reports) else 1
    return 2

if __name__=='__main__':raise SystemExit(main())
