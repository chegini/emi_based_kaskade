#!/usr/bin/env python3
"""One controlled grouped R12/P16 experiment against unchanged R16/P16."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
from run_emi_zero_preserving_check import assess


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    p.add_argument('--reference',type=Path,
                   default=Path('output/emi-zero-preserving-20261101T081725Z/baseline16'))
    args=p.parse_args()
    repo=Path(__file__).resolve().parents[1]
    reference=args.reference if args.reference.is_absolute() else repo/args.reference
    for case in ['uncompressed-cg','quantization-16','full-16']:
        if not (reference/case/'bddc-convergence-step1.csv').is_file():
            p.error('Missing preserved baseline CSV: '+str(reference/case))
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=(args.output or repo/'output'/f'emi-grouped-restriction-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True,exist_ok=False)
    common=[sys.executable,str(repo/'tests/run_emi_quantization_diagnostic.py'),
            '--ranks','8','--refine','2','--restriction-passes','1',
            '--tolerance','1e-7','--iterations','300','--full-compression','--cg-diagnostics']
    summaries={}
    for name,bits,flags in [
        ('baseline16',16,[]),
        ('group128-r12-p16',12,['--zero-preserving-restriction',
                              '--restriction-group-size','128','--prolongation-bits','16'])]:
        result=subprocess.run(common+['--bits',str(bits),'--output',str(root/name)]+flags,cwd=repo)
        rows={r['case']:r for r in json.loads((root/name/'summary.json').read_text())}
        summaries[name]=rows
        if any(rows.get(c,{}).get('status')!='PASS' for c in ['direct','uncompressed-cg']):
            raise RuntimeError('Reference failed; experiment is not comparable: '+name)
        if name=='baseline16':
            if result.returncode: raise RuntimeError('16-bit baseline failed')
            for case in ['uncompressed-cg','quantization-16','full-16']:
                csv='bddc-convergence-step1.csv'
                if (root/name/case/csv).read_bytes()!=(reference/case/csv).read_bytes():
                    raise RuntimeError('Disabled-group baseline changed: '+case)
            print('Disabled-group baseline parity passed.',flush=True)
    baseline=summaries['baseline16']['full-16']
    candidate=summaries['group128-r12-p16']['full-12']
    failures=assess(baseline,candidate)
    for name,case,r,pbits,group,zero in [
        ('baseline16','full-16',16,16,0,False),
        ('group128-r12-p16','full-12',12,16,128,True)]:
        for variant in [case,case.replace('full','quantization')]:
            row=summaries[name].get(variant,{})
            for key,value in [('restriction_bits',r),('prolongation_bits',pbits),
                              ('restriction_group_size',group),('restriction_passes',1),
                              ('zero_preserving_restriction',zero)]:
                if row.get(key)!=value: failures.append(name+': configuration mismatch '+key)
        csv='bddc-convergence-step1.csv'
        a=root/name/case/csv
        b=root/name/case.replace('full','quantization')/csv
        if not a.is_file() or not b.is_file() or a.read_bytes()!=b.read_bytes():
            failures.append(name+': quantization/full convergence histories differ')
    env=[json.loads((root/n/'environment.json').read_text()) for n in summaries]
    if env[0]['executable_sha256']!=env[1]['executable_sha256']:
        failures.append('Executable changed during experiment')
    old=baseline.get('logical_transferred_bytes')
    new=candidate.get('logical_transferred_bytes')
    report=dict(objective_passed=not failures,failures=failures,baseline=baseline,candidate=candidate,
                bytes_ratio=new/old if old and new is not None else None,
                note='Producer groups=128, R12/P16, one pass; CG unchanged. '
                     'All group metadata and slot padding included in logical counters. '
                     'Times include diagnostics and packet copies, not a performance benchmark.')
    (root/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('\nGrouped single-pass objective:', 'PASS' if not failures else 'NOT MET')
    print('Total full-codec logical bytes: baseline=',old,' candidate=',new)
    for failure in failures: print('  -',failure)
    print('Results:',root)
    return int(bool(failures))


if __name__=='__main__':
    sys.exit(main())
