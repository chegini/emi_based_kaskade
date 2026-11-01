#!/usr/bin/env python3
"""Replay one saved validation run with read-only CG diagnostics; no parameter sweep."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from run_emi_mixed_precision_validation import digest
from summarize_emi_cg_diagnostics import summarize


def observed_command(original, output):
    cmd=list(original)
    for flag,value in [('--bddc','1'),('--direct','0'),('--bddcUseCg','1'),
                       ('--bddcConvergenceReport','1'),('--bddcCgDiagnostics','0')]:
        if cmd.count(flag)!=1 or cmd[cmd.index(flag)+1]!=value:
            raise ValueError('Expected saved CG run with '+flag+' '+value)
    for flag,value in [('--dir',str(output)),('--bddcCgDiagnostics','1')]:
        if cmd.count(flag)!=1: raise ValueError('Missing/duplicate '+flag)
        cmd[cmd.index(flag)+1]=value
    return cmd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True,
                   help='one mode directory from run_emi_mixed_precision_validation.py')
    args=p.parse_args()
    repo=Path(__file__).resolve().parents[1]
    source=args.run.resolve()
    campaign=source.parents[1]
    metadata=json.loads((campaign/'environment.json').read_text())
    previous=json.loads((campaign/'summary.json').read_text())[source.parent.name]['runs'][source.name]
    old=shlex.split((source/'command.sh').read_text())
    models=[Path(x) for x in old if Path(x).name=='emiModel']
    if len(models)!=1 or digest(models[0])!=metadata['executable_sha256']:
        p.error('Executable differs from saved run; do not rebuild for this replay')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out=repo/'output'/f'emi-cg-replay-{source.parent.name}-{stamp}-{os.getpid()}'
    cmd=observed_command(old,out)
    histories=sorted(source.glob('bddc-convergence-step*.csv'))
    if not histories: p.error('No saved convergence CSVs to check observer parity')
    out.mkdir(parents=True,exist_ok=False)
    (out/'command.sh').write_text(shlex.join(cmd)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE',None)
    print('Diagnostic replay:',out,flush=True)
    print('Only --bddcCgDiagnostics and --dir changed; the previous solver failure may recur.',flush=True)
    with (out/'run.log').open('w') as stream:
        rc=subprocess.run(cmd,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT).returncode
    parity={f.name:(out/f.name).is_file() and f.read_bytes()==(out/f.name).read_bytes()
            for f in histories}
    same_files={f.name for f in out.glob('bddc-convergence-step*.csv')}==set(parity)
    diagnostics=[out/f.name.replace('bddc-convergence','bddc-cg-diagnostics') for f in histories]
    collected=all(f.is_file() and f.stat().st_size>0 for f in diagnostics)
    unchanged=digest(models[0])==metadata['executable_sha256']
    same_completion=(rc==0)==(previous['return_code']==0)
    valid=all(parity.values()) and same_files and collected and unchanged and same_completion
    report=dict(source=str(source),executable_sha256=metadata['executable_sha256'],
                previous_return_code=previous['return_code'],return_code=rc,
                convergence_csv_parity=parity,same_history_files=same_files,
                diagnostics_present=collected,executable_unchanged=unchanged,
                diagnostic_replay_valid=valid,
                note='Replay validity is NOT solver convergence. Timings include observer overhead.')
    (out/'replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Application exit code:',rc)
    print('Observer replay parity:', 'PASS' if valid else 'FAIL')
    if valid:
        for f in diagnostics: summarize(f)
    else:
        print('Do not interpret differences as the original failure until parity is resolved.')
    print('Results:',out)
    return int(not valid)


if __name__=='__main__':
    sys.exit(main())
