#!/usr/bin/env python3
"""Controlled restriction-only experiment; no change to CG, physics or tolerance."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys


def assess(baseline, candidate):
    """Require convergence, no worse raw solution error, and less total traffic."""
    failures = []
    for label, row in [('baseline',baseline),('candidate',candidate)]:
        if row.get('status') != 'PASS':
            failures.append(label+': '+row.get('reason','solver/comparison failed'))
        if row.get('restriction_passes') != 1:
            failures.append(label+': restriction pass count is not one')
        rel=row.get('relative_residual')
        if rel is None or not 0<=rel<=1e-7:
            failures.append(label+': true relative residual exceeds 1e-7 or is missing')
    if failures:
        return failures
    for key in ['relative L2 error','max absolute error']:
        ref=baseline.get('comparison_to_direct',{}).get(key)
        value=candidate.get('comparison_to_direct',{}).get(key)
        if ref is None or value is None or not 0<=value<=ref:
            failures.append('candidate raw exported-u '+key+' is missing or worse than baseline')
    ref=baseline.get('logical_transferred_bytes')
    value=candidate.get('logical_transferred_bytes')
    if ref is None or value is None or not 0<value<ref:
        failures.append('candidate total logical bytes do not beat the matched full-16 baseline')
    if value is None or not 0<value<2430000:
        failures.append('candidate does not beat the historical 2,430,000-byte target')
    return failures


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1]
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=args.output or repo/'output'/('emi-zero-preserving-'+stamp)
    root=root.resolve();root.mkdir(parents=True,exist_ok=False)
    runner=repo/'tests/run_emi_quantization_diagnostic.py'
    common=[sys.executable,str(runner),'--ranks','8','--refine','2',
            '--restriction-passes','1','--tolerance','1e-7','--iterations','300',
            '--full-compression','--cg-diagnostics']
    summaries={}
    for name,bits,extra in [('baseline16',16,[]),
                            ('zero12',12,['--zero-preserving-restriction'])]:
        command=common+['--bits',str(bits),'--output',str(root/name)]+extra
        result=subprocess.run(command,cwd=repo)
        path=root/name/'summary.json'
        if not path.is_file():
            print('Runner failed before writing summary:',root/name,file=sys.stderr)
            return 1
        summaries[name]={r['case']:r for r in json.loads(path.read_text())}
        if name=='baseline16' and result.returncode:
            print('Baseline failed; candidate not launched. Results:',root)
            return 1
        if any(summaries[name].get(case,{}).get('status')!='PASS'
               for case in ['direct','uncompressed-cg']):
            print('Reference failed; comparison not valid. Results:',root)
            return 1
    baseline=summaries['baseline16'].get('full-16',{})
    candidate=summaries['zero12'].get('full-12',{})
    failures=assess(baseline,candidate)
    for name,case in [('baseline16','full-16'),('zero12','full-12')]:
        rows=summaries[name]
        expected=name=='zero12'
        for key in [case,case.replace('full','quantization')]:
            if rows.get(key,{}).get('zero_preserving_restriction')!=expected:
                failures.append(name+': quantizer configuration not confirmed')
        if rows.get(case,{}).get('status')=='PASS':
            first=root/name/case/'bddc-convergence-step1.csv'
            second=root/name/case.replace('full','quantization')/'bddc-convergence-step1.csv'
            if not first.is_file() or not second.is_file() or first.read_bytes()!=second.read_bytes():
                failures.append(name+': quantization/full convergence histories differ')
    envs=[json.loads((root/name/'environment.json').read_text()) for name in summaries]
    if envs[0]['executable_sha256']!=envs[1]['executable_sha256']:
        failures.append('Executable changed between baseline and candidate')
    b=baseline.get('logical_transferred_bytes');c=candidate.get('logical_transferred_bytes')
    report={'objective_passed':not failures,'failures':failures,
            'baseline':baseline,'candidate':candidate,
            'candidate_bytes_vs_baseline':c/b if b and c is not None else None,
            'historical_byte_target':2430000,
            'note':'Diagnostic timings include observer work. Bytes are logical, not hardware traffic. '
                   'Accuracy compares raw exported-u L2 and max differences with direct; no mean removed.'}
    (root/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('\nSingle-pass zero-preserving objective:', 'PASS' if not failures else 'NOT MET')
    print('Baseline full-16 bytes:',b,'; candidate full-12 bytes:',c)
    if b and c is not None: print('Candidate / baseline bytes:',f'{c/b:.6f}')
    for failure in failures:print('  -',failure)
    print('Results:',root)
    return int(bool(failures))


if __name__=='__main__':
    sys.exit(main())
