#!/usr/bin/env python3
"""Capture existing trajectories without changing quantization; verify CSV parity."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
from analyze_emi_restriction_groups import analyze


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True, type=Path,
                        help='existing emi-zero-preserving directory with baseline16 and zero12')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    reference = args.reference.resolve()
    for cohort in ['baseline16','zero12']:
        env = json.loads((reference/cohort/'environment.json').read_text())
        expected = dict(mesh_refinements=2, bddc_ranks=8, steps=1,
                        true_relative_tolerance=1e-7, restriction_passes=1,
                        max_iterations=300)
        if any(env.get(k) != v for k,v in expected.items()):
            raise ValueError('Reference parameters differ: '+cohort)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = (args.output or repo/'output'/f'emi-restriction-ranges-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True, exist_ok=False)
    runner = repo/'tests/run_emi_quantization_diagnostic.py'
    common = [sys.executable, str(runner), '--ranks','8','--refine','2',
              '--restriction-passes','1','--tolerance','1e-7','--iterations','300',
              '--cg-diagnostics','--restriction-range-capture']
    parity = {}
    for cohort,bits,flags in [('baseline16',16,[]),
                              ('zero12',12,['--zero-preserving-restriction'])]:
        command = common+['--bits',str(bits),'--output',str(root/cohort)]+flags
        # Expected zero12 failure is retained; references must still pass.
        result = subprocess.run(command, cwd=repo)
        rows = {r['case']:r for r in json.loads((root/cohort/'summary.json').read_text())}
        if any(rows[name]['status'] != 'PASS' for name in ['direct','uncompressed-cg']):
            raise RuntimeError('Reference solve failed; ranges not comparable')
        if cohort == 'baseline16' and result.returncode:
            raise RuntimeError('16-bit baseline failed')
        case = f'quantization-{bits}'
        checks = {}
        for name in ['uncompressed-cg',case]:
            filename = 'bddc-convergence-step1.csv'
            checks[name] = ((root/cohort/name/filename).read_bytes()
                            == (reference/cohort/name/filename).read_bytes())
        parity[cohort] = checks
        (root/'capture-parity.json').write_text(json.dumps(parity,indent=2)+'\n')
        if not all(checks.values()):
            raise RuntimeError('Capture/reference convergence CSV differs; investigate before interpreting')
        capture = root/cohort/case/'restriction-buffers'
        buffers, groups, _ = analyze(capture,[32,64,128,256,512])
        if not buffers:
            raise RuntimeError('No buffer data: check rebuilt library header')
        subprocess.run([sys.executable,str(repo/'tests/analyze_emi_restriction_groups.py'),
                        str(capture),'--output',str(root/cohort/'range-analysis')],check=True)
        # Useful sanity counts, without treating sample totals as full-solve costs.
        print(cohort, 'captured subdomains:', len({r['subdomain'] for r in buffers}),
              'buffer snapshots:',len(buffers),'group rows:',len(groups))
    print('Capture parity passed. No candidate quantization was applied.')
    print('Range measurements and metadata estimates:',root)
    print('Do not use these I/O-instrumented timings for performance comparisons.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
