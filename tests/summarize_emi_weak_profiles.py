#!/usr/bin/env python3
"""Read completed weak-timing campaigns; exclude warm-ups and summarize phase maxima."""
import argparse
import json
from pathlib import Path
import re
import statistics

NUMBER = r'[0-9.eE+\-]+'


def phase_maxima(log, ranks):
    result = {}
    for name, pattern in (
        ('solver/coarse setup', rf'^  BDDC solver/coarse setup: ({NUMBER})'),
        ('BDDC solve calls', rf'^  BDDC solve calls: .*?wall seconds: ({NUMBER})'),
        ('residual assembly', rf'^  residual assembly: ({NUMBER})'),
        ('local factorization', rf'^  subdomain construction/factorization: ({NUMBER})')):
        values = [float(v) for v in re.findall(pattern, log, re.M)]
        if len(values) != ranks:
            raise ValueError(f'{name}: expected {ranks} rank records, found {len(values)}')
        result[name] = max(values)
    for name in ('global system assembly', 'material matrix exchange', 'local matrix extraction'):
        rows = re.findall(rf'^  {re.escape(name)}\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*$', log, re.M)
        if len(rows) != 1:
            raise ValueError(f'{name}: missing or duplicate phase row')
        result[name] = float(rows[0][2])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--job', required=True)
    p.add_argument('--output-dir', type=Path, default=Path('output'))
    a = p.parse_args()
    found = []
    for env_file in a.output_dir.glob('emi-weak-timings-weak*/environment.json'):
        env = json.loads(env_file.read_text())
        if str(env['slurm']['SLURM_JOB_ID']) == a.job:
            found.append((env['configuration']['expected_nodes'], env_file.parent, env))
    if not found:
        p.error('No matching timing campaigns found')
    print('Measured runs only; median [minimum, maximum] of per-run rank maxima, in seconds.')
    print('Phase maxima are not additive. Nested solver counters are not a wall-time partition.')
    for nodes, root, env in sorted(found):
        report = json.loads((root/'summary.json').read_text())
        if not report['complete']:
            raise ValueError('Incomplete campaign: '+str(root))
        grouped = {}
        for row in report['runs']:
            if row['warmup']:
                continue
            if row['status'] != 'PASS':
                raise ValueError('Failed measurement: '+str(row['folder']))
            folder = root/Path(row['folder']).name
            phases = phase_maxima((folder/'run.log').read_text(), env['configuration']['ranks'])
            phases['solver wrapper'] = row['solver_wrapper_max_seconds']
            phases['application'] = row['wall_seconds']
            for phase, value in phases.items():
                grouped.setdefault((row['case'], phase), []).append(value)
        print(f'\n{env["configuration"]["case"]}: {nodes} nodes; {root}')
        print('| Case | Phase | n | Median [min, max] s |')
        print('|---|---|---:|---|')
        for (case, phase), values in sorted(grouped.items()):
            print(f'| {case} | {phase} | {len(values)} | {statistics.median(values):.6f} [{min(values):.6f}, {max(values):.6f}] |')


if __name__ == '__main__':
    main()
