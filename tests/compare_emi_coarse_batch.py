#!/usr/bin/env python3
"""Compare new coarse-batch validation with archived scalar-reduction runs.

Checks campaign completeness, convergence, inputs, layout and snapshot integrity.
Reports numerical changes without inventing a scientific solution-error threshold.
"""
import argparse
import json
import math
from pathlib import Path

from analyze_emi_potential_diagnostics import read, error_metrics
from run_emi_histogram_start_check import require
from run_emi_mpi_routing_check import digest


REFERENCES = {
    32: 'emi-solver-study-weak_32_h16-before-20261108T203929Z-2844236',
    64: 'emi-solver-study-weak_64_h16-before-20261108T203929Z-3720694',
    128: 'emi-solver-study-weak_128_h16-before-20261108T203929Z-441739',
}
CASES = ('fcg-uncompressed', 'fcg-before-12-split-on', 'fcg-before-24-split-on')


def compare(before, after):
    configs = [json.loads((p/'environment.json').read_text()) for p in (before, after)]
    summaries = [json.loads((p/'summary.json').read_text()) for p in (before, after)]
    for config in configs:
        require(config['solvers'] == ['fcg'] and config['bits'] in ([12, 16], [12, 24]) and
                config['split'] == 'on' and config['transform_order'] == 'before',
                'Unexpected validation settings')
    for key in ('case', 'ranks', 'expected_nodes', 'iterations', 'steps',
                'input_sha256', 'slurm_btl', 'slurm_vader_single_copy'):
        require(configs[0][key] == configs[1][key], 'Configuration changed: '+key)
    for config, summary in zip(configs, summaries):
        require(summary.get('complete') and not summary.get('failed_cases') and
                not summary.get('fatal_error'), 'Incomplete/failed validation')
        for case in ('direct', 'fcg-uncompressed',
                     *(f'fcg-before-{b}-split-on' for b in config['bits'])):
            require(summary['runs'][case]['status'] == 'PASS', 'Failed case: '+case)
    old_initial = read(before/'direct/potential-state-step0.csv')['u']
    old_direct = read(before/'direct/potential-state-step1.csv')['u']
    new_direct = read(after/'direct/potential-state-step1.csv')['u']
    rows = []
    cases = ['fcg-uncompressed'] + [f'fcg-before-{b}-split-on' for b in configs[1]['bits']]
    for case in cases:
        old, new = summaries[0]['runs'].get(case), summaries[1]['runs'][case]
        control = old or summaries[0]['runs']['fcg-uncompressed']
        for key in ('cells', 'dofs', 'subdomains', 'coarse_interfaces'):
            require(control[key] == new[key], 'Problem changed: '+key)
        for root in ((before, after) if old else (after,)):
            require(read(root/case/'potential-state-step0.csv')['u'] == old_initial,
                    'Initial state changed: '+case)
        require(read(after/'direct/potential-state-step0.csv')['u'] == old_initial,
                'Direct initial state changed')
        u_old = read(before/case/'potential-state-step1.csv')['u'] if old else None
        u_new = read(after/case/'potential-state-step1.csv')['u']
        old_error = error_metrics(u_old, old_direct) if old else None
        new_error = error_metrics(u_new, new_direct)
        for error in ([old_error, new_error] if old else [new_error]):
            require(all(error[k] is not None and math.isfinite(error[k])
                        for k in ('relative_l2', 'centered_relative_l2')),
                    'Invalid solution-error comparison')
        for record in ([old, new] if old else [new]):
            require(len(record['relative_residuals']) == 1 and
                    0 <= record['relative_residuals'][0] <= 1e-7,
                    'True residual tolerance not met')
        rows.append(dict(case=case, baseline_available=old is not None,
            old_iterations=old['iterations_per_step'] if old else None,
            new_iterations=new['iterations_per_step'], old_error=old_error,
            new_error=new_error, change_from_scalar=error_metrics(u_new, u_old) if old else None,
            old_relative_residual=old['relative_residuals'][0] if old else None,
            new_relative_residual=new['relative_residuals'][0],
            old_payload=old['transmitted_sent_bytes'] if old else None, new_payload=new['transmitted_sent_bytes']))
    return dict(mesh=configs[0]['case'], baseline=str(before), candidate=str(after),
        baseline_summary_sha256=digest(before/'summary.json'),
        candidate_summary_sha256=digest(after/'summary.json'),
        baseline_executable_sha256=configs[0]['executable_sha256'],
        candidate_executable_sha256=configs[1]['executable_sha256'],
        direct_change=error_metrics(new_direct, old_direct), rows=rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign', type=Path, required=True)
    a = p.parse_args(argv)
    repo = Path(__file__).resolve().parents[1]
    root = a.campaign.resolve()
    results = [compare(repo/'output'/baseline, root/f'validation-weak_{size}_h16')
               for size, baseline in REFERENCES.items()]
    note = ('Convergence and consistency checks passed. Errors and trajectory changes '
            'require review; no scientific solution-error acceptance threshold is imposed.')
    (root/'comparison.json').write_text(json.dumps(dict(note=note, results=results),
                                                  indent=2, allow_nan=False)+'\n')
    print(note)
    print('| Mesh | Case | Iterations old/new | Relative error/direct old/new | '
          'Centered error old/new | Relative change from scalar | True residual new |')
    print('|---|---|---|---|---|---|---|')
    for result in results:
        for row in result['rows']:
            def number(value):
                return '--' if value is None else f'{value:.6g}'
            old_error = row['old_error'] or {}
            delta = row['change_from_scalar'] or {}
            old_iterations = row['old_iterations'][0] if row['old_iterations'] else None
            print(f"| {result['mesh']} | {row['case']} | "
                  f"{number(old_iterations)}/{row['new_iterations'][0]} | "
                  f"{number(old_error.get('relative_l2'))}/{number(row['new_error']['relative_l2'])} | "
                  f"{number(old_error.get('centered_relative_l2'))}/"
                  f"{number(row['new_error']['centered_relative_l2'])} | "
                  f"{number(delta.get('relative_l2'))} | "
                  f"{row['new_relative_residual']:.6g} |")
    print('-- means no matching scalar baseline is available; 16 bits is not compared with 24 bits.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
