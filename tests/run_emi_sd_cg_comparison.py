#!/usr/bin/env python3
"""Test existing steepest descent, CG and one-term FCG at matched full-coded 8/8."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

from analyze_emi_potential_diagnostics import read, error_metrics
from collect_emi_scaling_smoke import collect
from run_emi_flexible_cg_check import MARKER
from run_emi_low_bit_one_pass import value
from run_emi_mixed_precision_validation import digest, logical_bytes

# mode -> compression, solver, precision; 16-bit controls retain their original grid.
MODES = {
    'uncompressed-cg': (False, 'cg', 16),
    'uncompressed-sd': (False, 'sd', 16),
    'full-cg16': (True, 'cg', 16),
    'full-cg8': (True, 'cg', 8),
    'full-fcg8': (True, 'fcg', 8),
    'full-sd8': (True, 'sd', 8),
}


def set_value(cmd, flag, new):
    if cmd.count(flag) > 1:
        raise ValueError('Duplicate flag '+flag)
    if flag in cmd:
        cmd[cmd.index(flag)+1] = str(new)
    else:
        cmd += [flag, str(new)]


def make_command(saved, mode, out, iterations=300):
    """Keep the saved physical problem; explicitly normalize CG-only options."""
    compressed, solver, bits = MODES[mode]
    cmd = list(saved)
    required = {'--mpi': '1', '--sdc': '0', '--bddc': '1', '--direct': '0',
                '--bddcRelativeTolerance': '1e-7', '--bddcRestrictionPasses': '1',
                '--bddcRestrictGroupSize': '0', '--potentialDiagnostics': '1',
                '--bddcCompression': '1', '--bddcRestrictCompressionBits': '16',
                '--bddcProlongateCompressionBits': '16'}
    for flag, wanted in required.items():
        if value(cmd, flag) != wanted:
            raise ValueError('Incompatible saved control: '+flag)
    if iterations < 1:
        raise ValueError('Iteration limit must be positive')
    changes = {
        '--bddcUseCg': int(solver != 'sd'), '--bddcFlexibleCg': int(solver == 'fcg'),
        '--bddcCgOriginalResidual': int(solver != 'sd'),
        '--bddcConstantCorrection': 0, '--bddcCgDiagnostics': 0,
        '--bddcConvergenceReport': 0, '--bddcConditionEstimate': 0,
        '--bddcIterations': iterations, '--bddcCompression': int(compressed),
        '--bddcCompressionBits': bits, '--bddcRestrictCompressionBits': bits,
        '--bddcProlongateCompressionBits': bits,
        '--bddcRestrictZeroPreserving': int(compressed and bits == 8),
        '--bddcProlongateZeroPreserving': 0, '--dir': out,
    }
    for flag in ['--bddcGraphLifting', '--bddcGraphLiftingAfterQuantization',
                 '--bddcZigzag', '--bddcBitlength', '--bddcHuffman']:
        changes[flag] = int(compressed)
    for flag, new in changes.items():
        set_value(cmd, flag, new)
    return cmd


def save(root, records):
    lines = ['Full compression; one restriction pass; true relative tolerance 1e-7.',
             'case         mode              status  calls  logical bytes  saved vs CG16  raw rel L2  centered rel L2']
    for case, runs in records.items():
        baseline = runs.get('full-cg16', {})
        for mode, row in runs.items():
            passed = row['status'] == 'PASS'
            calls = (sum(row['metrics']['iterations_per_step']) if passed
                     else row.get('failure_solver_call', '--'))
            count = logical_bytes(row) if passed else None
            saving = (100*(1-count/logical_bytes(baseline))
                      if passed and baseline.get('status') == 'PASS' else None)
            row['saved_vs_cg16_percent'] = saving
            err = row.get('full_precision_vs_direct', {})
            fmt = lambda x: '--' if x is None else f'{x:.3e}'
            lines.append(f'{case:12s} {mode:17s} {row["status"]:6s} {str(calls):>5s} '
                         f'{str(count) if count is not None else "--":>14s} '
                         f'{f"{saving:.2f}%" if saving is not None else "--":>14s} '
                         f'{fmt(err.get("relative_l2")):>11s} '
                         f'{fmt(err.get("centered_relative_l2")):>16s}')
            for error in row['errors']:
                lines.append('  '+error)
    lines += ['', 'Passing calls are totals over all steps; failing calls refer to the failed step.',
              'All 8/8 runs use exact-zero restriction and ORIGINAL prolongation quantization.',
              'SD uses the existing restricted-residual step numerator; CG/FCG use original-residual rho.',
              'Post-solve constant correction and CG-only diagnostics are OFF in every run.',
              'PASS requires every step to converge and finite full-precision final comparisons.',
              'Raw/centered coefficient errors are reported separately, without an accuracy acceptance threshold.',
              'Logical bytes include intra-rank transfers; this is not measured network traffic.',
              'More iterations are allowed within the selected cap; compare TOTAL bytes against fresh CG16.',
              'Timings include setup/output and are saved in JSON; no speedup claim is made.',
              'No precision fallback, extra restriction pass, tolerance change or numerical-code modification.']
    (root/'summary.json').write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--reference-campaign', type=Path, required=True)
    p.add_argument('--cases', nargs='+', default=['10cells'])
    p.add_argument('--iterations', type=int, default=300,
                   help='same maximum solver calls per step for all modes; default 300')
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if args.iterations < 1 or len(set(args.cases)) != len(args.cases):
        p.error('Use a positive iteration limit and distinct cases')
    repo = args.repo.resolve(); ref = args.reference_campaign.resolve()
    model = repo/'emiModel'; sha = digest(model)
    oldenv = json.loads((ref/'environment.json').read_text())
    old = json.loads((ref/'summary.json').read_text())
    steps = oldenv['steps']; ranks = oldenv['bddc_ranks']
    seeds = {}; direct = {}
    for case in args.cases:
        for mode in ['direct', 'full-r16-p16']:
            if old.get(case, {}).get('runs', {}).get(mode, {}).get('status') != 'PASS':
                p.error('Missing passing reference: '+case+'/'+mode)
        seeds[case] = shlex.split((ref/case/'full-r16-p16/command.sh').read_text())
        for mode in MODES:
            cmd = make_command(seeds[case], mode, Path('validation-only'), args.iterations)
            if str(model) not in cmd or value(cmd, '-np') != str(ranks):
                p.error('Reference executable path or rank count differs')
            if value(cmd, '--maximumNumberOfTimeSteps') != str(steps):
                p.error('Reference step count differs')
        direct[case] = read(ref/case/'direct'/f'potential-state-step{steps}.csv')['u']
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = (args.output or repo/'output'/f'emi-sd-cg8-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root/'environment.json').write_text(json.dumps(dict(
        executable_sha256=sha, reference_campaign=str(ref), reference_environment=oldenv,
        cases=args.cases, steps=steps, ranks=ranks, tolerance=1e-7, restriction_passes=1,
        iterations=args.iterations, modes=MODES, constant_correction=False,
        restriction_quantizer_8='symmetric exact-zero', prolongation_quantizer='original',
        controls_reused=False, direct_reference_reused=True), indent=2)+'\n')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE', None)
    records = {}
    print('SD/CG comparison output:', root, flush=True)
    for case in args.cases:
        records[case] = {}
        for mode in MODES:
            if digest(model) != sha:
                raise RuntimeError('Executable changed during comparison')
            folder = root/case/mode; folder.mkdir(parents=True)
            cmd = make_command(seeds[case], mode, folder, args.iterations)
            (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
            print(f'Running {case}: {mode}, {steps} steps ...', flush=True)
            with (folder/'run.log').open('w') as stream:
                rc = subprocess.run(cmd, cwd=repo, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT).returncode
            text = (folder/'run.log').read_text(errors='replace')
            metrics = collect(text, steps, ranks, rc, 1e-7)
            row = dict(status='FAIL', metrics=metrics, return_code=rc, errors=list(metrics['errors']))
            if text.count(MARKER) != int(mode == 'full-fcg8'):
                row['errors'].append('Flexible-CG marker does not match requested solver')
            failures = re.findall(r'BDDC step (\d+) failed after (\d+) solver calls: ([^\r\n]+)', text)
            if failures:
                step, calls, reason = failures[-1]
                row.update(failure_step=int(step), failure_solver_call=int(calls), failure_reason=reason)
                row['errors'].insert(0, f'Step {step}: {reason}')
            if not row['errors']:
                try:
                    row['full_precision_vs_direct'] = error_metrics(
                        read(folder/f'potential-state-step{steps}.csv')['u'], direct[case])
                except (ValueError, OSError) as error:
                    row['errors'].append(str(error))
            row['status'] = 'FAIL' if row['errors'] else 'PASS'
            records[case][mode] = row
            save(root, records)
            (folder/'metrics.json').write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')
            print('  '+row['status'], flush=True)
            if row['status'] == 'FAIL' and mode in ('uncompressed-cg', 'uncompressed-sd', 'full-cg16'):
                print('Control failed; skipping remaining candidates for this case.', flush=True)
                break
    if digest(model) != sha:
        raise RuntimeError('Executable changed during comparison')
    print('\n'+(root/'summary.txt').read_text())
    print('Results:', root)
    return int(any(len(runs) != len(MODES) or any(r['status'] != 'PASS' for r in runs.values())
                   for runs in records.values()))


if __name__ == '__main__':
    raise SystemExit(main())
