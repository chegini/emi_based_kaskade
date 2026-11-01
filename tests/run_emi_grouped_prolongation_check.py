#!/usr/bin/env python3
"""One-step 10cells experiment: fixed P16 or P32 groups, R8/P8, one restriction pass."""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess

from analyze_emi_potential_diagnostics import read, error_metrics
from collect_emi_scaling_smoke import collect
from run_emi_flexible_cg_check import MARKER
from run_emi_low_bit_one_pass import value
from run_emi_mixed_precision_validation import digest, direct_errors, logical_bytes
from run_emi_sd_cg_comparison import make_command as sd_command, set_value

def modes_for(group_size):
    if group_size not in (16, 32):
        raise ValueError('Controlled experiment supports group sizes 16 and 32')
    return ('direct', 'full-cg16', 'full-fcg16', f'full-fcg8-pgroup{group_size}')


def group_marker(group_size):
    return (f'BDDC prolongation groups: producer-order size={group_size}; float32 scale/minimum; '
            'descriptors included in byte counters')


MODES = modes_for(32)
GROUP_MARKER = group_marker(32)


def make_command(saved, mode, out, group_size=32):
    modes = modes_for(group_size)
    if mode not in modes:
        raise ValueError('Unknown experiment mode')
    # Validate the diagnostic source; do not silently change its physical problem.
    expected = {'-np': '8', '--refine': '2', '--bddcIterations': '300',
                '--bddcCgOriginalResidual': '1', '--bddcFlexibleCg': '0',
                '--bddcConstantCorrection': '0', '--bddcUseCg': '1',
                '--bddcRestrictZeroPreserving': '0', '--bddcProlongateZeroPreserving': '0',
                '--bddcCompressionMinInterfaceSize': '0'}
    for flag in ('--bddcGraphLifting', '--bddcGraphLiftingAfterQuantization',
                 '--bddcZigzag', '--bddcBitlength', '--bddcHuffman'):
        expected[flag] = '1'
    for flag, wanted in expected.items():
        if value(saved, flag) != wanted:
            raise ValueError('Incompatible source control: '+flag)
    for flag in ('--bddcProlongateGroupSize', '--bddcCgDiagnostics',
                 '--bddcProlongationDiagnostics', '--bddcProlongationConstraintCall',
                 '--bddcProlongationGroupingCall'):
        if flag in saved and value(saved, flag) != '0':
            raise ValueError('Source must have grouping/observers disabled: '+flag)
    dt = float(value(saved, '--dt'))
    if not math.isfinite(dt) or dt <= 0:
        raise ValueError('Invalid source time step')
    candidate = mode == modes[-1]
    direct = mode == 'direct'
    cmd = sd_command(saved, 'full-fcg8' if candidate else 'full-cg16', out)
    changes = {'--finalTime': value(saved, '--dt'), '--maximumNumberOfTimeSteps': 1,
               '--bddcFlexibleCg': int(mode in modes[2:]),
               '--bddcProlongateGroupSize': group_size if candidate else 0,
               '--bddcProlongationDiagnostics': 0, '--bddcProlongationConstraintCall': 0,
               '--bddcProlongationGroupingCall': 0, '--potentialDiagnostics': 1}
    if direct:
        changes.update({'-np': 1, '--direct': 1, '--bddc': 0, '--directBackend': 'mumps-spd',
                        '--bddcCompression': 0, '--bddcRelativeTolerance': 0,
                        '--bddcMaterialAssembly': 0, '--bddcInterfaceAssembly': 0,
                        '--bddcCgOriginalResidual': 0, '--bddcCompressionReport': 0})
        for flag in ('--bddcGraphLifting', '--bddcGraphLiftingAfterQuantization',
                     '--bddcZigzag', '--bddcBitlength', '--bddcHuffman'):
            changes[flag] = 0
    for flag, new in changes.items():
        set_value(cmd, flag, new)
    return cmd


def first_residual(text):
    rows = re.findall(r'^BDDC global residual step 1: ([^\r\n]+)$', text, re.M)
    if len(rows) != 1:
        raise ValueError('Expected exactly one first-step true residual summary')
    return rows[0]


def save(root, records, group_size=32):
    lines = ['10cells, ONE time step, one restriction pass, true tolerance 1e-7.',
             'mode                      status  calls  logical bytes  saved vs CG16 / FCG16  raw rel L2 / centered / max abs']
    for mode, row in records.items():
        passed = row['status'] == 'PASS'
        count = logical_bytes(row) if passed and mode != 'direct' else None
        savings = {}
        for baseline in ('full-cg16', 'full-fcg16'):
            control = records.get(baseline, {})
            b = logical_bytes(control) if control.get('status') == 'PASS' else None
            savings[baseline] = 100*(1-count/b) if count is not None and b else None
        row['saved_percent'] = savings
        row['beats_both_byte_baselines'] = all(v is not None and v > 0 for v in savings.values())
        calls = (sum(row['metrics']['iterations_per_step']) if passed and mode != 'direct'
                 else row.get('failure_solver_call', '--'))
        fmt = lambda x: '--' if x is None else f'{x:.4e}'
        err = row.get('full_precision_vs_direct', {})
        ss = ' / '.join('--' if v is None else f'{v:.2f}%' for v in savings.values())
        ee = ' / '.join(fmt(err.get(k)) for k in ('relative_l2', 'centered_relative_l2', 'max_absolute_error'))
        lines.append(f'{mode:25s} {row["status"]:6s} {str(calls):>5s} {str(count):>14s} {ss:>23s}  {ee}')
        lines.extend('  '+e for e in row['errors'])
    lines += ['', 'Fresh direct, CG16 and FCG16 controls; same mesh, physics and initial state.',
              f'Candidate: R8 exact-zero; affine P8 groups of {group_size} producer-union values, excluding exact corners.',
              'Scale and minimum remain float32. Group counts, descriptors, padding and exact corners count toward bytes.',
              'No extra restriction pass, precision fallback, tolerance relaxation or post-solve constant correction.',
              'PASS checks convergence and finite full-precision errors; accuracy has no hidden acceptance threshold.',
              'Logical bytes include intra-rank deliveries; exclude MPI framing and codebook training.',
              'Solver/profile times are retained in JSON; this validation includes setup/output and is not a timing benchmark.',
              'This experiment stops after one step. No ten-step campaign is launched.']
    (root/'summary.json').write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--reference-campaign', type=Path, required=True,
                   help='Saved emi-sd-cg8 campaign containing 10cells/full-cg16')
    p.add_argument('--output', type=Path)
    p.add_argument('--group-size', type=int, choices=(16, 32), default=32,
                   help='fixed producer prolongation group size; default 32 preserves the previous experiment')
    args = p.parse_args()
    modes = modes_for(args.group_size)
    repo = args.repo.resolve(); ref = args.reference_campaign.resolve()
    source = ref/'10cells/full-cg16'
    saved = shlex.split((source/'command.sh').read_text())
    source_log = (source/'run.log').read_text(errors='replace')
    source_metrics = collect(source_log, int(value(saved, '--maximumNumberOfTimeSteps')), 8, 0, 1e-7)
    if source_metrics['errors']:
        p.error('Source CG16 control is not a complete passing run: '+str(source_metrics['errors']))
    source_first = first_residual(source_log)
    model = repo/'emiModel'; sha = digest(model)
    if str(model) not in saved:
        p.error('Source executable path differs from this repository')
    for mode in modes:
        make_command(saved, mode, Path('validation-only'), args.group_size)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = (args.output or repo/'output'/f'emi-grouped-prolongation-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root/'environment.json').write_text(json.dumps(dict(
        executable_sha256=sha, reference_campaign=str(ref), source_command=saved,
        source_command_sha256=digest(source/'command.sh'), source_log_sha256=digest(source/'run.log'),
        steps=1, bddc_ranks=8, direct_ranks=1, tolerance=1e-7, restriction_passes=1,
        prolongation_group_size=args.group_size, metadata='float32 scale/minimum', controls_reused=False,
        direct_reference_reused=False, constant_correction=False), indent=2)+'\n')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE', None)
    records = {}; reference = None
    print('Grouped prolongation one-step experiment:', root, flush=True)
    for mode in modes:
        if digest(model) != sha:
            raise RuntimeError('Executable changed during experiment')
        folder = root/'10cells'/mode; folder.mkdir(parents=True)
        cmd = make_command(saved, mode, folder, args.group_size)
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print('Running '+mode+' ...', flush=True)
        with (folder/'run.log').open('w') as stream:
            rc = subprocess.run(cmd, cwd=repo, env=env, stdout=stream, stderr=subprocess.STDOUT).returncode
        text = (folder/'run.log').read_text(errors='replace')
        metrics = {} if mode == 'direct' else collect(text, 1, 8, rc, 1e-7)
        errors = direct_errors(text, 1, rc) if mode == 'direct' else list(metrics['errors'])
        row = dict(status='FAIL', errors=errors, metrics=metrics, return_code=rc)
        if text.count(MARKER) != int(mode in modes[2:]):
            errors.append('Flexible-CG marker does not match requested solver')
        if (text.count(group_marker(args.group_size)) != int(mode == modes[-1]) or
                text.count('BDDC prolongation groups:') != int(mode == modes[-1])):
            errors.append('Grouped prolongation marker does not match requested mode')
        failures = re.findall(r'BDDC step (\d+) failed after (\d+) solver calls: ([^\r\n]+)', text)
        if failures:
            step, calls, reason = failures[-1]
            row.update(failure_step=int(step), failure_solver_call=int(calls))
            errors.insert(0, reason)
        if not errors:
            try:
                state = read(folder/'potential-state-step1.csv')['u']
                initial = read(folder/'potential-state-step0.csv')['u']
                if mode == 'direct':
                    reference = state; initial_reference = initial
                else:
                    if initial != initial_reference:
                        raise ValueError('Initial state differs from fresh direct solve')
                    if len(state) != metrics['dofs']:
                        raise ValueError('Snapshot size differs from logged DOF count')
                    row['full_precision_vs_direct'] = error_metrics(state, reference)
                    if any(v is None or not math.isfinite(v) for v in row['full_precision_vs_direct'].values()):
                        raise ValueError('Nonfinite/undefined solution comparison')
                    if mode == 'full-cg16':
                        row['source_first_step_parity'] = first_residual(text) == source_first
                        if not row['source_first_step_parity']:
                            errors.append('Fresh CG16 first-step residual/iterations differ from saved control')
            except (OSError, ValueError) as e:
                errors.append(str(e))
        row['status'] = 'FAIL' if errors else 'PASS'
        records[mode] = row
        save(root, records, args.group_size)
        (folder/'metrics.json').write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')
        print('  '+row['status'], flush=True)
        if errors and mode != modes[-1]:
            print('Control failed; candidate was not run.', flush=True)
            break
    if digest(model) != sha:
        raise RuntimeError('Executable changed during experiment')
    print('\n'+(root/'summary.txt').read_text())
    print('Results:', root)
    return int(len(records) != len(modes) or any(r['status'] != 'PASS' for r in records.values()))


if __name__ == '__main__':
    raise SystemExit(main())
