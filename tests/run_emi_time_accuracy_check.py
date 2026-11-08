#!/usr/bin/env python3
"""Pepe FGMRES trajectory checks on 2/4/8 nodes; raw, compressed and local bypass."""
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
from collect_emi_scaling_smoke import collect, NUMBER
from run_emi_histogram_start_check import require, mesh_inputs, snapshots_equal
from run_emi_mpi_routing_check import digest, make_command, routing_counts
from run_emi_routing_timings import reference_arguments


CASES = ('direct', 'fgmres-uncompressed-destination',
         'fgmres-compressed-split-destination', 'fgmres-compressed-split-local-bypass')


def accuracy_arguments(config, repo, image, kaskade, steps):
    require(1 <= steps <= 10, 'FGMRES split accuracy check supports 1..10 steps')
    require(config['expected_nodes'] in (2, 4, 8) and
            config['ranks'] == 2*config['expected_nodes'],
            'Require 2, 4 or 8 nodes and two ranks per node')
    require(config['steps'] == 1 and config['uncompressed_steps'] == 0 and
            config['split'] == 'on' and config['bits'] == 16 and
            config['transform_order'] == 'before',
            'Require a one-step 16-bit transform-first split reference without warm-up steps')
    a = reference_arguments(config, repo, image, kaskade)
    a.steps = steps
    a.solvers = ['fgmres']
    return a


def accuracy_command(a, folder, name):
    require(name in CASES, 'Unknown accuracy case')
    cmd = make_command(a, folder, ['mpirun'], 'fgmres',
                       name.startswith('fgmres-compressed'), name != 'direct',
                       direct=name == 'direct', bypass=name.endswith('-local-bypass'))
    # Ten physical steps can exceed the one-step launcher's ten-minute limit.
    require(cmd.count('--time=00:10:00') == 1, 'Unexpected Slurm step time limit')
    cmd[cmd.index('--time=00:10:00')] = '--time=00:40:00'
    return cmd


def direct_check(log, exit_code, steps):
    rows = re.findall(rf'^Direct global residual step (\d+): abs=({NUMBER}) rel=({NUMBER})', log, re.M)
    require(exit_code == 0 and 'EMI model completed' in log and
            [int(r[0]) for r in rows] == list(range(1, steps+1)),
            'Direct reference failed completion/residual-count check')
    require(all(math.isfinite(float(r[1])) and float(r[1]) >= 0 and
                math.isfinite(float(r[2])) and 0 <= float(r[2]) <= 1e-7 for r in rows),
            'Direct reference residual is nonfinite, negative or above tolerance')
    return [float(r[2]) for r in rows]


def trajectory_errors(folder, direct, raw, steps):
    result = []
    for step in range(1, steps+1):
        current = read(folder/f'potential-state-step{step}.csv')
        reference = read(direct/f'potential-state-step{step}.csv')
        baseline = read(raw/f'potential-state-step{step}.csv')
        # Each trajectory evolves its own state/RHS after the initial step.
        result.append(dict(step=step, time=step*0.01,
                           vs_direct=error_metrics(current['u'], reference['u']),
                           vs_uncompressed=error_metrics(current['u'], baseline['u'])))
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True,
                   help='Successful one-step routing/bypass validation at the same rank/node layout')
    p.add_argument('--slurm-container', type=Path, required=True)
    p.add_argument('--kaskade', type=Path, required=True)
    p.add_argument('--steps', type=int, default=10)
    p.add_argument('--output', type=Path)
    p.add_argument('--plan-only', action='store_true')
    x = p.parse_args(argv)
    repo = Path(__file__).resolve().parents[1]
    reference = x.reference.resolve()
    config = json.loads((reference/'environment.json').read_text())
    validation = json.loads((reference/'summary.json').read_text())
    require(not validation['errors'] and all(validation['runs'][n]['status'] == 'PASS' for n in CASES),
            'Reference must pass direct and all three FGMRES cases')
    a = accuracy_arguments(config, repo, x.slurm_container, x.kaskade, x.steps)
    print(f'Pepe: {a.expected_nodes} nodes, {a.ranks} ranks, {a.steps} steps, dt=0.01; FGMRES.', flush=True)
    print('16-bit transform-first; compressed cases split ON; histograms OFF; snapshots ON.', flush=True)
    print('Cases: '+', '.join(CASES), flush=True)
    if x.plan_only:
        for name in CASES:
            print(shlex.join(accuracy_command(a, Path('/planned-output')/name, name)))
        return 0
    require(bool(os.environ.get('SLURM_JOB_ID')) and
            int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) == a.expected_nodes,
            'Use an active allocation with exactly the reference node count')
    require(a.slurm_container.is_file(), 'Missing image')
    model = repo/'emiModel'
    model_hash = digest(model)
    require(model_hash == config['executable_sha256'], 'Executable changed; repeat one-step validation first')
    for filename in mesh_inputs(repo, 'pepe').values():
        require(Path(filename).is_file(), 'Missing mesh input: '+filename)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = (x.output or repo/'output'/f'emi-time-accuracy-pepe-n{a.expected_nodes}-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    for key in ('KASKADE_RESTRICTION_CAPTURE', 'KASKADE_FROZEN_TRANSFER_DIR',
                'KASKADE_TRANSFER_AUDIT_DIR', 'EMI_MATERIAL_TRACE'):
        env.pop(key, None)
    manifest = dict(configuration=vars(a), reference=str(reference),
        executable_sha256=model_hash, image_sha256=digest(a.slurm_container),
        reference_summary_sha256=digest(reference/'summary.json'),
        script_sha256=digest(Path(__file__)), dt=0.01,
        slurm={k: env.get(k) for k in ('SLURM_JOB_ID', 'SLURM_JOB_NODELIST', 'SLURM_JOB_NUM_NODES')},
        note='Short time-integration validation, not timing benchmark. Raw global-DOF errors; '
             'centered errors are supplementary. No membrane-voltage/Sobolev norm or hidden solution-error threshold. '
             'Steps 1,5,10 are checkpoints of the same trajectory. Direct uses one rank.')
    (root/'environment.json').write_text(json.dumps(manifest, default=str, indent=2)+'\n')
    report = dict(runs={}, errors=[], complete=False)

    def save():
        (root/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    print('Output: '+str(root), flush=True)
    for name in CASES:
        folder = root/name
        folder.mkdir()
        require(digest(model) == model_hash, 'Executable changed during campaign')
        cmd = accuracy_command(a, folder, name)
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print('Running '+name+' ...', flush=True)
        with (folder/'run.log').open('w') as out:
            rc = subprocess.run(cmd, cwd=repo, env=env, stdout=out, stderr=subprocess.STDOUT).returncode
        log = (folder/'run.log').read_text(errors='replace')
        row = dict(status='FAIL', exit_code=rc)
        report['runs'][name] = row
        try:
            if name == 'direct':
                row['relative_residuals'] = direct_check(log, rc, a.steps)
            else:
                row.update(collect(log, a.steps, a.ranks, rc, 1e-7))
                require(not row['errors'], '; '.join(row['errors']))
                row['routes'] = routing_counts(log, a.expected_nodes)
                require(sum(r['bytes'] for r in row['routes']) == row['transmitted_sent_bytes'],
                        'Routing byte sum differs from reported total')
                if name.endswith('-local-bypass'):
                    require('BDDC local compression bypass: on;' in log, 'Missing local-bypass confirmation')
            initial = read(folder/'potential-state-step0.csv')['u']
            require(initial == read(root/'direct/potential-state-step0.csv')['u'],
                    'Initial potential differs from direct reference')
            for step in range(1, a.steps+1):
                state = read(folder/f'potential-state-step{step}.csv')
                require(len(state['u']) == len(initial), 'Snapshot DOF count changed')
            snapshots_equal(reference/name, folder, (0, 1))
            row['one_step_reference_parity'] = 'EXACT_INITIAL_AND_FIRST_STEP_SNAPSHOTS'
            if name != 'direct':
                row['per_step_errors'] = trajectory_errors(folder, root/'direct',
                    root/'fgmres-uncompressed-destination', a.steps)
                if name != 'fgmres-uncompressed-destination':
                    require(read(folder/'potential-state-step1.csv')['rhs'] ==
                            read(root/'fgmres-uncompressed-destination/potential-state-step1.csv')['rhs'],
                            'First-step physical RHS differs from raw baseline')
                    baseline = report['runs']['fgmres-uncompressed-destination']
                    inter = lambda r: sum(v['bytes'] for v in r['routes'] if v['route'] == 'inter-node')
                    row['inter_node_saved_percent_vs_raw'] = 100*(1-inter(row)/inter(baseline))
            row['status'] = 'PASS'
            save()
            print('  PASS'+(' (errors recorded; no solution-error threshold)' if name != 'direct' else ''), flush=True)
        except (ValueError, OSError, KeyError) as exc:
            row['status'] = 'FAIL'
            report['errors'].append(name+': '+str(exc))
            save()
            print('  FAIL: '+str(exc), flush=True)
            return 1
    report['complete'] = True
    save()
    checkpoints = sorted({1, a.steps} | ({5} if a.steps >= 5 else set()))
    lines = ['# Pepe time-integration accuracy', '',
             'Errors are raw/centered relative Euclidean errors of global potential DOFs; no scientific threshold is imposed.', '',
             '| Case | Step | Iterations | True residual | Raw L2/direct | Centered L2/direct | Raw L2/uncompressed |',
             '|---|---:|---:|---:|---:|---:|---:|']
    for name in CASES[1:]:
        row = report['runs'][name]
        for step in checkpoints:
            e = row['per_step_errors'][step-1]
            lines.append(f'| {name} | {step} | {row["iterations_per_step"][step-1]} | '
                         f'{row["relative_residuals"][step-1]:.4e} | '
                         f'{e["vs_direct"]["relative_l2"]:.4e} | '
                         f'{e["vs_direct"]["centered_relative_l2"]:.4e} | '
                         f'{e["vs_uncompressed"]["relative_l2"]:.4e} |')
    (root/'accuracy.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)
    print('Results: '+str(root), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
