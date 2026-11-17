#!/usr/bin/env python3
"""Matched Pepe splitting checks; plans are safe locally, solves require SLURM."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_emi_partition_smoke as smoke
from analyze_emi_potential_diagnostics import read, error_metrics

pilot = smoke.pilot


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def cases(precision=16):
    return [dict(name=name, variant=variant, split=split, bits=bits, depth=1,
                 policy='P2', partition_mode='baseline', optimization=4)
            for name, variant, split, bits in [
                ('raw-fcg', 'raw', 'off', None),
                ('B-recenter-split-off', 'B-recenter', 'off', precision),
                ('B-split-on', 'B', 'on', precision),
                ('B-recenter-split-on', 'B-recenter', 'on', precision)]]


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--bits', type=int, choices=[12, 16], default=16)
    parser.add_argument('--executable', type=Path, default=Path('emiModel'))
    parser.add_argument('--kaskade', type=Path,
                        default=Path('/data/numerik/people/fchegini/kaskade7_new'))
    parser.add_argument('--slurm-container', type=Path, default=Path(
        '/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif'))
    parser.add_argument('--output', type=Path, default=Path('weight-aware-pepe-results'))
    return parser.parse_args(argv)


def configuration(a):
    cfg = smoke.arguments(['--executable', str(a.executable), '--kaskade', str(a.kaskade),
                           '--slurm-container', str(a.slurm_container),
                           '--output', str(a.output), '--byte-entropy'])
    cfg.steps = 3
    return cfg


def command(cfg, case, folder):
    settings, cmd = smoke.command(cfg, case, folder)
    compressed = case['bits'] is not None
    expected = dict(bddcUseCg=1, bddcFlexibleCg=1, bddcFgmres=0,
                    bddcCgOriginalResidual=1, bddcRestrictionPasses=1,
                    bddcCgDiagnostics=0, bddcConvergenceReport=1,
                    maximumNumberOfTimeSteps=3, bddcCompression=int(compressed),
                    bddcExcludeCornersFromCompression=1, bddcSerializedByteEntropy=1,
                    bddcGraphLifting=int(compressed), bddcGraphLiftingAfterQuantization=0,
                    bddcProlongateConstantSplit=int(case['split'] == 'on'),
                    bddcTransformFirstRecentering=int(case['variant'] == 'B-recenter'),
                    bddcHuffman=int(compressed), bddcBitlength=int(compressed),
                    bddcZigzag=int(compressed))
    if compressed:
        expected.update(bddcCompressionBits=case['bits'], bddcRestrictCompressionBits=case['bits'],
                        bddcProlongateCompressionBits=case['bits'], bddcGraphLiftingLevels=1,
                        bddcLosslessPacking=1)
    for flag, value in expected.items():
        require(cmd.count('--' + flag) == 1, 'Missing/duplicate option: ' + flag)
        require(cmd[cmd.index('--' + flag) + 1] == str(value), 'Wrong option: ' + flag)
    require(float(cmd[cmd.index('--bddcRelativeTolerance') + 1]) == 1e-7,
            'Tolerance differs')
    require(cmd[cmd.index('--bddcPartitionMode') + 1] == 'baseline', 'Ownership differs')
    return settings, cmd


def invariants(folder):
    problem = pilot.readcsv(folder / 'partition-problem.csv')[0]
    return dict(problem={k: v for k, v in problem.items() if k not in ['mode', 'planning_s']},
                owners=pilot.readcsv(folder / 'partition-ownership.csv'),
                ranks=pilot.readcsv(folder / 'partition-ranks.csv'),
                interfaces=pilot.readcsv(folder / 'partition-interfaces.csv'),
                coarse=digest(folder / 'partition-coarse.csv'),
                operator=digest(folder / 'potential-operator.csv'),
                initial=digest(folder / 'potential-state-step0.csv'))


def snapshot_residual(data):
    numerator = math.sqrt(math.fsum(x * x for x in data['residual']))
    denominator = math.sqrt(math.fsum(x * x for x in data['rhs']))
    return numerator / denominator if denominator else (0.0 if numerator == 0 else None)


def diagnose(folder, rc, cfg, case, reference):
    log = (folder / 'run.log').read_text(errors='replace')
    metrics = pilot.collect(log, cfg.steps, cfg.ranks, rc, 1e-7)
    require(not metrics['errors'], '; '.join(metrics['errors']))
    profiles, stages, totals = pilot.hb.profiles(folder, cfg.ranks, cfg.steps)
    transport, timing, completion = pilot.observations(folder, cfg.ranks, cfg.steps)
    pilot.validate_topology(folder, transport)
    entropy = pilot.serialized_entropy(folder, transport)
    smoke.check_placement(folder, case, cfg)
    require(not re.search(r'BDDC flexible CG breakdown:', log), 'FCG breakdown was reported')
    for record in completion:
        step = int(record['step'])
        require(record['status'] == 'PASS' and int(record['accepted_updates']) ==
                metrics['iterations_per_step'][step - 1], 'Rank completion disagrees')
    actual = sum(int(r['sent_bytes']) for r in transport if r['scope'] == 'inter_rank')
    require(actual == sum(r['packed_bytes'] for r in profiles), 'Byte counters disagree')
    if case['variant'] == 'raw':
        require(all(r['sent_bytes'] == r['raw_reference_bytes'] for r in transport),
                'Raw serialized bytes do not match raw-reference counters')
    per_step = []
    for step in range(1, cfg.steps + 1):
        data = read(folder / f'potential-state-step{step}.csv')
        residual = snapshot_residual(data)
        require(residual is not None and math.isfinite(residual) and residual <= 1e-7,
                f'Snapshot true residual exceeds tolerance at step {step}: {residual}')
        row = dict(configuration=case['name'], step=step, run_status='PASS',
                   accepted_iterations=metrics['iterations_per_step'][step - 1],
                   true_relative_residual=metrics['relative_residuals'][step - 1],
                   snapshot_true_relative_residual=residual,
                   relative_l2_vs_raw=None, max_absolute_vs_raw=None,
                   error_scope='RAW_REFERENCE' if reference is None else 'CONVERGED_STATE_VS_RAW')
        if reference is not None:
            error = error_metrics(data['u'], read(reference / f'potential-state-step{step}.csv')['u'])
            row.update(relative_l2_vs_raw=error['relative_l2'],
                       max_absolute_vs_raw=error['max_absolute_error'])
        for phase in ['restriction', 'prolongation']:
            samples = [r for r in transport if int(r['step']) == step and
                       r['phase'] == phase and r['scope'] == 'inter_rank']
            applications = {int(r['phase_applications']) for r in samples}
            require(len(applications) == 1, 'Ranks disagree about phase applications')
            row[phase + '_applications'] = applications.pop()
            row[phase + '_bytes'] = sum(int(r['sent_bytes']) for r in samples)
        per_step.append(row)
    summary = dict(iterations_per_step=metrics['iterations_per_step'],
                   accepted_iterations=sum(metrics['iterations_per_step']),
                   true_relative_residual_max=max(r['true_relative_residual'] for r in per_step),
                   relative_l2_vs_raw_max=max((r['relative_l2_vs_raw'] for r in per_step
                                               if r['relative_l2_vs_raw'] is not None), default=None),
                   max_absolute_vs_raw=max((r['max_absolute_vs_raw'] for r in per_step
                                             if r['max_absolute_vs_raw'] is not None), default=None),
                   serialized_inter_rank_bytes=actual,
                   serialized_inter_node_bytes=sum(int(r['sent_bytes']) for r in transport
                                                   if r['scope'] == 'inter_node'),
                   diagnostic_application_seconds=max(r['application_compute_s'] for r in totals),
                   diagnostic_fcg_seconds=max(sum(r['fcg_s'] for r in profiles if r['rank'] == rank)
                                               for rank in range(cfg.ranks)), **entropy)
    return summary, per_step


def main(argv=None):
    a = arguments(argv)
    cfg = configuration(a)
    configurations = cases(a.bits)
    commands = [command(cfg, case, a.output / case['name'])[1] for case in configurations]
    manifest = dict(model='EMI Pepe', nodes=2, ranks=4, timesteps=3, solver='FCG',
                    bit_depth=a.bits, graph_lifting_depth=1, ordering='transform-first',
                    corner_policy='exact', ownership='P2 baseline', configurations=configurations,
                    commands=commands, status='PREPARED_NOT_EXECUTED',
                    note='Diagnostic correctness runs, not a repeated performance benchmark. '
                         'Raw FCG agreement is not an independent exact-solution error. '
                         'No preassigned expected iteration count or solution-error threshold.')
    if not a.execute:
        print(json.dumps(manifest, indent=2, default=str))
        return 0
    require(os.environ.get('SLURM_JOB_ID') and
            (int(os.environ.get('SLURM_JOB_NUM_NODES', 0)),
             int(os.environ.get('SLURM_NTASKS', 0))) == (2, 4),
            'Execution requires a user-submitted 2-node/4-rank SLURM allocation')
    a.output = a.output.resolve()
    a.output.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[2]
    watched = [a.executable.resolve(), Path(__file__).resolve()]
    watched += [p for pattern in ['bddc*.hh', 'bddc*.hpp'] for p in (a.kaskade / 'mg').glob(pattern)]
    watched += [p for p in repo.glob('*.hh') if not p.name.startswith('.')]
    watched += [repo / 'EMI_model.cpp', repo / 'Makefile']
    watched += [p for p in (repo / 'tests').glob('*.py') if not p.name.startswith('.')]
    manifest.update(status='EXECUTING', job_id=os.environ['SLURM_JOB_ID'],
                    sha256={str(p): digest(p) for p in watched},
                    mesh_inputs_sha256={k: digest(v) for k, v in pilot.hb.mesh_inputs(repo, 'pepe').items()})
    (a.output / 'manifest.json').write_text(json.dumps(manifest, default=str, indent=2) + '\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith(
        ('KASKADE_', 'APPTAINERENV_KASKADE_', 'EMI_MESSAGE_ENTROPY_PREFIX'))}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    for key, value in dict(KASKADE_LIFTING_OPTIMIZATION='4', KASKADE_TIMESTEP_BYTES='1',
                           KASKADE_TRANSFER_TIMING='2', KASKADE_PACKET_TIMING='1').items():
        env[key] = value
        env['APPTAINERENV_' + key] = value
    runs, steps = [], []
    reference, fingerprint = None, None
    for case in configurations:
        folder = a.output / case['name']
        folder.mkdir()
        row = dict(configuration=case['name'], bits=case['bits'], run_status='FAILED', failure_reason=None,
                   folder=str(folder), solution_history_parity=None)
        runs.append(row)
        if case['variant'] != 'raw' and reference is None:
            row.update(run_status='SKIPPED', failure_reason='No validated matched raw FCG reference')
        else:
            settings, cmd = command(cfg, case, folder)
            (folder / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
            print('Running ' + case['name'], flush=True)
            try:
                rc, timeout, launcher = pilot.hb.study.execute(cmd, folder, settings, env)
                row.update(exit_code=rc, timed_out=timeout)
                require(not timeout, 'Application timeout')
                summary, case_steps = diagnose(folder, rc, cfg, case, reference)
                row.update(summary)
                current = invariants(folder)
                if reference is None:
                    reference, fingerprint = folder, current
                else:
                    require(current == fingerprint, 'Problem, initial state, ownership, topology or placement changed')
                if case['name'] == 'B-recenter-split-on':
                    peer = next(r for r in runs if r['configuration'] == 'B-split-on')
                    if peer['run_status'] == 'PASS':
                        exact = pilot.hb.exact_trajectory(Path(peer['folder']), folder, cfg.steps)
                        row['solution_history_parity'] = exact
                        require(exact, 'Reversible recentering changed states or residual histories')
                    else:
                        row['parity_status'] = 'UNAVAILABLE: split-on recenter-OFF reference failed'
                row['run_status'] = 'PASS'
                steps.extend(case_steps)
            except Exception as error:
                row['failure_reason'] = str(error)
                log = folder / 'run.log'
                if log.exists():
                    row['breakdown_lines'] = [line for line in log.read_text(errors='replace').splitlines()
                                              if 'breakdown' in line.lower() or 'BDDC global residual step' in line]
        (a.output / 'runs.json').write_text(json.dumps(runs, indent=2, allow_nan=False) + '\n')
        pilot.writecsv(a.output / 'summary.csv', runs)
        pilot.writecsv(a.output / 'per-step.csv', steps)
        print(json.dumps(row, indent=2, allow_nan=False), flush=True)
    require(all(digest(path) == value for path, value in manifest['sha256'].items()),
            'Sources or executable changed during the check')
    complete = all(r['run_status'] == 'PASS' for r in runs)
    manifest['status'] = 'COMPLETE' if complete else 'FAILED_OR_SKIPPED'
    (a.output / 'completion.json').write_text(json.dumps(dict(status=manifest['status']), indent=2) + '\n')
    (a.output / 'report.txt').write_text(
        'EMI Pepe weight-aware splitting check: 2 nodes, 4 ranks, 3 timesteps.\n'
        f'FCG, true tolerance 1e-7; compressed runs: {a.bits}-bit Method B, depth 1, full coding.\n'
        'Exact corners; unchanged P2 ownership and interface topology.\n'
        'See per-step.csv for convergence and global algebraic state error versus matched raw FCG.\n'
        'Byte entropy includes complete R/P application send buffers, pooled across ranks.\n'
        'Diagnostic timings have observer overhead; no speedup claim.\n' +
        '\n'.join(r['configuration'] + ': ' + r['run_status'] for r in runs) + '\n')
    print('Results: ' + str(a.output), flush=True)
    return 0 if complete else 1


if __name__ == '__main__':
    raise SystemExit(main())
