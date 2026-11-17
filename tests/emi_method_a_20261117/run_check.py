#!/usr/bin/env python3
"""Pepe Method A depth/split study. Default is plan-only; no job submission."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BASE = REPO / 'tests/emi_weight_aware_20261117/run_check.py'
spec = importlib.util.spec_from_file_location('weight_aware_check', BASE)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
pilot = base.pilot
require = base.require


def cases(bits):
    result = [dict(name='raw-fcg', variant='raw', bits=None, depth=1, split='off',
                   policy='P2', partition_mode='baseline', optimization=0)]
    result += [dict(name=f'A-b{bits}-d{depth}-split-{split}', variant='A', bits=bits,
                    depth=depth, split=split, policy='P2', partition_mode='baseline', optimization=5)
               for split in ['off', 'on'] for depth in [1, 2, 3, 4]]
    return result


def command(cfg, case, folder):
    settings, cmd = base.smoke.command(cfg, case, folder)
    compressed = case['variant'] == 'A'
    expected = dict(bddcUseCg=1, bddcFlexibleCg=1, bddcFgmres=0,
                    bddcCgOriginalResidual=1, bddcRestrictionPasses=1,
                    bddcCgDiagnostics=0, bddcConvergenceReport=1,
                    maximumNumberOfTimeSteps=3, bddcCompression=int(compressed),
                    bddcExcludeCornersFromCompression=1, bddcSerializedByteEntropy=1,
                    bddcGraphLifting=0,  # This flag enables floating transform-first lifting.
                    bddcGraphLiftingAfterQuantization=int(compressed),
                    bddcProlongateConstantSplit=int(case['split'] == 'on'),
                    bddcTransformFirstRecentering=0,
                    bddcHuffman=int(compressed), bddcBitlength=int(compressed))
    if compressed:
        expected.update(bddcCompressionBits=case['bits'], bddcRestrictCompressionBits=case['bits'],
                        bddcProlongateCompressionBits=case['bits'],
                        bddcGraphLiftingLevels=case['depth'], bddcLosslessPacking=1)
    for flag, value in expected.items():
        require(cmd.count('--' + flag) == 1, 'Missing/duplicate option: ' + flag)
        require(cmd[cmd.index('--' + flag) + 1] == str(value), 'Wrong option: ' + flag)
    require(float(cmd[cmd.index('--bddcRelativeTolerance') + 1]) == 1e-7, 'Tolerance differs')
    require(cmd[cmd.index('--bddcPartitionMode') + 1] == 'baseline', 'Ownership differs')
    return settings, cmd


def main(argv=None):
    a = base.arguments(argv)
    a.output = a.output.resolve()
    cfg = base.configuration(a)
    configurations = cases(a.bits)
    manifest = dict(model='EMI Pepe', nodes=2, ranks=4, timesteps=3, solver='FCG',
                    bit_depth=a.bits, depths=[1, 2, 3, 4], splitting=['off', 'on'],
                    ordering='quantize-first', recentering='Not applicable to Method A',
                    corner_policy='exact', ownership='P2 baseline', configurations=configurations,
                    commands=[command(cfg, case, a.output / case['name'])[1] for case in configurations],
                    status='PREPARED_NOT_EXECUTED',
                    note='Diagnostic correctness study, not a repeated performance benchmark.')
    if not a.execute:
        print(json.dumps(manifest, indent=2))
        return 0
    require(os.environ.get('SLURM_JOB_ID') and
            (int(os.environ.get('SLURM_JOB_NUM_NODES', 0)), int(os.environ.get('SLURM_NTASKS', 0))) == (2, 4),
            'Execution requires a user-submitted 2-node/4-rank SLURM allocation')
    a.output.mkdir(parents=True, exist_ok=False)
    watched = [a.executable.resolve(), Path(__file__).resolve(), BASE, REPO / 'EMI_model.cpp', REPO / 'Makefile']
    watched += list(REPO.glob('*.hh')) + list((REPO / 'tests').glob('*.py'))
    watched += [p for pattern in ['bddc*.hh', 'bddc*.hpp'] for p in (a.kaskade / 'mg').glob(pattern)]
    manifest.update(status='EXECUTING', job_id=os.environ['SLURM_JOB_ID'],
                    sha256={str(p): base.digest(p) for p in watched},
                    mesh_inputs_sha256={k: base.digest(v) for k, v in pilot.hb.mesh_inputs(REPO, 'pepe').items()})
    (a.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith(
        ('KASKADE_', 'APPTAINERENV_KASKADE_', 'EMI_MESSAGE_ENTROPY_PREFIX'))}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    for key, value in dict(KASKADE_TIMESTEP_BYTES='1', KASKADE_TRANSFER_TIMING='2', KASKADE_PACKET_TIMING='1').items():
        env[key] = value
        env['APPTAINERENV_' + key] = value
    runs, steps, entropy_rows = [], [], []
    reference, fingerprint = None, None
    depth_references = {}
    for case in configurations:
        folder = a.output / case['name']
        folder.mkdir()
        row = dict(configuration=case['name'], method=case['variant'], bits=case['bits'],
                   lifting_depth=case['depth'] if case['bits'] else None, splitting=case['split'],
                   run_status='FAILED', failure_reason=None, folder=str(folder),
                   depth1_solution_history_parity=None)
        runs.append(row)
        if case['variant'] != 'raw' and reference is None:
            row.update(run_status='SKIPPED', failure_reason='No validated matched raw reference')
        else:
            settings, cmd = command(cfg, case, folder)
            runenv = dict(env, KASKADE_LIFTING_OPTIMIZATION=str(case['optimization']),
                          APPTAINERENV_KASKADE_LIFTING_OPTIMIZATION=str(case['optimization']))
            (folder / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
            print('Running ' + case['name'], flush=True)
            try:
                rc, timeout, launcher = pilot.hb.study.execute(cmd, folder, settings, runenv)
                row.update(exit_code=rc, timed_out=timeout)
                require(not timeout, 'Application timeout')
                summary, case_steps = base.diagnose(folder, rc, cfg, case, reference)
                row.update(summary)
                current = base.invariants(folder)
                if reference is None:
                    reference, fingerprint = folder, current
                else:
                    require(current == fingerprint, 'Problem, topology, ownership or initial state changed')
                if case['variant'] == 'A':
                    if case['depth'] == 1:
                        depth_references[case['split']] = folder
                    elif case['split'] in depth_references:
                        exact = pilot.hb.exact_trajectory(depth_references[case['split']], folder, cfg.steps)
                        row['depth1_solution_history_parity'] = exact
                        require(exact, 'Integer lifting depth changed state or residual history')
                    else:
                        row['parity_status'] = 'UNAVAILABLE: depth-one reference failed'
                raw = runs[0]
                for scope in ['inter_rank', 'inter_node']:
                    count = row['serialized_' + scope + '_bytes']
                    row[scope + '_saved_percent'] = 100 * (1 - count / raw['serialized_' + scope + '_bytes'])
                    row[scope + '_bytes_per_accepted_iteration'] = count / row['accepted_iterations']
                row['run_status'] = 'PASS'
                steps.extend(case_steps)
            except Exception as error:
                row['failure_reason'] = str(error)
                log = folder / 'run.log'
                if log.exists():
                    row['breakdown_lines'] = [line for line in log.read_text(errors='replace').splitlines()
                                              if 'breakdown' in line.lower() or 'BDDC global residual step' in line]
        entropy_file = folder / 'serialized-message-byte-entropy.csv'
        if entropy_file.exists():
            entropy_rows.extend(dict(r, configuration=case['name'], configuration_status=row['run_status'])
                                for r in pilot.readcsv(entropy_file))
        (a.output / 'runs.json').write_text(json.dumps(runs, indent=2, allow_nan=False) + '\n')
        pilot.writecsv(a.output / 'summary.csv', runs)
        pilot.writecsv(a.output / 'per-step.csv', steps)
        pilot.writecsv(a.output / 'serialized-message-byte-entropy.csv', entropy_rows)
        print(json.dumps(row, indent=2, allow_nan=False), flush=True)
    unchanged = all(base.digest(path) == value for path, value in manifest['sha256'].items())
    complete = unchanged and all(r['run_status'] == 'PASS' for r in runs)
    (a.output / 'completion.json').write_text(json.dumps(dict(status='COMPLETE' if complete else 'FAILED_OR_SKIPPED',
                                                            sources_unchanged=unchanged), indent=2) + '\n')
    (a.output / 'report.txt').write_text(
        f'EMI Pepe Method A: {a.bits} bits, depths 1-4, splitting OFF/ON.\n'
        'Quantize first, reversible integer lifting, full coding, exact corners; no Method A recentering switch.\n'
        'Matched raw FCG reference; P2, 2 nodes, 4 ranks, 3 timesteps; true tolerance 1e-7.\n'
        'Entropy: full R/P serialized send buffers, pooled across ranks, separate scopes.\n'
        'Diagnostic observer timings; no repeated-performance or speedup claim.\n'
        + '\n'.join(r['configuration'] + ': ' + r['run_status'] for r in runs) + '\n')
    print('Results: ' + str(a.output), flush=True)
    return 0 if complete else 1


if __name__ == '__main__':
    raise SystemExit(main())
