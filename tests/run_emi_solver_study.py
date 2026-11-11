#!/usr/bin/env python3
"""One-step solver/precision screen: selectable transform order, all-interface compression.

No scientific solution-error threshold. Validation timings include snapshots.
Compressed failures are retained; invalid reference solves stop the campaign.
"""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import time

from analyze_emi_potential_diagnostics import read, error_metrics
from collect_emi_scaling_smoke import collect, NUMBER
from run_emi_histogram_start_check import command, mesh_inputs, require, LEGO_WEAK_CASES
from run_emi_mpi_routing_check import digest, routing_counts, slurm_container_command


def arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--executable', type=Path, help='Separate executable; default: REPO/emiModel')
    p.add_argument('--case', choices=['pepe', 'rizzo', 'robin', 'kermit', 'gonzo', 'animal', '10cells', *LEGO_WEAK_CASES], default='pepe')
    p.add_argument('--ranks', type=int, default=16)
    p.add_argument('--master-worker', type=int, choices=[0, 1], default=0,
                   help='Reserve rank 0 as coarse master; CG/FCG only')
    p.add_argument('--expected-nodes', type=int, default=8)
    p.add_argument('--slurm-container', type=Path)
    p.add_argument('--kaskade', type=Path, default=Path(os.environ.get('KASKADE7', '/kaskade7')))
    p.add_argument('--bits', type=int, nargs='+', choices=[8, 12, 16, 24, 32],
                   default=[8, 12, 16, 24, 32])
    p.add_argument('--solvers', nargs='+', choices=['cg', 'fcg', 'fgmres'],
                   default=['cg', 'fcg', 'fgmres'])
    p.add_argument('--split', choices=['off', 'on', 'both'], default='both',
                   help='Compressed cases only; uncompressed references always use split off')
    p.add_argument('--transform-order', choices=['before', 'after', 'both'], default='before',
                   help='One order per campaign; before remains the default')
    p.add_argument('--iterations', type=int, default=1000)
    p.add_argument('--case-minutes', type=int, default=20,
                   help='Bound each solve; timeout is reported separately from nonconvergence')
    p.add_argument('--plan-only', action='store_true')
    p.add_argument('--setup-diagnostic-only', action='store_true',
                   help='One selected uncompressed solver, with material progress; no direct/error comparison')
    p.add_argument('--output', type=Path)
    a = p.parse_args(argv)
    if min(a.ranks, a.expected_nodes, a.iterations, a.case_minutes) < 1 or a.ranks % a.expected_nodes:
        p.error('Positive values and equal ranks per node are required')
    if len(set(a.solvers)) != len(a.solvers) or len(set(a.bits)) != len(a.bits):
        p.error('Duplicate solver or precision')
    if a.setup_diagnostic_only and len(a.solvers) != 1:
        p.error('Setup diagnostic requires exactly one --solvers choice')
    if a.master_worker and (a.ranks < 2 or 'fgmres' in a.solvers):
        p.error('Master-worker requires at least two ranks and --solvers cg and/or fcg; FGMRES is unsupported')
    if a.master_worker and a.split != 'off' and not a.setup_diagnostic_only:
        p.error('Master-worker with CG/FCG constant splitting is unsupported; use --split off')
    a.repo, a.kaskade = a.repo.resolve(), a.kaskade.resolve()
    if a.executable:
        a.executable = a.executable.resolve()
    if a.slurm_container:
        a.slurm_container = a.slurm_container.resolve()
    a.steps, a.uncompressed_steps, a.histogram_steps = 1, 0, []
    a.slurm_btl, a.slurm_vader_single_copy = 'self,vader,tcp', 'none'
    return a


def planned_cases(a):
    if a.setup_diagnostic_only:
        return [dict(name=f'{a.solvers[0]}-uncompressed', solver=a.solvers[0],
                     bits=None, split='off', direct=False)]
    cases = [dict(name='direct', solver='fgmres', bits=None, split='off', direct=True)]
    cases += [dict(name=f'{s}-uncompressed', solver=s, bits=None, split='off', direct=False)
              for s in a.solvers]
    cases += [dict(name=f'{s}-{a.transform_order}-{b}-split-{split}', solver=s, bits=b,
                   split=split, direct=False)
              for b in a.bits for s in a.solvers
              for split in (('off', 'on') if a.split == 'both' else (a.split,))]
    return cases


def make_command(a, case, folder, launcher):
    cfg = copy.copy(a)
    cfg.bits, cfg.solver, cfg.split = case['bits'] or 16, case['solver'], case['split']
    compressed = case['bits'] is not None
    mode = 'direct' if case['direct'] else 'scheduled' if compressed else 'uncompressed'
    cmd = command(cfg, mode, folder, launcher)
    cmd[cmd.index('--masterWorker')+1] = str(0 if case['direct'] else getattr(a, 'master_worker', 0))
    # The opt-in gate applies only to split CG/FCG, never to the unsplit cases.
    cmd[cmd.index('--bddcCgConstantSplitExperiment')+1] = str(int(
        compressed and cfg.split == 'on' and cfg.solver in ('cg', 'fcg')))
    cmd += ['--bddcMpiRouting', str(int(not case['direct'])), '--bddcLocalCompressionBypass', '0']
    if a.setup_diagnostic_only and not a.slurm_container:
        model_index = cmd.index(str(a.repo/'emiModel'))
        cmd[model_index:model_index] = ['env', 'EMI_MATERIAL_TRACE=1']
    if a.slurm_container:
        payload = cmd[cmd.index(str(a.repo/'emiModel')):]
        nodes, ranks = (1, 1) if case['direct'] else (a.expected_nodes, a.ranks)
        cmd = slurm_container_command(a, payload, nodes, ranks, mpi='pmix_v5')
        cmd[cmd.index('--time=00:10:00')] = f'--time={a.case_minutes}'
        cmd.insert(1, '--cpu-bind=cores')
    if getattr(a, 'executable', None):
        cmd[cmd.index(str(a.repo/'emiModel'))] = str(a.executable)
    return cmd


def execute(cmd, folder, a, env):
    """Bound the launcher and its children; Slurm also enforces the step limit."""
    started = time.monotonic()
    with (folder/'run.log').open('w') as stream:
        process = subprocess.Popen(cmd, cwd=a.repo, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        timed_out = False
        try:
            rc = process.wait(timeout=60*a.case_minutes + (60 if a.slurm_container else 0))
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                rc = process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                rc = process.wait()
    return rc, timed_out, time.monotonic()-started


def failure_kind(log, timed_out):
    if timed_out or 'DUE TO TIME LIMIT' in log:
        return 'runtime limit'
    if 'local space omits a nonzero material coupling endpoint' in log:
        return 'material local-space support mismatch'
    for reason in ('iteration limit', 'solver breakdown', 'solver stagnation'):
        if reason in log.lower():
            return reason
    return 'invalid or incomplete run'


def trace_metrics(folder, solver):
    path = folder/('bddc-fgmres-step1.csv' if solver == 'fgmres' else 'bddc-convergence-step1.csv')
    if not path.exists():
        return {}
    with path.open() as stream:
        samples = list(csv.DictReader(stream))
    if not samples:
        return {}
    last = samples[-1]
    relative = float(last.get('true_relative_residual', last.get('relative_residual')))
    return dict(trace_accepted_updates=int(last['iteration']),
                trace_relative_residual=relative if math.isfinite(relative) else None)


def analyze(a, case, folder, rc, timed_out, references, baseline):
    log = (folder/'run.log').read_text(errors='replace')
    row = dict(case, status='FAIL', exit_code=rc, errors=[], reason=None)
    if case['direct']:
        values = re.findall(rf'^Direct global residual step 1: abs=({NUMBER}) rel=({NUMBER})', log, re.M)
        if (rc or timed_out or 'EMI model completed' not in log or len(values) != 1 or
                not math.isfinite(float(values[0][1])) or not 0 <= float(values[0][1]) <= 1e-7):
            row['errors'].append('Direct reference failed residual/completion check')
        else:
            row['relative_residuals'] = [float(values[0][1])]
    else:
        row.update(collect(log, 1, a.ranks, rc, 1e-7))
        try:
            row.update(trace_metrics(folder, case['solver']))
        except (OSError, ValueError, KeyError, TypeError) as exc:
            row['errors'].append('Invalid convergence trace: '+str(exc))
        # Complete accounting is required for PASS; absent failed-run counters stay absent.
        try:
            row['routes'] = routing_counts(log, a.expected_nodes)
            require(sum(r['bytes'] for r in row['routes']) == row['transmitted_sent_bytes'],
                    'Routing counters do not sum to total payload')
            for route in ('same-rank', 'same-node', 'inter-node'):
                row[route+'_bytes'] = sum(r['bytes'] for r in row['routes'] if r['route'] == route)
        except ValueError as exc:
            row['errors'].append(str(exc))
    if timed_out:
        row['errors'].append('Launcher runtime limit')
    # Only complete, converged states are used as solution-accuracy evidence.
    if not row['errors']:
        try:
            initial = read(folder/'potential-state-step0.csv')['u']
            final = read(folder/'potential-state-step1.csv')['u']
            if not case['direct'] and not a.setup_diagnostic_only:
                require(initial == references['initial'], 'Initial state differs from direct reference')
                row['error_vs_direct'] = error_metrics(final, references['direct'])
                if case['bits'] is not None:
                    row['error_vs_uncompressed'] = error_metrics(final, references[case['solver']])
                else:
                    row['error_vs_uncompressed'] = error_metrics(final, final)
                control = baseline if case['bits'] is not None else row
                for key, label in [('transmitted_sent_bytes', 'logical'), ('inter-node_bytes', 'inter_node')]:
                    denominator = control.get(key)
                    row[label+'_bytes_saved_percent'] = (
                        100*(1-row[key]/denominator) if denominator else None)
            if case['direct']:
                references['initial'], references['direct'] = initial, final
            elif case['bits'] is None:
                references[case['solver']] = final
        except (OSError, KeyError, ValueError) as exc:
            row['errors'].append(str(exc))
    row['status'] = 'FAIL' if row['errors'] else 'PASS'
    if row['errors']:
        row['reason'] = failure_kind(log, timed_out)
    return row


def flat_row(row):
    def last(key):
        return row.get(key, [None])[-1] if row.get(key) else None
    error = row.get('error_vs_direct', {})
    raw_error = row.get('error_vs_uncompressed', {})
    return dict(case=row['name'], status=row['status'], reason=row.get('reason'),
        iterations=last('iterations_per_step') if row.get('iterations_per_step') else row.get('trace_accepted_updates'),
        true_relative_residual=last('relative_residuals') if row.get('relative_residuals') else row.get('trace_relative_residual'),
        raw_relative_l2_direct=error.get('relative_l2'),
        centered_relative_l2_direct=error.get('centered_relative_l2'),
        max_absolute_error_direct=error.get('max_absolute_error'),
        raw_relative_l2_uncompressed=raw_error.get('relative_l2'),
        centered_relative_l2_uncompressed=raw_error.get('centered_relative_l2'),
        max_absolute_error_uncompressed=raw_error.get('max_absolute_error'),
        same_rank_bytes=row.get('same-rank_bytes'), same_node_bytes=row.get('same-node_bytes'),
        inter_node_bytes=row.get('inter-node_bytes'), logical_bytes=row.get('transmitted_sent_bytes'),
        logical_bytes_saved_percent=row.get('logical_bytes_saved_percent'),
        inter_node_bytes_saved_percent=row.get('inter_node_bytes_saved_percent'),
        solver_wrapper_seconds=row.get('solver_wrapper_max_seconds'), application_seconds=row.get('wall_seconds'),
        launcher_seconds=row.get('launcher_seconds'))


def save(root, report):
    (root/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    rows = [flat_row(r) for r in report['runs'].values()]
    if not rows:
        return
    with (root/'results.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fields = ['case', 'status', 'iterations', 'raw_relative_l2_direct',
              'logical_bytes_saved_percent', 'inter_node_bytes_saved_percent', 'solver_wrapper_seconds']
    lines = [f"# EMI solver and precision study: {report['mesh']} ({report['transform_order']})", '',
        'Validation timings include snapshots. Errors have no scientific acceptance threshold.',
        'Savings use each solver\'s complete uncompressed solve. Logical R/P payloads only.',
        'Failed-run iterations may be accepted updates from the final available trace; see summary.json.', '',
        '| '+' | '.join(fields)+' |', '| '+' | '.join(['---']*len(fields))+' |']
    if report.get('setup_diagnostic_only'):
        lines[:6] = [f"# EMI setup diagnostic: {report['mesh']}", '',
            'Uncompressed solve only. No direct-reference, accuracy or compression comparison.', '',
            'Timings include diagnostic progress and snapshots.', '']
    for row in rows:
        lines.append('| '+' | '.join('--' if row[k] is None else
            f'{row[k]:.6g}' if isinstance(row[k], float) else str(row[k]) for k in fields)+' |')
    (root/'summary.md').write_text('\n'.join(lines)+'\n')


def campaign(a):
    cases = planned_cases(a)
    print(f'Master-worker mode: {getattr(a, "master_worker", 0)}; direct reference remains serial.', flush=True)
    if a.setup_diagnostic_only:
        print('Setup diagnostic: one uncompressed solve; material progress ON; no direct/error comparison.', flush=True)
        print(f'{a.case}: {a.expected_nodes} nodes, {a.ranks} ranks; solver={a.solvers[0]}; compression/split OFF.', flush=True)
        print('One step; destination routing; snapshots ON; histograms/VTK OFF.', flush=True)
    else:
        print(f'{len(cases)} cases: direct + matched raw controls + compressed split={a.split}.', flush=True)
        print(f'{a.case}: {a.expected_nodes} nodes, {a.ranks} ranks; bits={a.bits}; solvers={a.solvers}.', flush=True)
        print(f'One step; transform order={a.transform_order}; destination routing; local bypass OFF; snapshots ON; histograms/VTK OFF.', flush=True)
    print('Planned runs: '+', '.join(c['name'] for c in cases), flush=True)
    if a.plan_only:
        return 0
    if a.slurm_container:
        require(bool(os.environ.get('SLURM_JOB_ID')), 'Run inside a Slurm allocation or use the batch script')
        require(a.slurm_container.is_file(), 'Missing container image')
        require(int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) >= a.expected_nodes, 'Insufficient allocated nodes')
        require(int(os.environ.get('SLURM_NTASKS', '0')) >= a.ranks, 'Insufficient allocated tasks')
    model = a.executable or a.repo/'emiModel'
    require(model.is_file() and os.access(model, os.X_OK), 'Build emiModel first')
    inputs = [Path(p) for p in mesh_inputs(a.repo, a.case).values()]
    for path in inputs:
        require(path.is_file(), 'Missing mesh input: '+str(path))
    model_hash = digest(model)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    kind = 'setup-diagnostic' if a.setup_diagnostic_only else 'solver-study'
    root = (a.output or a.repo/'output'/f'emi-{kind}-{a.case}-{a.transform_order}-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    for key in ('KASKADE_RESTRICTION_CAPTURE', 'KASKADE_FROZEN_TRANSFER_DIR', 'KASKADE_TRANSFER_AUDIT_DIR', 'EMI_MATERIAL_TRACE'):
        env.pop(key, None)
    sources = [a.repo/p for p in ('EMI_model.cpp', 'EMI_bddc.hh', 'EMI_local_support_diagnostic.hh', 'EMI_bddc_precision.hh',
        'tests/run_emi_solver_study.py', 'tests/run_emi_mpi_routing_check.py',
        'tests/run_emi_histogram_start_check.py', 'tests/collect_emi_scaling_smoke.py',
        'tests/analyze_emi_potential_diagnostics.py')]
    sources += [a.kaskade/'mg'/p for p in ('bddc.hh', 'bddc.hpp', 'bddcMpiRouting.hh',
                                          'bddcSpaceTransfer.hh', 'bddcSpaceTransfer.hpp')]
    coarse_header = a.kaskade/'mg/bddcCoarseReduction.hh'
    if coarse_header.is_file():
        sources.append(coarse_header)
    manifest = dict(vars(a), planned_cases=cases, executable_sha256=model_hash,
        input_sha256={str(p): digest(p) for p in inputs},
        source_sha256={str(p): digest(p) for p in sources},
        slurm={k: env.get(k) for k in ('SLURM_JOB_ID', 'SLURM_JOB_NODELIST', 'SLURM_JOB_NUM_NODES', 'SLURM_NTASKS')},
        note=('Setup diagnostic only; no solution-error comparison.' if a.setup_diagnostic_only else
              'One-step validation; not repeated performance measurements; no scientific error threshold.'))
    (root/'environment.json').write_text(json.dumps(manifest, indent=2, default=str)+'\n')
    print('Output: '+str(root), flush=True)
    if a.slurm_container:
        probe = slurm_container_command(a, ['/bin/bash', '-c',
            'test -r /kaskade7/mg/bddcMpiRouting.hh && test -x "$1" && '
            'printf "Node: %s\\n" "$SLURMD_NODENAME" && lscpu', 'check', str(model)],
            a.expected_nodes, a.expected_nodes)
        (root/'preflight-command.sh').write_text(shlex.join(probe)+'\n')
        with (root/'preflight.log').open('w') as stream:
            result = subprocess.run(probe, cwd=a.repo, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=660)
        require(result.returncode == 0, 'Container/CPU inventory preflight failed')
    launcher = ['mpirun'] if a.slurm_container else shlex.split(env.get('MPIEXEC', 'mpirun'))
    require(bool(launcher), 'Empty MPIEXEC')
    report = dict(mesh=a.case, transform_order=a.transform_order,
                  setup_diagnostic_only=a.setup_diagnostic_only,
                  complete=False, runs={}, failed_cases=[], fatal_error=None)
    references = {}
    for case in cases:
        require(digest(model) == model_hash, 'Executable changed during campaign')
        folder = root/case['name']
        folder.mkdir()
        cmd = make_command(a, case, folder, launcher)
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print('Running '+case['name']+' ...', flush=True)
        rc, timed_out, elapsed = execute(cmd, folder, a, env)
        baseline = report['runs'].get(f"{case['solver']}-uncompressed")
        row = analyze(a, case, folder, rc, timed_out, references, baseline)
        row['launcher_seconds'] = elapsed
        report['runs'][case['name']] = row
        if row['status'] != 'PASS':
            report['failed_cases'].append(case['name'])
            print('  FAIL: '+row['reason']+'; '+('; '.join(row['errors'])), flush=True)
            if row['reason'] == 'material local-space support mismatch':
                for line in (folder/'run.log').read_text(errors='replace').splitlines():
                    if 'local space omits a nonzero material coupling endpoint' in line:
                        print('  '+line.strip(), flush=True)
                        break
            if case['bits'] is None:
                report['fatal_error'] = ('Setup diagnostic failed: ' if a.setup_diagnostic_only else
                                         'Invalid reference: ')+case['name']
                save(root, report)
                return 2
        else:
            print('  PASS (uncompressed solve; no reference comparison)' if a.setup_diagnostic_only else
                  '  PASS (convergence and diagnostics; no solution-error threshold)', flush=True)
        save(root, report)
    report['complete'] = True
    save(root, report)
    print((root/'summary.md').read_text(), flush=True)
    print(f"Campaign complete: {len(report['failed_cases'])} failed configurations retained. Results: {root}", flush=True)
    return 1 if report['failed_cases'] else 0


if __name__ == '__main__':
    raise SystemExit(campaign(arguments()))
