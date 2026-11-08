#!/usr/bin/env python3
"""Repeat validated destination-routing cases on Slurm; snapshots off, residual checks on."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import random
import shlex
import statistics
import subprocess
import time

from collect_emi_scaling_smoke import collect
from run_emi_mpi_routing_check import arguments as routing_arguments
from run_emi_mpi_routing_check import digest, make_command, routing_counts
from run_emi_histogram_start_check import require, mesh_inputs


KINDS = ('uncompressed-destination', 'compressed-split-destination',
         'compressed-split-local-bypass')


def schedule(solvers, warmups, repetitions, seed):
    rng = random.Random(seed)
    cases = [f'{solver}-{kind}' for solver in solvers for kind in KINDS]
    result = []
    for warmup, count in ((True, warmups), (False, repetitions)):
        for repeat in range(1, count+1):
            order = cases.copy()
            rng.shuffle(order)
            result.extend(dict(case=name, warmup=warmup, repeat=repeat) for name in order)
    return result


def timing_command(a, folder, name):
    solver = name.split('-')[0]
    cmd = make_command(a, folder, ['mpirun'], solver,
                       '-uncompressed-' not in name, True,
                       bypass=name.endswith('-local-bypass'))
    require(cmd.count('--potentialDiagnostics') == 1, 'Missing/duplicate snapshot option')
    cmd[cmd.index('--potentialDiagnostics')+1] = '0'
    return cmd


def route_totals(row):
    return {route: sum(r['bytes'] for r in row['routes'] if r['route'] == route)
            for route in ('same-rank', 'same-node', 'inter-node')}


def aggregate(rows):
    result = {}
    for name in sorted({r['case'] for r in rows}):
        measured = [r for r in rows if r['case'] == name and not r['warmup'] and r['status'] == 'PASS']
        if not measured:
            continue
        data = dict(samples=len(measured),
                    iterations=[r['iterations_per_step'][0] for r in measured],
                    true_residuals=[r['relative_residuals'][0] for r in measured])
        for metric in ('solver_wrapper_max_seconds', 'wall_seconds', 'launcher_seconds',
                       'transmitted_sent_bytes'):
            values = [r[metric] for r in measured]
            data[metric] = dict(median=statistics.median(values), min=min(values), max=max(values))
        data['route_bytes'] = {route: statistics.median([route_totals(r)[route] for r in measured])
                               for route in ('same-rank', 'same-node', 'inter-node')}
        result[name] = data
    for name, data in result.items():
        base = result.get(name.split('-')[0]+'-uncompressed-destination')
        if not base:
            continue
        for metric in ('solver_wrapper_max_seconds', 'wall_seconds'):
            data[metric]['change_percent_vs_raw'] = 100*(data[metric]['median']/base[metric]['median']-1)
        denominator = base['route_bytes']['inter-node']
        data['inter_node_saved_percent_vs_raw'] = (
            100*(1-data['route_bytes']['inter-node']/denominator) if denominator else None)
    return result


def check_reference(row, reference):
    # Snapshot observation should not alter these quantities. A different trajectory
    # needs review before attributing a timing difference to the observer being off.
    for key in ('cells', 'dofs', 'subdomains', 'coarse_interfaces', 'iterations_per_step',
                'raw_sent_bytes', 'transmitted_sent_bytes', 'routes'):
        require(row[key] == reference[key], 'Timing/reference mismatch: '+key)


def reference_arguments(config, repo, image, kaskade):
    ranks, nodes = config['ranks'], config['expected_nodes']
    require(config['case'] == 'pepe' and ranks >= nodes >= 2 and ranks % nodes == 0,
            'Require validated multi-node Pepe with equal ranks per node')
    require(config['slurm_btl'] == 'self,vader,tcp' and config['slurm_vader_single_copy'] == 'none',
            'Reference must use Vader with single-copy disabled')
    return routing_arguments(['--repo', str(repo), '--case', 'pepe', '--ranks', str(ranks),
        '--expected-nodes', str(nodes), '--slurm-container', str(image),
        '--kaskade', str(kaskade), '--slurm-btl', 'self,vader,tcp',
        '--slurm-vader-single-copy', 'none', '--check-local-bypass',
        '--bits', str(config['bits']), '--iterations', str(config['iterations']),
        '--transform-order', config['transform_order'], '--solvers', *config['solvers']])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', required=True, type=Path)
    p.add_argument('--slurm-container', required=True, type=Path)
    p.add_argument('--kaskade', required=True, type=Path)
    p.add_argument('--repetitions', type=int, default=3)
    p.add_argument('--warmups', type=int, default=1)
    p.add_argument('--seed', type=int, default=20261107)
    p.add_argument('--output', type=Path)
    p.add_argument('--plan-only', action='store_true')
    x = p.parse_args(argv)
    require(x.repetitions >= 2 and x.warmups >= 1, 'Require >=2 repetitions and >=1 warm-up per case')
    reference = x.reference.resolve()
    config = json.loads((reference/'environment.json').read_text())
    validation = json.loads((reference/'summary.json').read_text())
    require(not validation['errors'], 'Reference campaign has failures')
    repo = Path(__file__).resolve().parents[1]
    a = reference_arguments(config, repo, x.slurm_container, x.kaskade)
    jobs = schedule(a.solvers, x.warmups, x.repetitions, x.seed)
    for job in jobs:
        require(validation['runs'][job['case']]['status'] == 'PASS', 'Reference case did not pass')
    print(f'{len(jobs)} solves: {x.warmups} warm-up(s) + {x.repetitions} measured runs per case.', flush=True)
    print('Destination routing; snapshots/histograms OFF; residual checks/profiles ON.', flush=True)
    print(f'Validated layout: {a.ranks} ranks on {a.expected_nodes} nodes '
          f'({a.ranks//a.expected_nodes} ranks per node).', flush=True)
    print(f'Pepe: {a.bits} bits, transform={a.transform_order}; transport=self,vader,tcp; single-copy=none.', flush=True)
    if x.plan_only:
        for j in jobs:
            print(('warmup' if j['warmup'] else 'measured'), j['repeat'], j['case'])
        return 0
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Run in an active Slurm allocation')
    require(int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) == a.expected_nodes,
            f'Use an allocation with exactly {a.expected_nodes} nodes, matching the reference')
    require(a.slurm_container.is_file(), 'Container image is missing')
    model = repo/'emiModel'
    executable_hash = digest(model)
    require(executable_hash == config['executable_sha256'], 'Executable differs from validated reference; validate again first')
    for filename in mesh_inputs(repo, a.case).values():
        require(Path(filename).is_file(), 'Missing mesh input: '+filename)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = (x.output or repo/'output'/f'emi-routing-timings-pepe-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    for key in ('KASKADE_RESTRICTION_CAPTURE', 'KASKADE_FROZEN_TRANSFER_DIR',
                'KASKADE_TRANSFER_AUDIT_DIR', 'EMI_MATERIAL_TRACE'):
        env.pop(key, None)
    manifest = dict(configuration=vars(a), reference=str(reference),
        reference_summary_sha256=digest(reference/'summary.json'),
        executable_sha256=executable_hash, image_sha256=digest(a.slurm_container),
        script_sha256=digest(Path(__file__)), seed=x.seed, schedule=jobs,
        slurm={k: env.get(k) for k in ('SLURM_JOB_ID', 'SLURM_JOB_NODELIST', 'SLURM_JOB_NUM_NODES')},
        note='Warm-ups excluded. Fresh process/setup per solve. Residual checks and solver traces retained. '
             'Snapshot output disabled; solution errors belong to the validated reference, not these repeats. '
             'No assumed one-time codebook training or subtraction of codec cost.')
    (root/'environment.json').write_text(json.dumps(manifest, default=str, indent=2)+'\n')
    report = dict(runs=[], aggregate={}, complete=False)

    def save():
        report['aggregate'] = aggregate(report['runs'])
        (root/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    print('Output: '+str(root), flush=True)
    save()
    for index, job in enumerate(jobs, 1):
        require(digest(model) == executable_hash, 'Executable changed during campaign')
        folder = root/f'{index:03d}-{job["case"]}'
        folder.mkdir()
        cmd = timing_command(a, folder, job['case'])
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print(f'[{index}/{len(jobs)}] '+('warmup ' if job['warmup'] else 'measured ')+job['case'], flush=True)
        start = time.monotonic()
        with (folder/'run.log').open('w') as out:
            rc = subprocess.run(cmd, cwd=repo, env=env, stdout=out, stderr=subprocess.STDOUT).returncode
        elapsed = time.monotonic()-start
        log = (folder/'run.log').read_text(errors='replace')
        row = dict(collect(log, 1, a.ranks, rc, 1e-7), **job, launcher_seconds=elapsed,
                   folder=str(folder))
        row['measurement_kind'] = 'timing_with_residual_checks_and_profiles'
        report['runs'].append(row)
        try:
            require(not row['errors'], '; '.join(row['errors']))
            for metric in ('wall_seconds', 'solver_wrapper_max_seconds'):
                require(row[metric] > 0, 'Nonpositive timing: '+metric)
            row['routes'] = routing_counts(log, a.expected_nodes)
            require(sum(r['bytes'] for r in row['routes']) == row['transmitted_sent_bytes'], 'Route byte sum differs')
            if job['case'].endswith('-local-bypass'):
                require('BDDC local compression bypass: on;' in log, 'Local bypass was not enabled')
            check_reference(row, validation['runs'][job['case']])
            row['status'] = 'PASS'
        except ValueError as exc:
            row['status'] = 'FAIL'
            row['errors'].append(str(exc))
            save()
            print('FAIL: '+str(exc), flush=True)
            return 1
        save()
        print(f'  PASS: iterations={row["iterations_per_step"]}; solver={row["solver_wrapper_max_seconds"]:.3f}s; '
              f'application={row["wall_seconds"]:.3f}s', flush=True)
    report['complete'] = True
    save()
    lines = ['# Repeated Pepe timings', '',
        'Warm-ups excluded. Seconds are median [minimum, maximum]. These runs retain residual checks and profiles.', '',
        '| Case | n | Iterations | Solver wrapper s | Application s | Inter-node bytes (median) |',
        '|---|---:|---|---|---|---:|']
    for name, data in report['aggregate'].items():
        def interval(metric):
            v = data[metric]
            return f'{v["median"]:.3f} [{v["min"]:.3f}, {v["max"]:.3f}]'
        lines.append(f'| {name} | {data["samples"]} | {data["iterations"]} | '
                     f'{interval("solver_wrapper_max_seconds")} | {interval("wall_seconds")} | '
                     f'{data["route_bytes"]["inter-node"]:g} |')
    (root/'timings.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)
    print('Results: '+str(root), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
