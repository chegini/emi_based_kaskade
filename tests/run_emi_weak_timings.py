#!/usr/bin/env python3
"""Repeat validated LEGO FCG cases; snapshots off, residual checks and profiles on."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import random
import shlex
import statistics
import subprocess

import run_emi_solver_study as study
from run_emi_histogram_start_check import LEGO_WEAK_CASES, mesh_inputs, require
from run_emi_routing_timings import check_reference, route_totals
from run_emi_mpi_routing_check import digest, routing_counts
from collect_emi_scaling_smoke import collect


def reference_arguments(config, validation, repo, image, kaskade):
    require(config['case'] in LEGO_WEAK_CASES, 'Require a standard LEGO weak mesh')
    require(validation.get('complete') and not validation.get('failed_cases') and
            not validation.get('fatal_error') and not validation.get('setup_diagnostic_only'),
            'Require a complete successful validation campaign')
    require(config['solvers'] == ['fcg'] and config['bits'] in ([12, 16], [12, 24]) and
            config['split'] == 'on' and config['transform_order'] == 'before',
            'Require validated FCG + split, before, 12/16 or 12/24 bits')
    require(config['slurm_btl'] == 'self,vader,tcp' and
            config['slurm_vader_single_copy'] == 'none', 'Transport differs from study')
    require(config['ranks'] == 2*config['expected_nodes'], 'Require two ranks per node')
    a = study.arguments(['--repo', str(repo), '--case', config['case'],
        '--ranks', str(config['ranks']), '--expected-nodes', str(config['expected_nodes']),
        '--solvers', 'fcg', '--bits', *map(str, config['bits']), '--split', 'on',
        '--transform-order', 'before', '--iterations', str(config['iterations']),
        '--slurm-container', str(image), '--kaskade', str(kaskade), '--case-minutes', '10'])
    if config.get('executable'):
        a.executable = Path(config['executable']).resolve()
    cases = [c for c in study.planned_cases(a) if not c['direct']]
    for name in ['direct'] + [c['name'] for c in cases]:
        require(validation['runs'][name]['status'] == 'PASS', 'Failed reference: '+name)
    return a, cases


def schedule(cases, warmups, repetitions, seed):
    rng = random.Random(seed)
    result = []
    for warmup, count in ((True, warmups), (False, repetitions)):
        for repeat in range(1, count+1):
            order = list(cases)
            rng.shuffle(order)
            result.extend(dict(case=c['name'], warmup=warmup, repeat=repeat) for c in order)
    return result


def timing_command(a, case, folder):
    cmd = study.make_command(a, case, folder, ['mpirun'])
    require(cmd.count('--potentialDiagnostics') == 1, 'Snapshot flag missing/duplicated')
    cmd[cmd.index('--potentialDiagnostics')+1] = '0'
    return cmd


def aggregate(rows):
    result = {}
    for name in sorted({r['case'] for r in rows}):
        measured = [r for r in rows if r['case'] == name and not r['warmup'] and r['status'] == 'PASS']
        if not measured:
            continue
        data = dict(samples=len(measured), iterations=[r['iterations_per_step'][0] for r in measured])
        for metric in ('solver_wrapper_max_seconds', 'wall_seconds', 'launcher_seconds'):
            values = [r[metric] for r in measured]
            data[metric] = dict(median=statistics.median(values), min=min(values), max=max(values))
        data['route_bytes'] = {route: statistics.median([route_totals(r)[route] for r in measured])
                               for route in ('same-rank', 'same-node', 'inter-node')}
        result[name] = data
    baseline = result.get('fcg-uncompressed')
    if baseline:
        for data in result.values():
            for metric in ('solver_wrapper_max_seconds', 'wall_seconds'):
                data[metric]['change_percent_vs_raw'] = 100*(data[metric]['median']/baseline[metric]['median']-1)
            denominator = baseline['route_bytes']['inter-node']
            data['inter_node_saved_percent_vs_raw'] = 100*(1-data['route_bytes']['inter-node']/denominator) if denominator else None
    return result


def check_cpu_inventory(text, nodes, expected):
    rows = [json.loads(line) for line in text.splitlines() if line.startswith('{')]
    require(len(rows) == nodes and len({r['host'] for r in rows}) == nodes,
            'CPU inventory must include exactly one record per allocated node')
    require(all(r['cpu'] == expected for r in rows), 'Allocated CPU model differs from requested model')
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--slurm-container', type=Path, required=True)
    p.add_argument('--kaskade', type=Path, required=True)
    p.add_argument('--cpu-model', required=True)
    p.add_argument('--warmups', type=int, default=1)
    p.add_argument('--repetitions', type=int, default=3)
    p.add_argument('--seed', type=int, default=20261108)
    p.add_argument('--plan-only', action='store_true')
    x = p.parse_args(argv)
    require(x.warmups >= 1 and x.repetitions >= 3, 'Require >=1 warm-up and >=3 measured repeats')
    reference = x.reference.resolve()
    config = json.loads((reference/'environment.json').read_text())
    validation = json.loads((reference/'summary.json').read_text())
    repo = Path(__file__).resolve().parents[1]
    a, cases = reference_arguments(config, validation, repo, x.slurm_container, x.kaskade)
    jobs = schedule(cases, x.warmups, x.repetitions, x.seed)
    print(f'{a.case}: {a.expected_nodes} nodes, {a.ranks} ranks; {len(jobs)} solves.', flush=True)
    print('Snapshots/histograms/VTK OFF; residual checks/profiles ON; local bypass OFF.', flush=True)
    if x.plan_only:
        for job in jobs:
            print(job)
        return 0
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Run inside a Slurm allocation')
    allocated = int(os.environ.get('SLURM_JOB_NUM_NODES', '0'))
    require(allocated >= a.expected_nodes and int(os.environ.get('SLURM_NTASKS', '0')) >= a.ranks,
            'Insufficient allocated nodes/tasks')
    model = a.executable or repo/'emiModel'
    require(digest(model) == config['executable_sha256'], 'Executable changed: revalidate first')
    require(a.slurm_container.is_file(), 'Container image missing')
    for path in mesh_inputs(repo, a.case).values():
        require(digest(Path(path)) == config['input_sha256'].get(path), 'Mesh input changed: '+path)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = repo/'output'/f'emi-weak-timings-{a.case}-{stamp}-{os.getpid()}'
    root.mkdir(parents=True, exist_ok=False)
    print('Output: '+str(root), flush=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', LC_ALL='C')
    for key in ('KASKADE_RESTRICTION_CAPTURE', 'KASKADE_FROZEN_TRANSFER_DIR', 'KASKADE_TRANSFER_AUDIT_DIR', 'EMI_MATERIAL_TRACE'):
        env.pop(key, None)
    probe_code = ('import json,socket; from pathlib import Path; '
        'cpu=next(s.split(":",1)[1].strip() for s in Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")); '
        'print(json.dumps(dict(host=socket.gethostname(),cpu=cpu)))')
    probe = ['srun', '--mpi=none', '--nodes='+str(allocated), '--ntasks='+str(allocated),
             '--ntasks-per-node=1', '--cpus-per-task=1', '--time=00:02:00',
             'python3', '-c', probe_code]
    inventory = subprocess.run(probe, cwd=repo, env=env, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, timeout=180)
    (root/'cpu-inventory.log').write_text(inventory.stdout)
    require(inventory.returncode == 0, 'CPU probe failed; see cpu-inventory.log')
    cpus = check_cpu_inventory(inventory.stdout, allocated, x.cpu_model)
    manifest = dict(configuration=vars(a), reference=str(reference), schedule=jobs,
        cpu_inventory=cpus, cpu_model=x.cpu_model, seed=x.seed,
        reference_summary_sha256=digest(reference/'summary.json'),
        executable_sha256=digest(model), image_sha256=digest(a.slurm_container),
        source_sha256={str(repo/f): digest(repo/f) for f in (
            'tests/run_emi_weak_timings.py', 'tests/run_emi_solver_study.py',
            'tests/run_emi_histogram_start_check.py', 'tests/run_emi_mpi_routing_check.py',
            'tests/run_emi_routing_timings.py', 'tests/collect_emi_scaling_smoke.py')},
        slurm={k: env.get(k) for k in ('SLURM_JOB_ID', 'SLURM_JOB_NODELIST', 'SLURM_JOB_NUM_NODES', 'SLURM_NTASKS')},
        note='Fresh setup per solve; errors belong to validation references. Warm-ups excluded. No new solution-error threshold.')
    (root/'environment.json').write_text(json.dumps(manifest, indent=2, default=str)+'\n')
    report = dict(complete=False, runs=[], aggregate={})
    def save():
        report['aggregate'] = aggregate(report['runs'])
        (root/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    save()
    lookup = {c['name']: c for c in cases}
    for index, job in enumerate(jobs, 1):
        require(digest(model) == config['executable_sha256'], 'Executable changed during campaign')
        folder = root/f'{index:03d}-{job["case"]}'
        folder.mkdir()
        cmd = timing_command(a, lookup[job['case']], folder)
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print(f'[{index}/{len(jobs)}] '+('warmup ' if job['warmup'] else 'measured ')+job['case'], flush=True)
        rc, timeout, elapsed = study.execute(cmd, folder, a, env)
        log = (folder/'run.log').read_text(errors='replace')
        row = dict(collect(log, 1, a.ranks, rc, 1e-7), **job, launcher_seconds=elapsed, folder=str(folder))
        report['runs'].append(row)
        try:
            require(not timeout and not row['errors'], 'Timeout or invalid solve: '+'; '.join(row['errors']))
            for metric in ('wall_seconds', 'solver_wrapper_max_seconds'):
                require(row[metric] > 0, 'Nonpositive timing: '+metric)
            row['routes'] = routing_counts(log, a.expected_nodes)
            require(sum(r['bytes'] for r in row['routes']) == row['transmitted_sent_bytes'], 'Route byte sum differs')
            check_reference(row, validation['runs'][job['case']])
            row['status'] = 'PASS'
        except (ValueError, KeyError) as exc:
            row['status'] = 'FAIL'
            row['errors'].append(str(exc))
            save()
            print('FAIL: '+str(exc)+'; inspect and revalidate if the trajectory changed on this hardware.', flush=True)
            return 1
        save()
        print(f'  PASS: iterations={row["iterations_per_step"]}; solver={row["solver_wrapper_max_seconds"]:.3f}s; application={row["wall_seconds"]:.3f}s', flush=True)
    report['complete'] = True
    save()
    lines = ['# Repeated LEGO timings: '+a.case, '',
        'Warm-ups excluded; seconds are median [minimum, maximum]. Residual checks/profiles retained.', '',
        '| Case | n | Iterations | Solver wrapper s | Application s | Inter-node bytes |',
        '|---|---:|---|---|---|---:|']
    for name, data in report['aggregate'].items():
        def interval(metric):
            v = data[metric]
            return f'{v["median"]:.3f} [{v["min"]:.3f}, {v["max"]:.3f}]'
        lines.append(f'| {name} | {data["samples"]} | {data["iterations"]} | '+
            f'{interval("solver_wrapper_max_seconds")} | {interval("wall_seconds")} | {data["route_bytes"]["inter-node"]:g} |')
    (root/'timings.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)
    print('Results: '+str(root), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
