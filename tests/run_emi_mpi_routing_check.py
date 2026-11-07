#!/usr/bin/env python3
"""Paired one-step EMI routing check; unchanged compression, no histograms."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess

from analyze_emi_potential_diagnostics import read, error_metrics
from collect_emi_scaling_smoke import collect, NUMBER
from run_emi_histogram_start_check import command, mesh_inputs, require, snapshots_equal


def arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--case', choices=['10cells', 'pepe'], default='10cells')
    p.add_argument('--ranks', type=int, default=2)
    p.add_argument('--expected-nodes', type=int, default=1)
    p.add_argument('--slurm-container', type=Path,
                   help='Launch each solve with host srun and this Apptainer image')
    p.add_argument('--slurm-btl', choices=['self,vader,tcp', 'self,tcp'],
                   default='self,vader,tcp',
                   help='Open MPI container transports; self,tcp bypasses vader for diagnosis')
    p.add_argument('--kaskade', type=Path,
                   default=Path(os.environ.get('KASKADE7', '/kaskade7')))
    p.add_argument('--solvers', nargs='+', choices=['cg', 'fcg', 'fgmres'],
                   default=['cg', 'fcg', 'fgmres'])
    p.add_argument('--bits', type=int, choices=[12, 16, 24], default=16)
    p.add_argument('--iterations', type=int, default=1000)
    p.add_argument('--transform-order', choices=['before', 'after', 'both'], default='before')
    p.add_argument('--plan-only', action='store_true')
    p.add_argument('--setup-diagnostic-only', action='store_true',
                   help='Run only uncompressed legacy CG with per-rank material progress; no parity claim')
    p.add_argument('--check-local-bypass', action='store_true',
                   help='Also validate exact local/mixed-owner DOFs with compressed remote transfers')
    p.add_argument('--output', type=Path)
    a = p.parse_args(argv)
    if min(a.ranks, a.expected_nodes, a.iterations) < 1 or a.expected_nodes > a.ranks:
        p.error('Require ranks >= expected nodes >= 1 and a positive iteration limit')
    if len(set(a.solvers)) != len(a.solvers):
        p.error('Duplicate solvers')
    if a.slurm_container and a.ranks % a.expected_nodes:
        p.error('Slurm routing check requires equal ranks per node')
    if not a.slurm_container and a.slurm_btl != 'self,vader,tcp':
        p.error('--slurm-btl requires --slurm-container')
    if a.setup_diagnostic_only and a.check_local_bypass:
        p.error('Setup-only and local-bypass validation are separate checks')
    a.repo = a.repo.resolve()
    a.kaskade = a.kaskade.resolve()
    if a.slurm_container:
        a.slurm_container = a.slurm_container.resolve()
    a.steps, a.uncompressed_steps, a.histogram_steps, a.split = 1, 0, [], 'on'
    return a


def planned_modes(a):
    if a.setup_diagnostic_only:
        return [('cg-uncompressed-legacy', 'cg', False, False)]
    modes = [('direct', 'fgmres', False, False)] + [
        (f'{solver}-{kind}-{route}', solver, kind == 'compressed-split', route == 'destination')
        for solver in a.solvers for kind in ('uncompressed', 'compressed-split')
        for route in ('legacy', 'destination')]
    if a.check_local_bypass:
        modes += [(f'{s}-compressed-split-local-bypass', s, True, True) for s in a.solvers]
    return modes


def make_command(a, folder, launcher, solver, compressed, routed, direct=False, bypass=False):
    # Reuse the explicit physical/solver controls and full-precision snapshot export.
    a.solver = solver
    cmd = command(a, 'direct' if direct else 'scheduled' if compressed else 'uncompressed',
                  folder, launcher)
    if compressed and solver in ('cg', 'fcg'):
        cmd[cmd.index('--bddcCgConstantSplitExperiment')+1] = '1'
    cmd += ['--bddcMpiRouting', str(int(routed and not direct))]
    cmd += ['--bddcLocalCompressionBypass', str(int(bypass and not direct))]
    if a.setup_diagnostic_only and not a.slurm_container:
        model_index = cmd.index(str(a.repo/'emiModel'))
        cmd[model_index:model_index] = ['env', 'EMI_MATERIAL_TRACE=1']
    if a.slurm_container:
        # Python orchestrates on the host; every application process is launched
        # by Slurm inside the container. Do not run mpirun inside an srun task.
        model = str(a.repo/'emiModel')
        payload = cmd[cmd.index(model):]
        nodes, ranks = (1, 1) if direct else (a.expected_nodes, a.ranks)
        cmd = slurm_container_command(a, payload, nodes, ranks, mpi='pmix_v5')
    return cmd


CONTAINER_ENV = '''set -eu
export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
export LANG=C LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
exec "$@"
'''


def slurm_container_command(a, payload, nodes, ranks, mpi='none'):
    # Preserve PMIx/Slurm variables: --cleanenv would remove launcher context.
    # Bind at the host path too so recorded absolute inputs/outputs remain valid.
    container_env = CONTAINER_ENV.replace('OMPI_MCA_btl=self,vader,tcp',
                                           'OMPI_MCA_btl='+a.slurm_btl)
    if a.setup_diagnostic_only:
        container_env = container_env.replace('exec "$@"', 'export EMI_MATERIAL_TRACE=1\nexec "$@"')
    return ['srun', '--mpi='+mpi, '--nodes='+str(nodes), '--ntasks='+str(ranks),
            '--ntasks-per-node='+str(ranks//nodes), '--cpus-per-task=1',
            '--distribution=block:block', '--kill-on-bad-exit=1', '--export=ALL',
            '--time=00:10:00', 'apptainer', 'exec',
            '--bind', f'{a.repo.parent}:{a.repo.parent}', '--bind', f'{a.repo.parent}:/BDDC',
            '--bind', f'{a.kaskade}:/kaskade7', '--pwd', str(a.repo), str(a.slurm_container),
            '/bin/bash', '--noprofile', '--norc', '-c', container_env, 'emi-container']+payload


def routing_counts(log, expected_nodes):
    topology = re.findall(r'^BDDC MPI routing: destination; shared-memory domains=(\d+);', log, re.M)
    require(topology == [str(expected_nodes)], 'Missing/unexpected MPI node topology')
    rows = re.findall(r'^BDDC routed payload: phase=(restrict|prolongate) '
                      r'route=(same-rank|same-node|inter-node) bytes=(\d+) packets=(\d+)$', log, re.M)
    expected = {(p, r) for p in ('restrict', 'prolongate')
                for r in ('same-rank', 'same-node', 'inter-node')}
    require(len(rows) == 6 and {(p, r) for p, r, _, _ in rows} == expected,
            'Missing/duplicate routing counters')
    result = [dict(phase=p, route=r, bytes=int(b), packets=int(n)) for p, r, b, n in rows]
    if expected_nodes == 1:
        require(all(r['bytes'] == 0 and r['packets'] == 0 for r in result if r['route'] == 'inter-node'),
                'Single-node run reports inter-node interface packets')
    else:
        require(all(r['bytes'] > 0 and r['packets'] > 0 for r in result if r['route'] == 'inter-node'),
                'Multi-node check did not exercise both R/P inter-node interface transfers')
    return result


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    a = arguments()
    modes = planned_modes(a)
    print('Checks: '+', '.join(m[0] for m in modes), flush=True)
    print('One physical step; histograms OFF; '+
          ('includes local-bypass numerical comparison.' if a.check_local_bypass else
           'routing only; compression unchanged.'), flush=True)
    if a.slurm_container:
        print('Container Open MPI transports: '+a.slurm_btl, flush=True)
    if a.setup_diagnostic_only:
        print('Setup diagnostic: per-rank material progress ON; no direct or routing parity comparison.', flush=True)
    if a.plan_only:
        return 0
    if a.slurm_container:
        require(bool(os.environ.get('SLURM_JOB_ID')), 'Run within an active Slurm allocation')
        require(a.slurm_container.is_file(), 'Missing Apptainer image: '+str(a.slurm_container))
        require(int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) >= a.expected_nodes,
                'Allocation has fewer nodes than requested')
    for filename in mesh_inputs(a.repo, a.case).values():
        require(Path(filename).is_file(), 'Missing mesh input: '+filename)
    model = a.repo/'emiModel'
    model_hash = digest(model)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = (a.output or a.repo/'output'/f'emi-mpi-routing-{a.case}-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    for key in ('KASKADE_RESTRICTION_CAPTURE', 'KASKADE_FROZEN_TRANSFER_DIR', 'KASKADE_TRANSFER_AUDIT_DIR'):
        env.pop(key, None)
    launcher = ['mpirun'] if a.slurm_container else shlex.split(env.get('MPIEXEC', 'mpirun'))
    require(bool(launcher), 'Empty MPIEXEC')
    sources = [a.repo/n for n in ('EMI_model.cpp', 'EMI_bddc.hh', 'EMI_material_assembly.hh',
                'tests/run_emi_mpi_routing_check.py', 'tests/run_emi_histogram_start_check.py')]
    kaskade = a.kaskade
    sources += [kaskade/'mg'/n for n in ('bddc.hh', 'bddc.hpp', 'bddcMpiRouting.hh',
                                       'bddcMpi.hh', 'bddcSpaceTransfer.hh', 'bddcSpaceTransfer.hpp')]
    if a.slurm_container:
        sources.append(a.repo/'tests/build_emi_mpi_routing_cluster.sh')
    manifest = dict(vars(a), executable_sha256=model_hash,
                    slurm={key: env.get(key) for key in
                           ('SLURM_JOB_ID', 'SLURM_JOB_NODELIST', 'SLURM_JOB_NUM_NODES')},
                    source_sha256={str(p): digest(p) for p in sources},
                    note='Validation with snapshot I/O; no performance or network-bandwidth claim.')
    (root/'environment.json').write_text(json.dumps(manifest, indent=2, default=str)+'\n')
    report = dict(runs={}, checks={}, errors=[])

    def save():
        (root/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    print('Output: '+str(root), flush=True)
    if a.slurm_container:
        # Confirm image/source/executable visibility on ALL requested nodes first.
        probe = slurm_container_command(a, ['/bin/bash', '-c',
            'test -r /kaskade7/mg/bddcMpiRouting.hh && test -x "$1" && '
            'printf "Container ready on %s\\n" "$SLURMD_NODENAME"',
            'check', str(model)], a.expected_nodes, a.expected_nodes)
        (root/'preflight-command.sh').write_text(shlex.join(probe)+'\n')
        with (root/'preflight.log').open('w') as stream:
            result = subprocess.run(probe, cwd=a.repo, env=env,
                                    stdout=stream, stderr=subprocess.STDOUT)
        require(result.returncode == 0, 'Container preflight failed; see '+str(root/'preflight.log'))
    for name, solver, compressed, routed in modes:
        folder = root/name
        folder.mkdir()
        require(digest(model) == model_hash, 'Executable changed during validation')
        bypass = name.endswith('-local-bypass')
        cmd = make_command(a, folder, launcher, solver, compressed, routed, name == 'direct', bypass)
        (folder/'command.sh').write_text(shlex.join(['env', 'OPENBLAS_NUM_THREADS=1',
                                                   'OMP_NUM_THREADS=1']+cmd)+'\n')
        print('Running '+name+' ...', flush=True)
        with (folder/'run.log').open('w') as stream:
            rc = subprocess.run(cmd, cwd=a.repo, env=env, stdout=stream, stderr=subprocess.STDOUT).returncode
        log = (folder/'run.log').read_text(errors='replace')
        row = dict(exit_code=rc, status='FAIL')
        report['runs'][name] = row
        try:
            if name == 'direct':
                residuals = re.findall(rf'^Direct global residual step 1: abs=({NUMBER}) rel=({NUMBER})', log, re.M)
                require(rc == 0 and 'EMI model completed' in log and len(residuals) == 1 and
                        math.isfinite(float(residuals[0][1])) and 0 <= float(residuals[0][1]) <= 1e-7,
                        'Direct reference failed')
            else:
                row.update(collect(log, 1, a.ranks, rc, 1e-7))
                require(not row['errors'], '; '.join(row['errors']))
                if not a.setup_diagnostic_only:
                    row['error_vs_direct'] = error_metrics(read(folder/'potential-state-step1.csv')['u'],
                        read(root/'direct/potential-state-step1.csv')['u'])
                if compressed:
                    row['error_vs_uncompressed'] = error_metrics(read(folder/'potential-state-step1.csv')['u'],
                        read(root/f'{solver}-uncompressed-legacy/potential-state-step1.csv')['u'])
            for step in (0, 1):
                read(folder/f'potential-state-step{step}.csv')
            if routed:
                row['routes'] = routing_counts(log, a.expected_nodes)
                require(sum(r['bytes'] for r in row['routes']) == row['transmitted_sent_bytes'],
                        'Routing counters do not sum to total logical payload bytes')
                if bypass:
                    require('BDDC local compression bypass: on;' in log,
                            'Missing local-bypass policy confirmation')
                    if a.expected_nodes == 1:
                        require('nodes=1; backend=uncompressed; all steps use raw transfers' in log,
                                'Single-node bypass did not select the uncompressed backend')
                        row['effective_transfer_backend'] = 'uncompressed'
                    else:
                        row['effective_transfer_backend'] = 'compressed-with-exact-local-dofs'
                    row['error_vs_compressed_all_interfaces'] = error_metrics(
                        read(folder/'potential-state-step1.csv')['u'],
                        read(root/f'{solver}-compressed-split-destination/potential-state-step1.csv')['u'])
                    if a.expected_nodes > 1:
                        # Quantization changes with placement. Do NOT require parity with
                        # the all-compressed case or hide its solution error.
                        report['checks'][name] = 'CONVERGED_WITH_NODE_TOPOLOGY_AND_BYTE_ACCOUNTING'
                        row['status'] = 'PASS'
                        print('  PASS (mixed-node solution differences recorded)', flush=True)
                        save()
                        continue
                legacy_name = (f'{solver}-uncompressed-destination' if bypass else
                               name.removesuffix('destination')+'legacy')
                baseline = root/legacy_name
                snapshots_equal(baseline, folder, (0, 1))
                trace = 'bddc-fgmres-step1.csv' if solver == 'fgmres' else 'bddc-convergence-step1.csv'
                require((baseline/trace).read_bytes() == (folder/trace).read_bytes(), 'Solver trace changed')
                for key in ('iterations_per_step', 'relative_residuals', 'raw_sent_bytes', 'transmitted_sent_bytes'):
                    require(row[key] == report['runs'][legacy_name][key], 'Routing changed '+key)
                report['checks'][name] = 'EXACT_SNAPSHOTS_TRACES_ITERATIONS_AND_BYTES'
            row['status'] = 'PASS'
            print('  PASS', flush=True)
        except (ValueError, OSError, KeyError) as exc:
            row['status'] = 'FAIL'
            report['errors'].append(name+': '+str(exc))
            print('  FAIL: '+str(exc), flush=True)
            save()
            return 1
        save()
    if a.setup_diagnostic_only:
        print('Single uncompressed solve passed; no routing parity comparison. Results: '+str(root), flush=True)
        return 0
    print('All routing parity checks passed. Results: '+str(root), flush=True)
    if a.check_local_bypass:
        print('Local bypass checked: mixed-owner DOFs remain exact across nodes; errors are reported without a hidden acceptance threshold.', flush=True)
    else:
        print('Compression is still applied to local interfaces. Local bypass is opt-in.', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
