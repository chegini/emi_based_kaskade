#!/usr/bin/env python3
"""Validate raw starting steps, continued compressed BDDC and read-only histograms."""
import argparse
from collections import Counter, defaultdict
import csv
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


def arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--case', choices=['10cells', 'pepe'], default='10cells')
    p.add_argument('--ranks', type=int, default=2)
    p.add_argument('--steps', type=int, default=3)
    p.add_argument('--uncompressed-steps', type=int, default=2)
    p.add_argument('--histogram-steps', type=int, nargs='+',
                   help='One-based physical steps; default: final step')
    p.add_argument('--bits', type=int, choices=[4, 6, 8, 12, 16, 20, 24], default=16)
    p.add_argument('--solver', choices=['cg', 'fcg', 'fgmres'], default='fgmres')
    p.add_argument('--split', choices=['off', 'on'], default='on')
    p.add_argument('--transform-order', choices=['before', 'after', 'both'], default='before')
    p.add_argument('--iterations', type=int, default=1000)
    p.add_argument('--plan-only', action='store_true')
    p.add_argument('--output', type=Path)
    a = p.parse_args(argv)
    if a.ranks < 1 or a.steps < 1 or a.iterations < 1:
        p.error('Ranks, steps and iteration limit must be positive')
    if not 0 <= a.uncompressed_steps < a.steps:
        p.error('This switch check requires 0 <= uncompressed steps < total steps')
    if a.split == 'on' and (a.solver != 'fgmres' or a.steps > 10):
        p.error('Multi-step split validation currently supports FGMRES with at most 10 steps; use --split off otherwise')
    a.histogram_steps = sorted(a.histogram_steps or [a.steps])
    if (len(set(a.histogram_steps)) != len(a.histogram_steps) or
            any(s <= a.uncompressed_steps or s > a.steps for s in a.histogram_steps)):
        p.error('Histogram steps must be distinct compressed physical steps within this run')
    return a


def mesh_inputs(repo, case):
    base = (repo/'input/10Cells3d_10extra' if case == '10cells'
            else repo.parent/'input_emi_mesh/pepe_combi_domi')
    suffixes = (['_mesh.vtu', '_list_extracellular.txt', '_list_intracellular.txt', '_early_excited.txt']
                if case == '10cells' else
                ['.vtu', '_extracellular.txt', '_intracellular.txt', '_excited.txt'])
    return dict(zip(['input', 'extra_set', 'intra_set', 'excited'],
                    [str(base)+s for s in suffixes]))


def command(a, name, out, launcher):
    direct = name == 'direct'
    compressed = name not in ('direct', 'uncompressed')
    capture = name == 'scheduled-histograms'
    warmup = (a.steps if name == 'scheduled-all-raw' else
              a.uncompressed_steps if name.startswith('scheduled') else 0)
    solver = a.solver
    # Same mesh, physical problem, initial state, dt and solver-family tolerance.
    # Defaults are made explicit for every setting used by the observer gate.
    values = dict(mpi=1, sdc=0, refine=0, **mesh_inputs(a.repo, a.case),
        sigma_i=3, sigma_e=20, C_m=1, R=0.1, dt=0.01,
        finalTime=format(a.steps*0.01, '.12g'), maximumNumberOfTimeSteps=a.steps,
        nThreads=1, quiet=0, vtk=0, diagnostics=0, profile=1, potentialDiagnostics=1,
        bddc=int(not direct), direct=int(direct),
        bddcUseCg=int(not direct and solver != 'fgmres'),
        bddcFlexibleCg=int(not direct and solver == 'fcg'),
        bddcCgOriginalResidual=int(not direct and solver != 'fgmres'),
        bddcFgmres=int(not direct and solver == 'fgmres'), bddcFgmresRestart=30,
        bddcProlongateConstantSplit=int(compressed and a.split == 'on'),
        bddcCgConstantSplitExperiment=0, bddcConstantCorrection=0,
        bddcRelativeTolerance=0 if direct else 1e-7, bddcIterations=a.iterations,
        bddcVerbose=0, bddcConvergenceReport=int(not direct and solver != 'fgmres'),
        bddcCompressionReport=1, bddcCompression=int(compressed),
        bddcCompressionBits=a.bits, bddcRestrictCompressionBits=a.bits,
        bddcProlongateCompressionBits=a.bits, bddcRestrictionPasses=1,
        bddcRestrictZeroPreserving=0, bddcProlongateZeroPreserving=0,
        bddcRestrictGroupSize=0, bddcProlongateGroupSize=0,
        bddcCompressionMinInterfaceSize=0, bddcConditionEstimate=0,
        bddcCgDiagnostics=0, bddcProlongationDiagnostics=0,
        bddcOperatorAudit=0, bddcPreconditionerProbe=0,
        masterWorker=0, bddcSplitInterfaces=0, bddcInterfaceTypes=7,
        bddcExcludeCornersFromCompression=1,
        bddcMaterialAssembly=int(not direct), bddcInterfaceAssembly=0 if direct else 2,
        bddcGraphLifting=int(compressed and a.transform_order in ('before', 'both')),
        bddcGraphLiftingAfterQuantization=int(compressed and a.transform_order in ('after', 'both')),
        bddcGraphLiftingBothOrders=int(compressed and a.transform_order == 'both'),
        bddcZigzag=int(compressed), bddcBitlength=int(compressed), bddcHuffman=int(compressed),
        bddcHistogram=0, bddcUncompressedSteps=warmup, dir=str(out))
    if direct:
        values['directBackend'] = 'mumps-spd'
    if capture:
        values['bddcHistogramSteps'] = ','.join(map(str, a.histogram_steps))
    result = launcher+['--allow-run-as-root', '-np', str(1 if direct else a.ranks),
                       '-x', 'OPENBLAS_NUM_THREADS', '-x', 'OMP_NUM_THREADS', str(a.repo/'emiModel')]
    for key, value in values.items():
        result += ['--'+key, str(value)]
    return result


def require(condition, message):
    if not condition:
        raise ValueError(message)


def collect_bddc_run(a, name, log, exit_code):
    # A type switch constructs a second solver on the SAME interface topology.
    # Keep the ordinary collector strict; opt into two setups only for this run.
    raw = a.steps if name == 'scheduled-all-raw' else (
        a.uncompressed_steps if name.startswith('scheduled') else 0)
    phases = []
    if raw:
        phases.append((1, raw, 'uncompressed'))
        if raw < a.steps:
            phases.append((raw+1, a.steps, 'compressed'))
    row = collect(log, a.steps, a.ranks, exit_code, 1e-7,
                  expected_solver_setups=max(1, len(phases)))
    observed = re.findall(r'^BDDC transfer phase: steps (\d+)\.\.(\d+); '
                          r'(uncompressed|compressed); state retained; physical time=', log, re.M)
    if [(int(first), int(last), mode) for first, last, mode in observed] != phases:
        row['errors'].append('Transfer phases do not match the requested physical-step schedule')
    row['status'] = 'failed' if row['errors'] else 'converged'
    return row


def snapshots_equal(first, second, steps):
    for step in steps:
        name = f'potential-state-step{step}.csv'
        require((first/name).read_bytes() == (second/name).read_bytes(),
                f'Full-precision state/RHS/correction/residual differs at step {step}: {first.name}/{second.name}')


def summarize_histograms(folder, selected):
    hist = folder/'symbol-histograms'
    require(hist.is_dir(), 'No selected-step histogram directory')
    require({p.name for p in hist.iterdir()} == {f'step{s}' for s in selected},
            'Unexpected or missing physical-step histogram output')
    counts = defaultdict(Counter)
    packet_counts = Counter()
    fallback_counts = Counter()
    for step in selected:
        base = hist/f'step{step}'
        metadata = dict(line.split('=', 1) for line in (base/'metadata.txt').read_text().splitlines())
        require(int(metadata['physical_step']) == step, 'Histogram metadata step mismatch')
        require(math.isclose(float(metadata['start_time']), (step-1)*0.01, abs_tol=1e-14),
                'Histogram time reset at switch')
        for path in sorted(base.glob('subdomain*/*.csv')):
            phase = path.name.split('-exchange')[0]
            require(phase in ('restrict', 'prolongate'), 'Unknown histogram phase')
            with path.open() as stream:
                rows = list(csv.DictReader(stream))
            if not rows:  # A packet containing only exact values has no quantized symbols.
                continue
            encoded = {r['huffman_encoded'] for r in rows}
            require(encoded in ({'0'}, {'1'}), 'Mixed encoding flags in one packet')
            packet_counts[step, phase] += 1
            fallback_counts[step, phase] += int(encoded == {'0'})
            totals = Counter()
            for row in rows:
                stage, role = row['stage'], row['role']
                require(stage in ('quantized_input', 'raw', 'unsigned', 'bitlength'), 'Unknown histogram stage')
                count, symbol = int(row['count']), int(row['symbol'])
                require(count > 0, 'Invalid symbol count')
                totals[stage] += count
                counts[step, phase, stage, role][symbol] += count
            require(len(totals) == 4 and len(set(totals.values())) == 1,
                    'Histogram stage populations disagree')
        for phase in ('restrict', 'prolongate'):
            require(packet_counts[step, phase] > 0, f'No {phase} symbols captured at step {step}')
    rows, metrics = [], []
    for (step, phase, stage, role), distribution in sorted(counts.items()):
        total = sum(distribution.values())
        entropy = -math.fsum((n/total)*math.log2(n/total) for n in distribution.values())
        metrics.append(dict(step=step, phase=phase, stage=stage, role=role,
                            symbols=total, distinct_symbols=len(distribution),
                            empirical_entropy_bits_per_symbol=entropy,
                            packet_count=packet_counts[step, phase],
                            fallback_packet_count=fallback_counts[step, phase]))
        rows.extend(dict(step=step, phase=phase, stage=stage, role=role, symbol=s, count=n)
                    for s, n in sorted(distribution.items()))
    with (folder/'histogram-counts.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['step', 'phase', 'stage', 'role', 'symbol', 'count'])
        writer.writeheader()
        writer.writerows(rows)
    return metrics


def main(argv=None):
    a = arguments(argv)
    a.repo = a.repo.resolve()
    names = ['direct', 'uncompressed', 'scheduled-all-raw', 'scheduled', 'scheduled-histograms']
    print('Checks: '+', '.join(names), flush=True)
    print(f'Raw steps: 1..{a.uncompressed_steps}; compressed thereafter; histogram steps: {a.histogram_steps}', flush=True)
    if a.plan_only:
        return 0
    model = a.repo/'emiModel'
    require(model.is_file() and os.access(model, os.X_OK), 'Build emiModel first')
    for path in mesh_inputs(a.repo, a.case).values():
        require(Path(path).is_file(), 'Missing mesh input: '+path)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = a.output or a.repo/'output'/f'emi-histogram-start-{a.case}-{stamp}-{os.getpid()}'
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    for key in ('KASKADE_RESTRICTION_CAPTURE', 'KASKADE_FROZEN_TRANSFER_DIR', 'KASKADE_TRANSFER_AUDIT_DIR'):
        env.pop(key, None)
    launcher = shlex.split(env.get('MPIEXEC', 'mpirun'))
    require(bool(launcher), 'Empty MPIEXEC')
    digest = hashlib.sha256(model.read_bytes()).hexdigest()
    source_paths = [a.repo/name for name in ('EMI_model.cpp', 'EMI_bddc.hh', 'EMI_bddc_step_schedule.hh',
                                            'tests/run_emi_histogram_start_check.py')]
    source_paths += [Path(env.get('KASKADE7', '/kaskade7'))/'mg'/name
                     for name in ('bddc.hh', 'bddcSpaceTransfer.hh', 'bddcSpaceTransfer.hpp', 'bddcSymbolHistogram.hh')]
    manifest = dict(vars(a), hostname=os.uname().nodename, executable_sha256=digest, dt=0.01,
                    planned_runs=names, source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                                                      for p in source_paths if p.is_file()},
                    note='Validation timings include snapshots/histograms and one solver rebuild at the switch; not a performance benchmark.')
    (root/'environment.json').write_text(json.dumps(manifest, indent=2, default=str)+'\n')
    report = dict(runs={}, checks={}, errors=[])
    def save():
        (root/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print('Output: '+str(root), flush=True)
    for name in names:
        folder = root/name
        folder.mkdir()
        cmd = command(a, name, folder, launcher)
        (folder/'command.sh').write_text(shlex.join(['env', 'OPENBLAS_NUM_THREADS=1', 'OMP_NUM_THREADS=1']+cmd)+'\n')
        require(hashlib.sha256(model.read_bytes()).hexdigest() == digest, 'Executable changed during validation')
        print('Running '+name+' ...', flush=True)
        with (folder/'run.log').open('w') as stream:
            rc = subprocess.run(cmd, cwd=a.repo, env=env, stdout=stream, stderr=subprocess.STDOUT).returncode
        log = (folder/'run.log').read_text(errors='replace')
        if name == 'direct':
            residuals = re.findall(rf'^Direct global residual step (\d+): abs=({NUMBER}) rel=({NUMBER})', log, re.M)
            errors = []
            if (rc or 'EMI model completed' not in log or
                    [int(r[0]) for r in residuals] != list(range(1, a.steps+1)) or
                    any(not math.isfinite(float(r[1])) or float(r[1]) < 0 or
                        not math.isfinite(float(r[2])) or not 0 <= float(r[2]) <= 1e-7 for r in residuals)):
                errors.append('Direct reference failed residual/completion check')
            row = dict(exit_code=rc, errors=errors)
        else:
            row = collect_bddc_run(a, name, log, rc)
        report['runs'][name] = row
        try:
            for step in range(a.steps+1):
                read(folder/f'potential-state-step{step}.csv')
            if name != 'direct':
                row['errors_vs_direct'] = [dict(step=step, **error_metrics(
                    read(folder/f'potential-state-step{step}.csv')['u'],
                    read(root/'direct'/f'potential-state-step{step}.csv')['u']))
                    for step in range(1, a.steps+1)]
        except (OSError, ValueError, KeyError) as exc:
            row['errors'].append(str(exc))
        save()
        print('  '+('FAIL: '+'; '.join(row['errors']) if row['errors'] else 'PASS'), flush=True)
        if row['errors']:
            report['errors'].append('Stopped after failed '+name+'; see its run.log')
            save()
            return 1
    try:
        snapshots_equal(root/'scheduled-all-raw', root/'uncompressed', range(a.steps+1))
        report['checks']['all_raw_schedule'] = 'EXACT_SNAPSHOTS'
        snapshots_equal(root/'scheduled', root/'uncompressed', range(a.uncompressed_steps+1))
        report['checks']['uncompressed_prefix'] = 'EXACT_SNAPSHOTS'
        next_snapshot = f'potential-state-step{a.uncompressed_steps+1}.csv'
        require(read(root/'scheduled'/next_snapshot)['rhs'] == read(root/'uncompressed'/next_snapshot)['rhs'],
                'First compressed step does not have the same physical RHS as its raw reference')
        report['checks']['switch_rhs'] = 'EXACT_FIRST_COMPRESSED_STEP_RHS'
        snapshots_equal(root/'scheduled', root/'scheduled-histograms', range(a.steps+1))
        prefix = 'bddc-fgmres-step' if a.solver == 'fgmres' else 'bddc-convergence-step'
        expected = {f'{prefix}{s}.csv' for s in range(1, a.steps+1)}
        for name in ('scheduled', 'scheduled-histograms'):
            require({p.name for p in (root/name).glob(prefix+'*.csv')} == expected,
                    'Missing or unexpected physical-step solver traces')
        for name in expected:
            require((root/'scheduled'/name).read_bytes() == (root/'scheduled-histograms'/name).read_bytes(),
                    'Histogram observer changed convergence trace')
        for left, right in [('scheduled', 'scheduled-histograms'), ('uncompressed', 'scheduled-all-raw')]:
            for key in ('iterations_per_step', 'relative_residuals', 'raw_sent_bytes', 'transmitted_sent_bytes'):
                require(report['runs'][left][key] == report['runs'][right][key], 'Observer/schedule changed '+key)
        require(not (root/'scheduled'/'symbol-histograms').exists(), 'Observer-off run wrote histograms')
        report['checks']['observer'] = 'EXACT_SNAPSHOTS_TRACES_ITERATIONS_AND_BYTES'
        report['histograms'] = summarize_histograms(root/'scheduled-histograms', a.histogram_steps)
        report['errors_vs_uncompressed'] = [dict(step=step, **error_metrics(
            read(root/'scheduled'/f'potential-state-step{step}.csv')['u'],
            read(root/'uncompressed'/f'potential-state-step{step}.csv')['u']))
            for step in range(1, a.steps+1)]
    except (OSError, ValueError, KeyError) as exc:
        report['errors'].append(str(exc))
    save()
    lines = ['case | iterations per step | logical bytes | solver wrapper seconds | final relative L2/direct']
    for name in names:
        row = report['runs'][name]
        lines.append(' | '.join(map(str, [name, row.get('iterations_per_step', '--'),
            row.get('transmitted_sent_bytes', '--'), row.get('solver_wrapper_max_seconds', '--'),
            row.get('errors_vs_direct', [{}])[-1].get('relative_l2', '--')])))
    lines.append('Checks: '+json.dumps(report['checks']))
    lines.extend('FAIL: '+error for error in report['errors'])
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)
    print('Results: '+str(root))
    print('Solution errors are reported, without a hidden acceptance threshold. Logical byte counters include intra-rank transfers.')
    return int(bool(report['errors']))


if __name__ == '__main__':
    raise SystemExit(main())
