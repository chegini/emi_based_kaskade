#!/usr/bin/env python3
"""Direct reference and uncompressed CG baseline, then uniform quantization and optional full coding."""
import argparse
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
import sys

from collect_emi_scaling_smoke import collect, NUMBER


def comparison(repo, first, reference, log):
    command = [sys.executable, str(repo/'compare_emi_vtu.py'), str(first), str(reference)]
    with log.open('w') as stream:
        result = subprocess.run(command, cwd=repo, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError('Solution comparison failed; see '+str(log))
    text = log.read_text()
    if not re.search(r'^same coordinates:\s+yes\s*$', text, re.MULTILINE):
        raise RuntimeError('Comparison did not confirm matching coordinates')
    fields = {}
    for label in ['relative L2 error', 'max absolute error',
                  'mean signed offset (first - reference)',
                  'relative L2 after removing mean offset']:
        match = re.search(r'^'+re.escape(label)+r':\s+(\S+)', text, re.MULTILINE)
        if not match or not math.isfinite(float(match.group(1))):
            raise RuntimeError('Missing/nonfinite comparison field: '+label)
        fields[label] = float(match.group(1))
    return fields


def write_summary(root, records):
    (root/'summary.json').write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')
    lines = ['case                 status       Rpasses  iters     true rel     L2 vs direct   L2 vs BDDC   net bytes saved']
    for row in records:
        def number(value):
            return f'{value:.3e}' if value is not None else '--'
        direct = row.get('comparison_to_direct', {}).get('relative L2 error')
        baseline = row.get('comparison_to_uncompressed', {}).get('relative L2 error')
        lines.append(f"{row['case']:20s} {row['status']:12s} {str(row.get('restriction_passes','--')):>7s} {str(row.get('iterations','--')):>5s} "
                     f"{number(row.get('relative_residual')):>12s} "
                     f"{number(direct):>16s} {number(baseline):>12s} "
                     f"{str(round(row['net_logical_byte_savings_vs_uncompressed_percent'],2))+'%' if row.get('net_logical_byte_savings_vs_uncompressed_percent') is not None else '--':>15s}")
    lines += ['', 'Uniform bit depths apply to BOTH restriction and prolongation.',
              'Rpasses counts restriction exchanges per solver update; additional passes refine feedback and cost extra communication.',
              'quantization-* rows disable graph lifting, zigzag, bit-length coding and Huffman.',
              'full-* rows enable integer graph lifting after quantization, zigzag, bit-length coding and Huffman.',
              'Net bytes saved compares total logical transferred bytes with uncompressed CG, including all iterations and refinement passes.',
              'Negative net savings means more bytes were transferred. These counters include intra-rank exchanges and are not hardware network measurements.',
              'Solution differences are raw exported-u Euclidean norms.',
              'Mean-offset diagnostics are saved separately; no offset is removed from the main errors.',
              'PASS means solver residual/completion checks and comparisons completed.',
              'The displayed solution differences require scientific review; no hidden error threshold is applied.']
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--bits', type=int, nargs='+', default=[8,12,16,24,32])
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--refine', type=int, default=0,
                        help='uniform mesh refinement levels, applied to every case (default: 0)')
    parser.add_argument('--tolerance', type=float, default=1e-7)
    parser.add_argument('--iterations', type=int, default=300)
    parser.add_argument('--restriction-passes', type=int, default=1)
    parser.add_argument('--full-compression', action='store_true',
                        help='also run full post-quantization coding at each selected bit depth')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if (not args.bits or len(set(args.bits)) != len(args.bits)
            or any(bits < 1 or bits > 32 for bits in args.bits)):
        parser.error('--bits must be distinct integers in 1..32')
    if args.ranks < 1 or args.iterations < 1:
        parser.error('--ranks and --iterations must be positive')
    if args.refine < 0:
        parser.error('--refine must be nonnegative')
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error('--tolerance must be finite and positive')
    if args.restriction_passes < 1 or args.restriction_passes > 16:
        parser.error('--restriction-passes must be in 1..16')
    repo = args.repo.resolve()
    model = repo/'emiModel'
    if not model.is_file() or not os.access(model, os.X_OK):
        parser.error('Missing executable: '+str(model))
    base = repo/'input/10Cells3d_10extra'
    for suffix in ['_mesh.vtu','_list_extracellular.txt','_list_intracellular.txt','_early_excited.txt']:
        path = Path(str(base)+suffix)
        if not path.is_file():
            parser.error('Missing input: '+str(path))
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = args.output or Path('output')/f'10cells-quantization-diagnostic-{stamp}-{os.getpid()}'
    if not root.is_absolute():
        root = repo/root
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    launcher = shlex.split(env.get('MPIEXEC','mpirun'))
    if not launcher:
        parser.error('MPIEXEC must identify an MPI launcher')
    common = [str(model), '--mpi','1','--sdc','0','--refine',str(args.refine),
              '--input',str(base)+'_mesh.vtu',
              '--extra_set',str(base)+'_list_extracellular.txt',
              '--intra_set',str(base)+'_list_intracellular.txt',
              '--excited',str(base)+'_early_excited.txt',
              '--sigma_i','3','--sigma_e','20','--C_m','1','--R','0.1',
              '--dt','0.01','--finalTime','0.01','--maximumNumberOfTimeSteps','1',
              '--nThreads','1','--quiet','0','--vtk','1','--diagnostics','0','--profile','1',
              '--bddcConditionEstimate','0','--bddcHistogram','0']
    digest = hashlib.sha256()
    with model.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    metadata = {'purpose':'10-cell direct reference, quantization and optional full-coding diagnostic',
                'bits':args.bits, 'bddc_ranks':args.ranks, 'direct_ranks':1,
                'mesh_refinements':args.refine,
                'steps':1, 'true_relative_tolerance':args.tolerance,
                'restriction_passes':args.restriction_passes,
                'max_iterations':args.iterations, 'executable_sha256':digest.hexdigest(),
                'hostname':os.uname().nodename, 'uniform_precision':True,
                'full_encoding_enabled':args.full_compression,
                'byte_counter_scope':'logical application transfers including intra-rank exchanges'}
    (root/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    records = []
    print('Diagnostic output: '+str(root), flush=True)

    def run(name, bits=None, direct=False, full=False):
        folder = root/name
        folder.mkdir()
        ranks = 1 if direct else args.ranks
        options = ['--bddc','0','--direct','1','--directBackend','mumps-spd',
                   '--bddcCompression','0','--bddcConvergenceReport','0'] if direct else [
                   '--bddc','1','--direct','0','--bddcUseCg','1',
                   '--bddcRelativeTolerance',str(args.tolerance),
                   '--bddcIterations',str(args.iterations), '--bddcVerbose','0',
                   '--bddcRestrictionPasses',str(args.restriction_passes if bits is not None else 1),
                   '--bddcConvergenceReport','1','--bddcCompressionReport','1',
                   '--bddcCompression','1' if bits is not None else '0',
                   '--bddcCompressionBits',str(bits if bits is not None else 32),
                   '--bddcRestrictCompressionBits',str(bits if bits is not None else 32),
                   '--bddcProlongateCompressionBits',str(bits if bits is not None else 32),
                   '--masterWorker','0','--bddcSplitInterfaces','0','--bddcInterfaceTypes','7',
                   '--bddcExcludeCornersFromCompression','1',
                   '--bddcMaterialAssembly','1','--bddcInterfaceAssembly','2']
        command = launcher+['--allow-run-as-root','-np',str(ranks),
                            '-x','OPENBLAS_NUM_THREADS','-x','OMP_NUM_THREADS']
        coding = '1' if full else '0'
        encoding_options = ['--bddcGraphLifting',coding,'--bddcGraphLiftingAfterQuantization',coding,
                            '--bddcZigzag',coding,'--bddcBitlength',coding,'--bddcHuffman',coding]
        command += common+options+encoding_options+['--dir',str(folder)]
        (folder/'command.sh').write_text(shlex.join(command)+'\n')
        print('Running '+name+' ...', flush=True)
        with (folder/'run.log').open('w') as stream:
            result = subprocess.run(command,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT)
        text = (folder/'run.log').read_text(errors='replace')
        row = {'case':name, 'restriction_bits':bits, 'prolongation_bits':bits,
               'ranks':ranks, 'return_code':result.returncode, 'full_encoding':full,
               'status':'FAIL','reason':'','log':str(folder/'run.log')}
        if not direct:
            row['restriction_passes'] = args.restriction_passes if bits is not None else 1
        if direct:
            matches = re.findall(rf'Direct global residual step 1: abs=({NUMBER}) rel=({NUMBER})',text)
            relative = float(matches[0][1]) if len(matches)==1 else float('nan')
            if math.isfinite(relative):
                row['relative_residual'] = relative
            if (result.returncode==0 and re.search(r'^EMI model completed\s*$',text,re.M)
                    and len(matches)==1 and math.isfinite(relative) and 0<=relative<=args.tolerance
                    and (folder/'emiLast.vtu').is_file()):
                row['status']='PASS'
            else:
                row['reason']='Direct reference did not complete with a finite residual below target'
        else:
            metrics = collect(text,1,ranks,result.returncode,args.tolerance)
            (folder/'metrics.json').write_text(json.dumps(metrics,indent=2,allow_nan=False)+'\n')
            row['solver_metrics'] = metrics
            row['logical_transferred_bytes'] = metrics['transmitted_sent_bytes']
            baseline_record = next((r for r in records if r['case']=='uncompressed-cg'),None)
            baseline_bytes = (baseline_record.get('logical_transferred_bytes')
                              if baseline_record is not None else metrics['transmitted_sent_bytes'])
            if baseline_bytes and metrics['transmitted_sent_bytes'] is not None:
                row['net_logical_byte_savings_vs_uncompressed_percent'] = (
                    100*(1-metrics['transmitted_sent_bytes']/baseline_bytes))
            history = folder/'bddc-convergence-step1.csv'
            if history.is_file():
                try:
                    with history.open() as stream:
                        samples = list(csv.DictReader(stream))
                    if samples:
                        final = samples[-1]
                        relative = float(final['relative_residual'])
                        row['accepted_updates'] = int(final['iteration'])
                        row['iterations'] = row['accepted_updates']
                        row['iterations_source'] = 'accepted updates in convergence CSV'
                        row['relative_residual'] = relative if math.isfinite(relative) else None
                except (OSError,ValueError,KeyError,csv.Error) as error:
                    row['history_error'] = str(error)
            if metrics['iterations_per_step']:
                row['iterations'] = metrics['iterations_per_step'][-1]
                row['iterations_source'] = 'solver calls in summary'
            if 'relative_residual' not in row and metrics['relative_residuals']:
                row['relative_residual'] = metrics['relative_residuals'][-1]
            if metrics['status']=='converged':
                row['status']='PASS'
            elif 'iteration limit' in text:
                row['reason']='iteration limit'
            elif 'solver breakdown' in text:
                row['reason']='solver breakdown'
            elif 'solver stagnation' in text:
                row['reason']='solver stagnation'
            else:
                row['reason']='incomplete or invalid solver run'
            if row['status']=='PASS':
                solution = folder/('emiBDDCLastCompression.vtu' if bits is not None else 'emiBDDCLast.vtu')
                try:
                    row['comparison_to_direct'] = comparison(repo,solution,root/'direct/emiLast.vtu',folder/'compare-direct.log')
                    if bits is not None:
                        row['comparison_to_uncompressed'] = comparison(repo,solution,root/'uncompressed-cg/emiBDDCLast.vtu',folder/'compare-uncompressed.log')
                except (RuntimeError,OSError,ValueError) as error:
                    row['status']='FAIL'; row['reason']=str(error)
        records.append(row)
        write_summary(root,records)
        suffix = '' if row['status']=='PASS' else ': '+row['reason']
        print('  '+row['status']+suffix,flush=True)
        return row['status']=='PASS'

    if not run('direct',direct=True):
        print('Stopped before BDDC tests: inspect '+str(root/'direct/run.log'))
        return 1
    if not run('uncompressed-cg'):
        print('Stopped before quantization sweep: inspect the uncompressed baseline.')
        return 1
    for bits in args.bits:
        run(f'quantization-{bits}',bits=bits)
        if args.full_compression:
            run(f'full-{bits}',bits=bits,full=True)
    print('\n'+(root/'summary.txt').read_text())
    print('Summary and full logs: '+str(root))
    return int(any(row['status']!='PASS' for row in records))


if __name__ == '__main__':
    sys.exit(main())
