#!/usr/bin/env python3
"""Matched EMI direct/CG/FCG/FGMRES comparison with transform and constant-split choices."""
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
    lines += ['', 'case                 interface bytes   solver wrapper s']
    for row in records:
        seconds=row.get('solver_metrics',{}).get('solver_wrapper_max_seconds')
        lines.append(f"{row['case']:20s} {str(row.get('logical_transferred_bytes') if row.get('logical_transferred_bytes') is not None else '--'):>15s} "
                     f"{format(seconds,'.3f') if seconds is not None else '--':>18s}")
    lines += ['', 'case                           centered L2/direct  max error/direct  centered L2/BDDC  max error/BDDC']
    for row in records:
        d=row.get('comparison_to_direct',{})
        b=row.get('comparison_to_uncompressed',{})
        def error_value(values,key):
            x=values.get(key)
            return format(x,'.3e') if x is not None else '--'
        lines.append(f"{row['case']:30s} "
            f"{error_value(d,'relative L2 after removing mean offset'):>18s} "
            f"{error_value(d,'max absolute error'):>17s} "
            f"{error_value(b,'relative L2 after removing mean offset'):>17s} "
            f"{error_value(b,'max absolute error'):>15s}")
    for row in records:
        if row.get('restriction_bits') is not None:
            lines.append(f"Precision {row['case']}: restriction={row['restriction_bits']}, "
                         f"prolongation={row['prolongation_bits']}, "
                         f"restriction_group_size={row.get('restriction_group_size',0)}")
    lines += ['', 'Directional bit depths and producer restriction group sizes are recorded above and in summary.json.',
              'Zero-preserving restriction mode is recorded per case in summary.json; prolongation retains its original quantizer.',
              'Rpasses counts restriction exchanges per solver update; additional passes refine feedback and cost extra communication.',
              'quantization-* rows disable graph lifting, zigzag, bit-length coding and Huffman.',
              'full-* and after-* rows enable integer graph lifting after quantization and full coding.',
              'before-* rows apply floating graph lifting before affine quantization and full coding.',
              'Solvers: cg=CG, fcg=flexible CG, fgmres=flexible GMRES (restart 30); -split adds constant splitting.',
              'Uncompressed controls have no split; each compressed case uses its matched solver-family baseline.',
              'Net bytes saved compares total logical transferred bytes with the same uncompressed solver, including all iterations and refinement passes.',
              'Negative net savings means more bytes were transferred. These counters include intra-rank exchanges and are not hardware network measurements.',
              'Solution differences are raw exported-u Euclidean norms.',
              'Raw L2 errors retain mean offsets; the separate centered table removes the arithmetic mean error.',
              'PASS means solver residual/completion checks and comparisons completed.',
              'The displayed solution differences require scientific review; no hidden error threshold is applied.']
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--case', choices=['10cells','pepe'], default='10cells')
    parser.add_argument('--bits', type=int, nargs='+', default=[8,12,16,24,32])
    parser.add_argument('--solvers', nargs='+', choices=['cg','fcg','fgmres'], default=['cg'])
    parser.add_argument('--split', choices=['off','on','both'], default='off',
                        help='Constant-split prolongation on compressed cases only')
    parser.add_argument('--skip-quantization', action='store_true',
                        help='Omit the quantization-only control; keep selected full-codec orders')
    parser.add_argument('--plan-only', action='store_true',
                        help='Print validated experiment names without building/running')
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--refine', type=int, default=0,
                        help='uniform mesh refinement levels, applied to every case (default: 0)')
    parser.add_argument('--tolerance', type=float, default=1e-7)
    parser.add_argument('--iterations', type=int, default=300)
    parser.add_argument('--restriction-passes', type=int, default=1)
    parser.add_argument('--prolongation-bits', type=int,
                        help='override prolongation only; rows are labelled mixed precision')
    parser.add_argument('--restriction-group-size', type=int, default=0,
                        help='producer group size; requires zero-preserving restriction')
    parser.add_argument('--full-compression', action='store_true',
                        help='also run full post-quantization coding at each selected bit depth')
    parser.add_argument('--transform-orders', nargs='+', choices=['before','after'],
                        help='before=transform then quantize; after=quantize then integer lift')
    parser.add_argument('--cg-original-residual', type=int, choices=[0,1], default=0,
                        help='CG scalar convention; 1 matches the previous Pepe split campaign')
    parser.add_argument('--zero-preserving-restriction', action='store_true',
                        help='experimental symmetric exact-zero restriction; prolongation unchanged')
    parser.add_argument('--cg-diagnostics', action='store_true',
                        help='record restriction error, CG scalar error and adjacent conjugacy; adds validation work')
    parser.add_argument('--restriction-range-capture', action='store_true',
                        help='read-only pre-quantization buffer capture; adds file I/O, not MPI')
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
    if args.zero_preserving_restriction and min(args.bits)<2:
        parser.error('zero-preserving restriction requires 2..32 bits')
    if args.prolongation_bits is not None and not 1<=args.prolongation_bits<=32:
        parser.error('--prolongation-bits must be in 1..32')
    if args.restriction_group_size<0 or (args.restriction_group_size and not args.zero_preserving_restriction):
        parser.error('restriction groups require a nonnegative size and zero-preserving restriction')
    if args.transform_orders:
        if args.full_compression or len(set(args.transform_orders))!=len(args.transform_orders):
            parser.error('Use distinct --transform-orders without --full-compression')
        if 'before' in args.transform_orders and (
                args.zero_preserving_restriction or args.restriction_group_size):
            parser.error('Transform-first requires ungrouped affine quantization')
    if len(set(args.solvers))!=len(args.solvers):
        parser.error('Duplicate solvers')
    if args.skip_quantization and not (args.full_compression or args.transform_orders):
        parser.error('--skip-quantization requires a full-codec comparison')
    if 'fgmres' in args.solvers and (args.restriction_passes!=1 or args.cg_diagnostics):
        parser.error('FGMRES requires one restriction pass and no CG diagnostics')
    if args.split!='off':
        if (args.zero_preserving_restriction or args.restriction_group_size or
                args.cg_diagnostics or args.prolongation_bits is not None):
            parser.error('Split comparison requires uniform affine precision without grouped/CG diagnostics')
        if any(x in args.solvers for x in ('cg','fcg')):
            if (args.restriction_passes!=1 or args.tolerance!=1e-7 or
                    any(b not in (4,6,8,12,16,20,24,32) for b in args.bits) or
                    ('cg' in args.solvers and not args.cg_original_residual)):
                parser.error('CG/FCG split requires original-residual CG, one R pass, tolerance 1e-7, '
                             'and supported experimental precision')
    # Preserve old single-CG/off directory names for existing diagnostic users.
    legacy_names=args.solvers==['cg'] and args.split=='off'
    pipelines=([] if args.skip_quantization else ['quantization'])
    pipelines+=['full'] if args.full_compression else (args.transform_orders or [])
    split_choices={'off':[False],'on':[True],'both':[False,True]}[args.split]
    planned=[(f'{mode}-{bits}' if legacy_names else
              f'{mode}-{solver}{bits}'+('-split' if split else ''),
              bits,solver,split,mode)
             for bits in args.bits for solver in args.solvers
             for mode in pipelines for split in split_choices]
    print('Planned runs: direct, '+', '.join(
        ['uncompressed-'+solver for solver in args.solvers]+[x[0] for x in planned]),flush=True)
    if args.plan_only:
        return 0
    repo = args.repo.resolve()
    model = repo/'emiModel'
    if not model.is_file() or not os.access(model, os.X_OK):
        parser.error('Missing executable: '+str(model))
    base = (repo/'input/10Cells3d_10extra' if args.case=='10cells'
            else repo.parent/'input_emi_mesh/pepe_combi_domi')
    suffixes = (['_mesh.vtu','_list_extracellular.txt','_list_intracellular.txt','_early_excited.txt']
                if args.case=='10cells' else
                ['.vtu','_extracellular.txt','_intracellular.txt','_excited.txt'])
    inputs = {flag:str(base)+suffix for flag,suffix in
              zip(['--input','--extra_set','--intra_set','--excited'],suffixes)}
    for path in inputs.values():
        if not Path(path).is_file():
            parser.error('Missing input: '+path)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = args.output or Path('output')/f'{args.case}-quantization-diagnostic-{stamp}-{os.getpid()}'
    if not root.is_absolute():
        root = repo/root
    root.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE', None)
    launcher = shlex.split(env.get('MPIEXEC','mpirun'))
    if not launcher:
        parser.error('MPIEXEC must identify an MPI launcher')
    common = [str(model), '--mpi','1','--sdc','0','--refine',str(args.refine),
              '--input',inputs['--input'],'--extra_set',inputs['--extra_set'],
              '--intra_set',inputs['--intra_set'],'--excited',inputs['--excited'],
              '--sigma_i','3','--sigma_e','20','--C_m','1','--R','0.1',
              '--dt','0.01','--finalTime','0.01','--maximumNumberOfTimeSteps','1',
              '--nThreads','1','--quiet','0','--vtk','1','--diagnostics','0','--profile','1',
              '--bddcConditionEstimate','0','--bddcHistogram','0']
    digest = hashlib.sha256()
    with model.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    metadata = {'purpose':'Matched EMI direct reference, quantization and transform-order diagnostic',
                'mesh_case':args.case,'mesh_input_files':inputs,
                'transform_orders':args.transform_orders,'cg_original_residual':args.cg_original_residual,
                'solvers':args.solvers,'constant_split':args.split,'skip_quantization':args.skip_quantization,
                'planned_cases':[x[0] for x in planned],
                'baseline_policy':'Same solver family, mesh, ranks, physics and tolerance; no compression/split',
                'bits':args.bits, 'bddc_ranks':args.ranks, 'direct_ranks':1,
                'mesh_refinements':args.refine,
                'steps':1, 'true_relative_tolerance':args.tolerance,
                'restriction_passes':args.restriction_passes,
                'max_iterations':args.iterations, 'executable_sha256':digest.hexdigest(),
                'hostname':os.uname().nodename, 'uniform_precision':args.prolongation_bits is None,
                'prolongation_bits_override':args.prolongation_bits,
                'restriction_group_size':args.restriction_group_size,
                'full_encoding_enabled':args.full_compression,
                'cg_diagnostics':args.cg_diagnostics,
                'restriction_range_capture':args.restriction_range_capture,
                'zero_preserving_restriction':args.zero_preserving_restriction,
                'byte_counter_scope':'logical application transfers including intra-rank exchanges'}
    if args.transform_orders:
        kaskade=Path(env.get('KASKADE7','/kaskade7'))
        sources=[repo/name for name in
                 ('EMI_model.cpp','EMI_bddc.hh','EMI_bddc_precision.hh',
                  'tests/run_emi_quantization_diagnostic.py')]
        sources += [kaskade/'mg'/name for name in
                    ('bddcSpaceTransfer.hh','bddcSpaceTransfer.hpp','bddcLiftingPlan.hh')]
        metadata['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sources if p.is_file()}
    (root/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    records = []
    print('Diagnostic output: '+str(root), flush=True)

    def run(name, bits=None, direct=False, full=False, order=None, solver='cg', split=False):
        folder = root/name
        folder.mkdir()
        ranks = 1 if direct else args.ranks
        options = ['--bddc','0','--direct','1','--directBackend','mumps-spd',
                   '--bddcCompression','0','--bddcConvergenceReport','0'] if direct else [
                   '--bddc','1','--direct','0','--bddcUseCg',str(int(solver!='fgmres')),
                   '--bddcCgOriginalResidual',str(0 if solver=='fgmres' else
                                                   1 if solver=='fcg' else args.cg_original_residual),
                   '--bddcFlexibleCg',str(int(solver=='fcg')),
                   '--bddcFgmres',str(int(solver=='fgmres')),'--bddcFgmresRestart','30',
                   '--bddcProlongateConstantSplit',str(int(split)),
                   '--bddcRelativeTolerance',str(args.tolerance),
                   '--bddcIterations',str(args.iterations), '--bddcVerbose','0',
                   '--bddcRestrictionPasses',str(args.restriction_passes if bits is not None else 1),
                   '--bddcConvergenceReport',str(int(solver!='fgmres')),'--bddcCompressionReport','1',
                   '--bddcCompression','1' if bits is not None else '0',
                   '--bddcCompressionBits',str(bits if bits is not None else 32),
                   '--bddcRestrictCompressionBits',str(bits if bits is not None else 32),
                   '--bddcProlongateCompressionBits',str(
                       (args.prolongation_bits or bits) if bits is not None else 32),
                   '--masterWorker','0','--bddcSplitInterfaces','0','--bddcInterfaceTypes','7',
                   '--bddcExcludeCornersFromCompression','1',
                   '--bddcMaterialAssembly','1','--bddcInterfaceAssembly','2']
        command = launcher+['--allow-run-as-root','-np',str(ranks),
                            '-x','OPENBLAS_NUM_THREADS','-x','OMP_NUM_THREADS']
        case_env = dict(env)
        capture_prefix = ''
        if args.restriction_range_capture and bits is not None:
            capture = str(folder/'restriction-buffers')
            case_env['KASKADE_RESTRICTION_CAPTURE'] = capture
            command += ['-x','KASKADE_RESTRICTION_CAPTURE']
            capture_prefix = 'env '+shlex.quote('KASKADE_RESTRICTION_CAPTURE='+capture)+' '
        coding = '1' if full or order else '0'
        # Preserve the legacy full-compression command; explicit orders are unambiguous.
        lifting = '1' if full or order=='before' else '0'
        after = '1' if full or order=='after' else '0'
        encoding_options = ['--bddcGraphLifting',lifting,'--bddcGraphLiftingAfterQuantization',after,
                            '--bddcZigzag',coding,'--bddcBitlength',coding,'--bddcHuffman',coding]
        command += common+options+encoding_options+['--dir',str(folder)]
        if split and solver in ('cg','fcg'):
            command += ['--bddcCgConstantSplitExperiment','1']

        if args.zero_preserving_restriction and bits is not None:
            command += ['--bddcRestrictZeroPreserving','1']
        if args.restriction_group_size and bits is not None:
            command += ['--bddcRestrictGroupSize',str(args.restriction_group_size)]
        if args.cg_diagnostics and not direct:
            command += ['--bddcCgDiagnostics','1']
        (folder/'command.sh').write_text(capture_prefix+shlex.join(command)+'\n')
        print('Running '+name+' ...', flush=True)
        with (folder/'run.log').open('w') as stream:
            result = subprocess.run(command,cwd=repo,env=case_env,stdout=stream,stderr=subprocess.STDOUT)
        text = (folder/'run.log').read_text(errors='replace')
        row = {'case':name, 'restriction_bits':bits,
               'prolongation_bits':(args.prolongation_bits or bits) if bits is not None else None,
               'restriction_group_size':args.restriction_group_size if bits is not None else 0,
               'ranks':ranks, 'return_code':result.returncode, 'full_encoding':bool(full or order),
               'transform_order':order,'solver':'direct' if direct else solver,'constant_split':split,
               'uncompressed_reference':None if direct else 'uncompressed-'+solver,
               'status':'FAIL','reason':'','log':str(folder/'run.log')}
        row['zero_preserving_restriction'] = bool(args.zero_preserving_restriction and bits is not None)
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
            baseline_record = next((r for r in records if r['case']=='uncompressed-'+solver),None)
            baseline_bytes = (baseline_record.get('logical_transferred_bytes')
                              if baseline_record is not None else metrics['transmitted_sent_bytes'])
            if baseline_bytes and metrics['transmitted_sent_bytes'] is not None:
                row['net_logical_byte_savings_vs_uncompressed_percent'] = (
                    100*(1-metrics['transmitted_sent_bytes']/baseline_bytes))
            history = folder/('bddc-fgmres-step1.csv' if solver=='fgmres' else 'bddc-convergence-step1.csv')
            if history.is_file():
                try:
                    with history.open() as stream:
                        samples = list(csv.DictReader(stream))
                    if samples:
                        final = samples[-1]
                        relative = float(final.get('true_relative_residual',final.get('relative_residual')))
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
                        row['comparison_to_uncompressed'] = comparison(repo,solution,root/('uncompressed-'+solver)/'emiBDDCLast.vtu',folder/'compare-uncompressed.log')
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
    baseline_ok={solver:run('uncompressed-'+solver,solver=solver) for solver in args.solvers}
    for name,bits,solver,split,mode in planned:
        if not baseline_ok[solver]:
            records.append(dict(case=name,status='SKIP',solver=solver,constant_split=split,
                                reason='Matched uncompressed baseline failed'))
            write_summary(root,records)
            continue
        run(name,bits=bits,solver=solver,split=split,full=mode=='full',
            order=mode if mode in ('before','after') else None)
    print('\n'+(root/'summary.txt').read_text())
    print('Summary and full logs: '+str(root))
    return int(any(row['status']!='PASS' for row in records))


if __name__ == '__main__':
    sys.exit(main())
