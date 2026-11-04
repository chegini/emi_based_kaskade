#!/usr/bin/env python3
"""Experimental paired CG/FCG constant-split tests; one step, selectable uniform affine precision."""
import argparse
from datetime import datetime,timezone
import csv,hashlib,json,math,os
from pathlib import Path
import re,shlex,subprocess,sys


REFERENCE_CAMPAIGNS={
    '10cells':'output/emi-fgmres-20261101T205726Z-16111',
    'weak_8_h16':'output/emi-fgmres-20261101T210711Z-16464',
    'strong_h16':'output/emi-fgmres-20261101T211157Z-16694',
}


def reference_directory(case, campaign=None):
    if case=='pepe':return None
    return Path(campaign or REFERENCE_CAMPAIGNS[case]).resolve()/case


def validate_scaling_source(case, command):
    if case not in ('weak_8_h16','strong_h16'):return
    expected={'--input':'_mesh.vtu','--extra_set':'_list_extracellular.txt',
              '--intra_set':'_list_intracellular.txt','--excited':'_early_excited.txt'}
    for flag,suffix in expected.items():
        if command.count(flag)!=1 or Path(command[command.index(flag)+1]).name!=case+'_kaskade'+suffix:
            raise ValueError('Saved source does not match '+case+': '+flag)
    for flag,value in {'--refine':0,'-np':8,'--sigma_i':3,'--sigma_e':20,
                       '--C_m':1,'--R':0.1,'--dt':0.01}.items():
        if command.count(flag)!=1 or float(command[command.index(flag)+1])!=value:
            raise ValueError('Unexpected scaling configuration: '+flag)


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def make_command(saved,out,mode,pipeline="full"):
    if pipeline not in ("full","quantization-only"):
        raise ValueError("Unknown compression pipeline")
    cmd=list(saved)
    def set_value(flag,new):
        if cmd.count(flag)!=1:raise ValueError('Expected exactly one '+flag)
        cmd[cmd.index(flag)+1]=str(new)
    if cmd.count('--bddcSymbolHistogramStep')>1:raise ValueError('Duplicate histogram option')
    if '--bddcSymbolHistogramStep' in cmd:
        i=cmd.index('--bddcSymbolHistogramStep');del cmd[i:i+2]
    if '--bddcCgConstantSplitExperiment' in cmd:
        raise ValueError('Source command must not already opt into the split experiment')
    split=mode.endswith('-split')
    base_mode=mode.removesuffix('-split')
    uncompressed=base_mode.startswith('uncompressed-')
    if uncompressed:
        solver=base_mode.removeprefix('uncompressed-');bits=16
    else:
        match=re.fullmatch(r'(cg|fcg|fgmres)(4|6|8|12|16|20|24|32)',base_mode)
        if not match:raise ValueError('Unknown mode')
        solver=match[1];bits=int(match[2])
    if solver not in ('cg','fcg','fgmres'):raise ValueError('Unknown solver')
    if split and uncompressed:raise ValueError('Splitting requires compressed transfers')
    settings={'--dir':out,'--bddcHistogram':0,
        '--bddcUseCg':int(solver!='fgmres'), '--bddcFlexibleCg':int(solver=='fcg'),
        '--bddcFgmres':int(solver=='fgmres'), '--bddcCgOriginalResidual':int(solver!='fgmres'),
        '--bddcConvergenceReport':int(solver!='fgmres'),
        '--bddcProlongateConstantSplit':int(split)}
    settings['--bddcCompression']=int(not uncompressed)
    for flag in ('--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                 '--bddcZigzag','--bddcBitlength','--bddcHuffman'):
        settings[flag]=int(not uncompressed and pipeline=="full")
    for flag in ('--bddcCompressionBits','--bddcRestrictCompressionBits','--bddcProlongateCompressionBits'):
        settings[flag]=bits
    for flag,new in settings.items():set_value(flag,new)
    # Only these split candidates opt into the new application validation gate.
    # Baselines and FGMRES retain their previous commands and default flag=0.
    if split and solver in ('cg','fcg'):
        cmd+=['--bddcCgConstantSplitExperiment','1']
    return cmd


def experiment_modes(solvers,bits,split):
    suffixes={'off':('',),'on':('-split',),'both':('','-split')}[split]
    return ['uncompressed-'+s for s in solvers]+[
        f'{s}{b}{suffix}' for b in bits for s in solvers for suffix in suffixes]


def residual_tail(out):
    # Trace files remain useful when MPI abort truncates stdout; do not turn a
    # partial trace into a completed-run or successful-convergence claim.
    paths=list(out.glob('bddc-fgmres-step1.csv')) or list(out.glob('bddc-convergence-step1.csv'))
    if not paths:return None
    points=[]
    with paths[0].open() as f:
        for row in csv.DictReader(f):
            raw=row.get('true_relative_residual',row.get('relative_residual'))
            try:relative=float(raw);iteration=int(row['iteration'])
            except (TypeError,ValueError,KeyError):continue
            if math.isfinite(relative) and relative>=0:points.append((iteration,relative))
    if not points:return None
    best=min(points,key=lambda x:x[1])
    return {'file':paths[0].name,'rows':len(points),'last_iteration':points[-1][0],
            'last_true_relative_residual':points[-1][1],
            'best_iteration':best[0],'best_true_relative_residual':best[1]}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=('10cells','pepe','weak_8_h16','strong_h16'),default='10cells')
    p.add_argument('--reference-campaign',type=Path,help='Override the saved reference campaign selected for this case')
    p.add_argument('--source-run',type=Path,help='Optional Pepe CG source command; other cases use their saved reference campaigns')
    p.add_argument('--repo',type=Path,default=Path.cwd())
    p.add_argument('--kaskade',type=Path,default=Path('/kaskade7'))
    p.add_argument('--solvers',nargs='+',choices=('cg','fcg','fgmres'),default=['cg','fcg'])
    p.add_argument('--bits',nargs='+',type=int,choices=(4,6,8,12,16,20,24,32),default=[16,12],
                   help='Uniform R/P bit depths, in run order (default: 16 12)')
    p.add_argument('--pipeline',choices=('full','quantization-only'),default='full',
                   help='Full coding pipeline or affine quantization with all transform/coding stages off')
    p.add_argument('--split',choices=('off','on','both'),default='both',
                   help='Constant-split prolongation selection (default: paired off/on)')
    p.add_argument('--validate-only',action='store_true',help='Check arguments without building or running experiments')
    args=p.parse_args();repo=args.repo.resolve()
    if args.case!='pepe' and args.source_run:p.error('--source-run is only for Pepe')
    reference=reference_directory(args.case,args.reference_campaign)
    source=(reference/'full-cg16' if args.case!='pepe' else
            (args.source_run or Path('output/emi-step-histograms-20261104T143129Z-89/pepe-cg-r16-p16-off')).resolve())
    if len(set(args.solvers))!=len(args.solvers):p.error('Duplicate solvers')
    if len(set(args.bits))!=len(args.bits):p.error('Duplicate bit depths')
    if args.validate_only:return 0
    sys.path.insert(0,str(repo/'tests'))
    from collect_emi_scaling_smoke import collect
    from run_emi_fgmres_check import saved_control
    from analyze_emi_potential_diagnostics import read,error_metrics
    def git(where,*args):return subprocess.check_output(['git','-C',str(where),*args],text=True).strip()
    if git(repo,'branch','--show-current')!='experiment/cg-fcg-constant-split' or git(args.kaskade,'branch','--show-current')!='master':
        p.error('This check requires EMI experiment/cg-fcg-constant-split and Kaskade master; rebuild EMI first')
    saved=shlex.split((source/'command.sh').read_text())
    try:validate_scaling_source(args.case,saved)
    except ValueError as exc:p.error(str(exc))
    source_command=list(saved)
    def value(flag):return saved[saved.index(flag)+1]
    old_controls={}
    control_names={'cg16':'full-cg16','uncompressed-fgmres':'uncompressed-fgmres',
                   'fgmres16':'full-fgmres16'}
    if args.pipeline=='quantization-only':
        control_names={'uncompressed-fgmres':'uncompressed-fgmres'}
    control_names={m:n for m,n in control_names.items()
                   if m in experiment_modes(args.solvers,args.bits,args.split)}
    if args.case!='pepe':
        for mode,name in control_names.items():
            control=reference/name
            previous=shlex.split((control/'command.sh').read_text())
            previous_steps=int(previous[previous.index('--maximumNumberOfTimeSteps')+1])
            old_controls[mode]=saved_control(control,name,previous_steps)
        for flag,wanted in {'--refine':('2' if args.case=='10cells' else '0'),'--sigma_i':'3','--sigma_e':'20',
                            '--C_m':'1','--R':'0.1','--dt':'0.01'}.items():
            if float(value(flag))!=float(wanted):p.error('Unexpected '+args.case+' configuration: '+flag)
    # Preserve the reference mesh, refinement and physics; test its first step.
    for flag,new in {'--maximumNumberOfTimeSteps':'1','--finalTime':value('--dt'),'--vtk':'0'}.items():
        saved[saved.index(flag)+1]=new

    for flag,wanted in {'--bddcUseCg':'1','--bddcFgmres':'0','--bddcFlexibleCg':'0',
                        '--bddcProlongateConstantSplit':'0','--maximumNumberOfTimeSteps':'1',
                        '--bddcRelativeTolerance':'1e-7','--refine':('2' if args.case=='10cells' else '0'),
                        '--bddcRestrictionPasses':'1','--bddcRestrictZeroPreserving':'0',
                        '--bddcProlongateZeroPreserving':'0','--bddcRestrictGroupSize':'0',
                        '--bddcProlongateGroupSize':'0','--bddcFgmresRestart':'30',
                        '--bddcIterations':'300','--bddcCgDiagnostics':'0',
                        '--bddcProlongationDiagnostics':'0','--bddcConditionEstimate':'0',
                        '--bddcConstantCorrection':'0','--masterWorker':'0'}.items():
        if value(flag)!=wanted:p.error('Unexpected source configuration: '+flag)
    model=repo/'emiModel'
    if str(model) not in saved:p.error('Source executable path differs from repository')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    pipeline_tag='' if args.pipeline=='full' else '-quantization-only'
    root=repo/'output'/f'emi-cg-split-check-{args.case}{pipeline_tag}-{stamp}-{os.getpid()}';root.mkdir(exist_ok=False)
    # Match the saved reference runner's BLAS/OpenMP thread settings.
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE',None)
    model_hash=digest(model)
    manifest={'emi_commit':git(repo,'rev-parse','HEAD'),'kaskade_commit':git(args.kaskade,'rev-parse','HEAD'),
              'executable_sha256':model_hash,'source_run':str(source),
              'source_command':source_command,'source_command_sha256':digest(source/'command.sh'),
              'case':args.case,'refinements':int(value('--refine')),
              'historical_controls':{mode:str(reference/name) for mode,name in control_names.items()} if args.case!='pepe' else {},
              'thread_environment':{k:env[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')},
              'tolerance':1e-7,'steps':1,'capture':False,'constant_split':args.split,'pipeline':args.pipeline,
              'solvers':args.solvers,'bits':args.bits,'iteration_limit':300,'fgmres_restart':30,
              'restriction_grid':'affine','prolongation_grid':'affine',
              'communication_scope':'logical interface bytes plus an outer collective payload model; excludes coarse/codebook/setup traffic and MPI framing; not measured network traffic',
              'purpose':'Opt-in CG/FCG split validation; unchanged recurrences/breakdown checks; no convergence or speedup guarantee'}
    manifest['source_sha256']={name:digest(repo/name) for name in
        ('EMI_model.cpp','EMI_bddc_precision.hh','EMI_bddc.hh','tests/run_emi_cg_split_check.py')}
    (root/'emi-working-tree.patch').write_text(git(repo,'diff','HEAD','--','EMI_model.cpp','EMI_bddc_precision.hh'))
    (root/'environment.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Matched '+args.case+' solver check:',root,flush=True)
    print('One step; refine='+value('--refine')+'; tolerance=1e-7; one R pass; affine R/P; pipeline='+args.pipeline+'; split='+args.split,flush=True)
    if args.pipeline=='quantization-only':
        print('Transforms/mapping/bit-length/Huffman OFF. Quantized values use fixed-size transmission slots; bytes are not a packed b-bit stream.',flush=True)
    results={}
    modes=experiment_modes(args.solvers,args.bits,args.split)
    uncompressed_solutions={}
    for mode in modes:
        base_mode=mode.removesuffix('-split')
        solver=base_mode.removeprefix('uncompressed-') if base_mode.startswith('uncompressed-') else re.sub(r'\d+$','',base_mode)
        if not mode.startswith('uncompressed-') and results['uncompressed-'+solver]['errors']:
            print('Skipping '+mode+': its uncompressed control failed.',flush=True)
            continue
        if digest(model)!=model_hash:raise RuntimeError('Executable changed during campaign')
        out=root/mode;out.mkdir();cmd=make_command(saved,out,mode,args.pipeline)
        replay=['env']+[k+'='+v for k,v in manifest['thread_environment'].items()]+cmd
        (out/'command.sh').write_text(shlex.join(replay)+'\n')
        print('Running '+mode+' ...',flush=True)
        with (out/'run.log').open('w') as f:
            rc=subprocess.run(cmd,cwd=repo,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
        log=(out/'run.log').read_text(errors='replace')
        row=collect(log,1,int(value('-np')),rc,1e-7)
        row['pipeline']=args.pipeline
        row['constant_split']=mode.endswith('-split')
        row['solver_converged']=not row['errors']
        if mode.endswith('-split') and solver in ('cg','fcg'):
            row['experiment_marker_present']='BDDC CG/FCG constant-split experiment:' in log
            if row['solver_converged'] and not row['experiment_marker_present']:
                row['errors'].append('Missing explicit CG/FCG split experiment marker; verify rebuilt executable')
        if mode in old_controls:
            row['historical_parity']='NOT_CHECKED'
            if row['solver_converged']:
                try:
                    snapshot=saved_control(out,control_names[mode],1)
                    equal=snapshot==old_controls[mode]
                    row['historical_parity']='PASS' if equal else 'MISMATCH'
                    if not equal:row['errors'].append('Saved first-step solution/residual/FGMRES-trace parity mismatch')
                except (OSError,ValueError) as exc:
                    row['historical_parity']='ERROR';row['errors'].append(str(exc))
        row['breakdown_lines']=re.findall(r'^BDDC breakdown:.*$',log,re.M)
        row['residual_trace']=residual_tail(out)
        outer=re.findall(r'^BDDC outer solver logical payload bytes: (\d+) collectives=\d+',log,re.M)
        row['outer_payload_model_bytes']=int(outer[0]) if len(outer)==1 else None
        row['accounted_payload_bytes']=(row['transmitted_sent_bytes']+row['outer_payload_model_bytes']
            if not row['errors'] and row['transmitted_sent_bytes'] is not None and row['outer_payload_model_bytes'] is not None else None)
        if row['solver_converged']:
            try:
                solution=read(out/'potential-state-step1.csv')['u']
                if mode.startswith('uncompressed-'):uncompressed_solutions[solver]=solution
                elif solver in uncompressed_solutions:
                    row['accuracy_reference']='uncompressed-'+solver
                    row['vs_uncompressed_same_solver']=error_metrics(solution,uncompressed_solutions[solver])
            except (OSError,ValueError) as exc:row['errors'].append('Accuracy comparison: '+str(exc))
        if mode.endswith('-split') and base_mode in results and not row['errors'] and not results[base_mode]['errors']:
            control=results[base_mode]
            row['comparison_with_unsplit']={
                'control':base_mode,
                'iteration_difference':sum(row['iterations_per_step'])-sum(control['iterations_per_step']),
                'interface_saving_percent':100*(1-row['transmitted_sent_bytes']/control['transmitted_sent_bytes']),
                'accounted_saving_percent':100*(1-row['accounted_payload_bytes']/control['accounted_payload_bytes'])
                    if row['accounted_payload_bytes'] is not None and control['accounted_payload_bytes'] else None}
        row['validation_status']='FAIL' if row['errors'] else 'PASS'
        results[mode]=row
        (out/'metrics.json').write_text(json.dumps(row,indent=2)+'\n')
        (root/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
        print(('PASS' if not row['errors'] else 'FAIL')+': iterations='+str(row['iterations_per_step'])+
              '; residuals='+str(row['relative_residuals']),flush=True)
        if 'historical_parity' in row:
            print('  Saved '+args.case+' first-step parity: '+row['historical_parity'],flush=True)
        if row['solver_converged'] and row['errors']:
            print('  Solver converged, but a validation check failed.',flush=True)
        for error in row['errors'][:3]:print('  '+error,flush=True)
        for line in row['breakdown_lines']:print('  '+line,flush=True)
        if row['errors'] and row['residual_trace']:
            t=row['residual_trace']
            print('  Last recorded residual: '+str(t['last_true_relative_residual'])+
                  ' at iteration '+str(t['last_iteration']),flush=True)
        if not row['errors']:
            print('  Interface bytes='+str(row['transmitted_sent_bytes'])+
                  '; outer payload model='+str(row['outer_payload_model_bytes'])+
                  '; solver wrapper seconds='+str(row['solver_wrapper_max_seconds']),flush=True)
        if 'vs_uncompressed_same_solver' in row:
            error=row['vs_uncompressed_same_solver']
            print('  Solution vs '+row['accuracy_reference']+': raw relative L2='+str(error['relative_l2'])+
                  '; centered='+str(error['centered_relative_l2'])+'; mean error='+str(error['mean_error']),flush=True)
        if 'comparison_with_unsplit' in row:
            comparison=row['comparison_with_unsplit']
            print('  Split vs unsplit: iteration difference='+str(comparison['iteration_difference'])+
                  '; interface saving='+str(comparison['interface_saving_percent'])+'%; accounted saving='+
                  str(comparison['accounted_saving_percent'])+'%',flush=True)
    print('Results:',root)
    raise SystemExit(int(any(r['errors'] for r in results.values())))
if __name__=='__main__':main()
