#!/usr/bin/env python3
"""Controlled EMI right-FGMRES(30) experiment with fresh direct/CG16 controls."""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess

from analyze_emi_potential_diagnostics import read, error_metrics
from collect_emi_scaling_smoke import collect
from run_emi_grouped_prolongation_check import make_command as control_command, first_residual
from run_emi_mixed_precision_validation import digest, direct_errors
from run_emi_sd_cg_comparison import set_value, value

MODES=('direct','full-cg16','uncompressed-fgmres','full-fgmres16','full-fgmres8')
SPLIT_MODE='full-fgmres8-split'
SPLIT_BITS={f'full-fgmres{bits}-split':bits for bits in (8,6,4)}
MIXED_BITS={f'full-fgmres-r{r}-p{p}-split':(r,p) for r,p in ((8,6),(6,8),(8,4),(4,8))}
SPLIT_PAIRS={**{m:(b,b) for m,b in SPLIT_BITS.items()},**MIXED_BITS}
CASES=('10cells','weak_8_h16','strong_h16')
INPUT_SUFFIXES={'--input':'_mesh.vtu','--extra_set':'_list_extracellular.txt',
                '--intra_set':'_list_intracellular.txt','--excited':'_early_excited.txt'}
SPLIT_MARKER='BDDC constant-split prolongation: float32 c; 8-byte header included; existing setup weights'
MARKER=('BDDC FGMRES: root-owned right preconditioning; restart=30; one restriction pass; '
        'fresh RHS feedback per application; no CG update')


def make_command(saved, mode, out, steps=1, case="10cells", mesh_dir=None, four_bit_limit=300):
    if four_bit_limit not in (300,500,1000):raise ValueError('Unsupported controlled 4-bit limit')
    if steps not in (1,10):raise ValueError("Controlled validation supports 1 or 10 steps")
    if case not in CASES:raise ValueError('Unknown validation case')
    if case!='10cells' and mesh_dir is None:raise ValueError('Scaling cases require a mesh directory')
    if mode not in MODES and mode not in SPLIT_PAIRS: raise ValueError('Unknown mode')
    cmd=control_command(saved,'direct' if mode=='direct' else 'full-cg16',out)
    fgmres='fgmres' in mode
    changes={'--maximumNumberOfTimeSteps':steps,
             '--finalTime':str(Decimal(value(saved,'--dt'))*steps),
             '--bddcFgmres':int(fgmres),'--bddcFgmresRestart':30,
             '--bddcPreconditionerProbe':0,'--bddcProlongateConstantSplit':int(mode in SPLIT_PAIRS)}
    if fgmres:
        compressed=mode!='uncompressed-fgmres';bits=SPLIT_BITS.get(mode,8 if mode=='full-fgmres8' or mode in MIXED_BITS else 16)
        rbits,pbits=SPLIT_PAIRS.get(mode,(bits,bits))
        changes.update({'--bddcUseCg':0,'--bddcCgOriginalResidual':0,
                        '--bddcFlexibleCg':0,'--bddcCompression':int(compressed),
                        '--bddcCompressionBits':bits,'--bddcRestrictCompressionBits':rbits,
                        '--bddcProlongateCompressionBits':pbits,
                        '--bddcRestrictZeroPreserving':int(rbits in (8,6,4)),
                        '--bddcProlongateZeroPreserving':0,'--bddcProlongateGroupSize':0})
        for flag in ('--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                     '--bddcZigzag','--bddcBitlength','--bddcHuffman'):
            changes[flag]=int(compressed)
    if case!='10cells':
        # Reuse the established weak/strong inputs at refine=0; retain all physics
        # and solver/codec settings from the validated 10cells template.
        for flag,wanted in {'--sigma_i':'3','--sigma_e':'20','--C_m':'1',
                            '--R':'0.1','--dt':'0.01'}.items():
            if float(value(saved,flag))!=float(wanted):
                raise ValueError('Scaling physical parameters differ: '+flag)
        base=Path(mesh_dir)/(case+'_kaskade')
        changes.update({flag:str(base)+suffix for flag,suffix in INPUT_SUFFIXES.items()})
        changes['--refine']=0
    for flag,new in changes.items():set_value(cmd,flag,new)
    if mode=='full-fgmres4-split':set_value(cmd,'--bddcIterations',four_bit_limit)
    return cmd


def check_trace_prefix(previous, current):
    """A larger cap must reproduce every saved observation before continuing."""
    with previous.open() as f:old=list(csv.DictReader(f))
    with current.open() as f:new=list(csv.DictReader(f))
    if not old or len(new)<len(old) or new[:len(old)]!=old:
        raise ValueError('Extended-cap residual prefix differs from saved run')
    return len(old)


def inspect_accounting(text, mode, metrics):
    rows=re.findall(r'^BDDC outer solver logical payload bytes: (\d+) collectives=(\d+)\s*$',text,re.M)
    if len(rows)!=1:raise ValueError('Expected one outer communication ledger')
    payload,calls=map(int,rows[0]); n=metrics['dofs']
    iterations=metrics['iterations_per_step']; k=sum(iterations); steps=len(iterations)
    peers=metrics['mpi_ranks']-1
    if mode=='full-cg16':
        expected=peers*(k+steps)*(2*n*8+16); expected_calls=2*(k+steps)
    else:
        expected=peers*(k*(2*n*8+12)+steps*(4+56+n*8+16));expected_calls=4*k+4*steps
    if payload!=expected or calls!=expected_calls:
        raise ValueError(f'Outer communication ledger mismatch: {(payload,calls)} vs {(expected,expected_calls)}')
    return dict(outer_logical_bytes=payload,outer_collectives=calls,
                interface_logical_bytes=metrics['transmitted_sent_bytes'],
                accounted_logical_bytes=payload+metrics['transmitted_sent_bytes'])


def inspect_fgmres(folder,text,metrics):
    rows=re.findall(r'^BDDC FGMRES statistics step (\d+): applications=(\d+) cycles=(\d+) '
                    r'root_matvecs=(\d+) root_peak_basis_bytes=(\d+)\s*$',text,re.M)
    iterations=metrics['iterations_per_step']
    if [int(r[0]) for r in rows]!=list(range(1,len(iterations)+1)):
        raise ValueError('Missing, duplicate or unordered FGMRES statistics')
    per_step=[]
    for row,expected in zip(rows,iterations):
        step,applications,cycles,matvecs,basis=map(int,row)
        with (folder/f'bddc-fgmres-step{step}.csv').open() as f:trace=list(csv.DictReader(f))
        if applications!=expected or [int(r['iteration']) for r in trace]!=list(range(applications+1)):
            raise ValueError(f'FGMRES application/CSV iteration mismatch at step {step}')
        for observation in trace:
            for key in ('estimated_absolute_residual','true_absolute_residual','true_relative_residual'):
                if not math.isfinite(float(observation[key])) or float(observation[key])<0:
                    raise ValueError('Nonfinite/negative FGMRES residual')
        if float(trace[-1]['true_relative_residual'])>1e-7:
            raise ValueError(f'FGMRES CSV failed true tolerance at step {step}')
        per_step.append(dict(step=step,applications=applications,cycles=cycles,
                             root_matvecs=matvecs,root_peak_basis_bytes=basis))
    return dict(applications=sum(iterations),cycles=sum(r['cycles'] for r in per_step),
                root_matvecs=sum(r['root_matvecs'] for r in per_step),
                root_peak_basis_bytes=max(r['root_peak_basis_bytes'] for r in per_step),
                per_step=per_step)


def saved_control(folder, mode, steps=1):
    """Validate all requested steps, then return the first-step parity snapshot."""
    text=(folder/'run.log').read_text(errors='replace')
    errors=direct_errors(text,steps,0) if mode=='direct' else collect(text,steps,8,0,1e-7)['errors']
    if errors:raise ValueError('Saved control is not a completed run: '+str(folder))
    result={'initial':read(folder/'potential-state-step0.csv')['u'],
            'final':read(folder/'potential-state-step1.csv')['u']}
    if mode!='direct':result['residual']=first_residual(text)
    if 'fgmres' in mode:
        with (folder/'bddc-fgmres-step1.csv').open() as f:result['trace']=list(csv.DictReader(f))
    return result


def validate_mixed_scaling_reference(ref):
    """Require the completed ten-step 10cells R6/P8 run before new meshes."""
    mode='full-fgmres-r6-p8-split';folder=ref/'10cells'/mode
    cmd=shlex.split((folder/'command.sh').read_text())
    expected={'--maximumNumberOfTimeSteps':10,'--bddcRestrictCompressionBits':6,
              '--bddcProlongateCompressionBits':8,'--bddcRestrictionPasses':1,
              '--bddcFgmres':1,'--bddcFgmresRestart':30,'--bddcIterations':300,
              '--bddcRelativeTolerance':1e-7,'--bddcProlongateConstantSplit':1,
              '--bddcRestrictZeroPreserving':1,'--bddcProlongateZeroPreserving':0,
              '--bddcConstantCorrection':0,'--bddcCompression':1}
    expected.update({flag:1 for flag in ('--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                                      '--bddcZigzag','--bddcBitlength','--bddcHuffman')})
    for flag,wanted in expected.items():
        if float(value(cmd,flag))!=wanted:
            raise ValueError('R6/P8 ten-step reference setting differs: '+flag)
    saved_control(folder,mode,10)
    text=(folder/'run.log').read_text(errors='replace')
    inspect_fgmres(folder,text,collect(text,10,8,0,1e-7))
    return dict(path=str(folder),command_sha256=digest(folder/'command.sh'),
                log_sha256=digest(folder/'run.log'),steps=10)


def residual_progress(folder, step=1):
    """Retain partial progress on failure; do not turn it into a convergence claim."""
    path=folder/f'bddc-fgmres-step{step}.csv'
    if not path.exists():return None
    with path.open() as f:rows=list(csv.DictReader(f))
    valid=[]
    for r in rows:
        try:
            it=int(r['iteration']);v=float(r['true_relative_residual'])
            if math.isfinite(v) and v>=0:valid.append((it,v))
        except (KeyError,TypeError,ValueError):continue
    if not valid:return None
    best=min(valid,key=lambda r:r[1])
    return dict(last_iteration=valid[-1][0],last_true_relative=valid[-1][1],
                best_iteration=best[0],best_true_relative=best[1])


def save(root,records,steps=1,case="10cells"):
    lines=[f'{case}, {steps} time step(s); FGMRES(30), 8 ranks, tolerance 1e-7, one restriction pass.',
           'mode                 status  calls   interface B       outer B     accounted B   saved vs CG16  raw rel L2']
    for mode,row in records.items():
        for metric in ('interface_logical_bytes','accounted_logical_bytes'):
            row['saved_'+metric+'_percent']={}
            for control in ('full-cg16','full-fgmres16','uncompressed-fgmres',SPLIT_MODE):
                b=records.get(control,{})
                saving=(100*(1-row[metric]/b[metric]) if row['status']=='PASS' and
                        b.get('status')=='PASS' and metric in row and metric in b and b[metric] else None)
                row['saved_'+metric+'_percent'][control]=saving
        row['beats_both_accounted_baselines']=(mode in ('full-fgmres8',*SPLIT_PAIRS) and all(
            row['saved_accounted_logical_bytes_percent'][control] is not None and
            row['saved_accounted_logical_bytes_percent'][control]>0
            for control in ('full-cg16','full-fgmres16')))
        it=sum(row.get('metrics',{}).get('iterations_per_step',[])) if mode!='direct' else '--'
        counts=[str(row.get(k,'--')) for k in ('interface_logical_bytes','outer_logical_bytes','accounted_logical_bytes')]
        s=row['saved_accounted_logical_bytes_percent']['full-cg16']
        e=row.get('full_precision_vs_direct',{}).get('relative_l2')
        lines.append(f'{mode:20s} {row["status"]:6s} {str(it):>5s} {counts[0]:>13s} {counts[1]:>13s} '
                     f'{counts[2]:>15s} {str(round(s,2)) if s is not None else "--":>15s} '
                     f'{f"{e:.4e}" if e is not None else "--":>11s}')
        lines.extend('  '+e for e in row['errors'])
        metrics=row.get('metrics',{})
        if row['status']=='PASS' and metrics:
            lines.append(f"  Iterations per step: {metrics['iterations_per_step']}; "
                         f"solver wrapper max={metrics['solver_wrapper_max_seconds']} s; "
                         f"application wall={metrics['wall_seconds']} s (validation timing)")
        err=row.get('full_precision_vs_direct')
        if err:
            lines.append(f"  centered rel L2={err['centered_relative_l2']:.6e}; max abs={err['max_absolute_error']:.6e}; mean error={err['mean_error']:.6e}")
        if row.get('saved_control_parity'):lines.append('  Saved first-step residual/solution parity: PASS')
        progress=row.get('residual_progress')
        if progress and row['status']!='PASS':
            lines.append(f"  Last true rel={progress['last_true_relative']:.6e} at {progress['last_iteration']}; best={progress['best_true_relative']:.6e} at {progress['best_iteration']}")
        if mode in ('full-fgmres8',*SPLIT_PAIRS) and row['status']=='PASS':
            for control in ('full-cg16','full-fgmres16','uncompressed-fgmres',SPLIT_MODE):
                iv=row['saved_interface_logical_bytes_percent'][control]
                av=row['saved_accounted_logical_bytes_percent'][control]
                if iv is not None and av is not None:
                    lines.append(f'  Saved vs {control}: interface={iv:.2f}%; accounted={av:.2f}%')
    lines+=['','Fresh direct, uncompressed FGMRES, CG16 and FGMRES16 controls; no controls are reused.',
            'FGMRES stores distinct preconditioned directions; no CG recurrence or rho positivity guard.',
            'Each BDDC application starts from zero with fresh RHS feedback and includes the interior solve.',
            'Uniform split modes use the named bit depth in BOTH directions; mixed modes explicitly name R and P bits.',
            'All split modes use exact-zero restriction, affine prolongation and float32 scale/minimum.',
            'Suffix -split adds the float32 producer constant before-weighting split; its 8-byte packet headers are already included in interface bytes.',
            'No precision fallback, extra restriction pass, post-solve constant correction or tolerance relaxation.',
            'PASS requires true convergence and finite raw/centered errors, with no hidden accuracy threshold.',
            'Accounted bytes = logical interface deliveries + outer solver collective payload model.',
            'Outer model: Bcast/Reduce=(P-1)*payload; Allreduce=2*(P-1)*payload. Includes validation assembly/norms in CG.',
            'Coarse-grid communication, codebook training, internal BDDC control reductions, MPI framing and setup traffic',
            'are NOT included in this sum. It is not a measurement of all network bytes.',
            'Root owns Arnoldi/matvecs: zero distributed orthogonalization/matvec messages, but centralized CPU/memory.',
            'Basis memory, matvec counts, solver and end-to-end times are saved; this is not a scaling benchmark.',
            f'Tested duration: {steps} time step(s); previous campaigns are untouched.']
    per_step=[]
    for mode,row in records.items():
        for item in row.get('per_step_errors',[]):
            per_step.append(dict(mode=mode,**item))
    if per_step:
        lines+=['', 'Per-step comparisons with fresh direct (centering is diagnostic only):',
                'mode                  step calls true rel      raw rel L2    centered L2    mean error']
        for r in per_step:
            lines.append(f"{r['mode']:22s} {r['step']:3d} {r['iterations']:5d} "
                         f"{r['true_relative_residual']:.6e} {r['relative_l2']:.6e} "
                         f"{r['centered_relative_l2']:.6e} {r['mean_error']:.6e}")
        with (root/'per-step-errors.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(per_step[0]))
            writer.writeheader();writer.writerows(per_step)
    (root/'summary.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--reference-campaign',type=Path,required=True)
    p.add_argument('--output',type=Path)
    p.add_argument('--constant-split',action='store_true',
                   help='test split 8/8 instead of repeating failed original 8/8; requires saved one-step FGMRES controls')
    p.add_argument('--steps',type=int,choices=(1,10),default=1,
                   help='10 requires a successful constant-split one-step reference campaign')
    p.add_argument('--case',choices=CASES,default='10cells',
                   help='weak/strong cases use fresh controls and the successful ten-step 10cells reference')
    p.add_argument('--mesh-dir',type=Path,default=Path('../input_emi_mesh_blocks/weak_strong_scaling'))
    p.add_argument('--bits',type=int,nargs='+',choices=(8,6,4),
                   help='10cells split screening; 6/6 also supports ten steps; fresh 8/8 control precedes candidates')
    p.add_argument('--pairs',nargs='+',choices=('8/6','6/8','8/4','4/8'),
                   help='directional screen; ten-step 6/8 uses a one-step screen for 10cells or completed ten-step 10cells reference for weak/strong')
    p.add_argument('--four-bit-limit',type=int,choices=(300,500,1000),default=300,
                   help='change ONLY the 4/4 one-step cap; controls remain at 300')
    p.add_argument('--screening-campaign',type=Path,
                   help='saved one-step candidate campaign for first-step/prefix parity')
    args=p.parse_args()
    if args.pairs and (args.bits is not None or not args.constant_split or args.four_bit_limit!=300):
        p.error('--pairs requires --constant-split; no --bits or cap change')
    if args.pairs and args.steps==10 and args.pairs!=['6/8']:
        p.error('Mixed ten-step validation requires --pairs 6/8')
    if args.pairs and args.steps==10 and args.case=='10cells' and not args.screening_campaign:
        p.error('10cells mixed ten-step validation requires its successful --screening-campaign')
    if args.pairs and len(set(args.pairs))!=len(args.pairs):p.error('Use distinct directional pairs')
    if args.bits is not None and (not args.constant_split or args.case!='10cells'):
        p.error('--bits screening requires --constant-split --case 10cells')
    if args.bits and args.steps==10 and 4 in args.bits:p.error('4-bit validation remains one step')
    if args.four_bit_limit!=300 and (not args.bits or 4 not in args.bits or args.steps!=1 or not args.screening_campaign):
        p.error('Extended cap requires --bits 4 --steps 1 and --screening-campaign')
    if args.bits and 6 in args.bits and args.steps==10 and not args.screening_campaign:
        p.error('6-bit ten-step validation requires its successful --screening-campaign')
    if args.bits is not None and len(set(args.bits))!=len(args.bits):p.error('Use distinct bit depths')
    if args.steps==10 and not args.constant_split:p.error('--steps 10 requires --constant-split')
    if args.case!='10cells' and (args.steps!=10 or not args.constant_split):
        p.error('Weak/strong validation requires --constant-split --steps 10')
    repo=args.repo.resolve();ref=args.reference_campaign.resolve()
    mesh_dir=(repo/args.mesh_dir).resolve()
    if args.case!='10cells':
        for suffix in INPUT_SUFFIXES.values():
            path=mesh_dir/(args.case+'_kaskade'+suffix)
            if not path.is_file():p.error('Missing input '+str(path))
    source=ref/'10cells/full-cg16';saved=shlex.split((source/'command.sh').read_text())
    previous=(source/'run.log').read_text(errors='replace')
    source_steps=int(value(saved,'--maximumNumberOfTimeSteps'))
    if args.case!='10cells' and source_steps!=10:
        p.error('Weak/strong validation requires a completed ten-step 10cells reference')
    if collect(previous,source_steps,8,0,1e-7)['errors']:
        p.error('Source CG16 control did not pass')
    model=repo/'emiModel';sha=digest(model)
    if str(model) not in saved:p.error('Source executable path differs')
    modes=MODES[:-1]+(SPLIT_MODE,) if args.constant_split else MODES
    if args.bits is not None:
        modes=MODES[:-1]+tuple(m for m,bits in SPLIT_BITS.items() if bits==8 or bits in args.bits)
    candidate_modes=set(SPLIT_BITS)-{SPLIT_MODE} if args.bits is not None else {modes[-1]}
    if args.pairs:
        selected={tuple(map(int,pair.split('/'))) for pair in args.pairs}
        candidates=tuple(m for m,bits in MIXED_BITS.items() if bits in selected)
        modes=MODES[:-1]+(SPLIT_MODE,)+candidates;candidate_modes=set(candidates)
    old_controls={}
    if args.constant_split:
        for mode in (MODES[:-1]+(SPLIT_MODE,) if args.steps==10 or args.bits is not None or args.pairs else MODES[:-1]):
            old_folder=ref/'10cells'/mode
            old_controls[mode]=saved_control(old_folder,mode,source_steps)
            if args.case!='10cells' and 'fgmres' in mode:
                old_text=(old_folder/'run.log').read_text(errors='replace')
                inspect_fgmres(old_folder,old_text,collect(old_text,source_steps,8,0,1e-7))
    screening=args.screening_campaign.resolve() if args.screening_campaign else None
    if screening and args.bits and 6 in args.bits and args.steps==10:
        old_controls['full-fgmres6-split']=saved_control(screening/'10cells/full-fgmres6-split','full-fgmres6-split',1)
    mixed_scaling_reference=None
    if args.pairs and args.case!='10cells':
        mixed_scaling_reference=validate_mixed_scaling_reference(ref)
    if args.pairs and args.steps==10 and args.case=='10cells':
        mode='full-fgmres-r6-p8-split'
        old_controls[mode]=saved_control(screening/'10cells'/mode,mode,1)
    if args.four_bit_limit!=300:
        if not (screening/'10cells/full-fgmres4-split/bddc-fgmres-step1.csv').is_file():
            p.error('Missing saved 4-bit trace')
    for mode in modes:make_command(saved,mode,Path('validation-only'),args.steps,args.case,mesh_dir,args.four_bit_limit)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=(args.output or repo/'output'/f'emi-fgmres-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True,exist_ok=False)
    (root/'environment.json').write_text(json.dumps(dict(executable_sha256=sha,source_command=saved,
        source_log_sha256=digest(source/'run.log'),reference_campaign=str(ref),reference_steps=source_steps,steps=args.steps,ranks=8,
        case=args.case,mesh_dir=str(mesh_dir) if args.case!='10cells' else None,
        saved_first_step_parity_required=args.case=='10cells',
        mixed_scaling_reference=mixed_scaling_reference,
        restart=30,tolerance=1e-7,restriction_passes=1,controls_reused=False,
        constant_split=args.constant_split,requested_bits=args.bits,
        requested_pairs=args.pairs,directional_bits={m:list(SPLIT_PAIRS[m]) for m in modes if m in SPLIT_PAIRS},
        four_bit_limit=args.four_bit_limit,screening_campaign=str(screening) if screening else None,
        split_bit_depths=[SPLIT_BITS[m] for m in modes if m in SPLIT_BITS],
        metadata='float32 scale/minimum; split candidate adds float32 c and uint32 tag',
        layout='root-owned Arnoldi; distributed BDDC'),indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1');env.pop('KASKADE_RESTRICTION_CAPTURE',None)
    records={};reference={};initial_reference=None
    print(f'FGMRES {args.case} {args.steps}-step experiment:',root,flush=True)
    for mode in modes:
        if digest(model)!=sha:raise RuntimeError('Executable changed during experiment')
        folder=root/args.case/mode;folder.mkdir(parents=True)
        cmd=make_command(saved,mode,folder,args.steps,args.case,mesh_dir,args.four_bit_limit);(folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print('Running '+mode+' ...',flush=True)
        with (folder/'run.log').open('w') as f:
            rc=subprocess.run(cmd,cwd=repo,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
        text=(folder/'run.log').read_text(errors='replace')
        metrics={} if mode=='direct' else collect(text,args.steps,8,rc,1e-7)
        errors=direct_errors(text,args.steps,rc) if mode=='direct' else list(metrics['errors'])
        row=dict(status='FAIL',errors=errors,return_code=rc,metrics=metrics,
                 restriction_bits=int(value(cmd,'--bddcRestrictCompressionBits')),
                 prolongation_bits=int(value(cmd,'--bddcProlongateCompressionBits')),
                 constant_split=mode in SPLIT_PAIRS)
        if text.count(MARKER)!=int('fgmres' in mode):errors.append('FGMRES marker mismatch')
        if text.count(SPLIT_MARKER)!=int(mode in SPLIT_PAIRS):errors.append('Constant-split marker mismatch')
        if 'fgmres' in mode:
            row['residual_progress_per_step']={str(step):progress for step in range(1,args.steps+1)
                if (progress:=residual_progress(folder,step)) is not None}
            row['residual_progress']=next(reversed(row['residual_progress_per_step'].values()),None)
        if mode=='full-fgmres4-split' and args.four_bit_limit!=300:
            try:
                row['saved_trace_prefix_rows']=check_trace_prefix(
                    screening/'10cells'/mode/'bddc-fgmres-step1.csv',folder/'bddc-fgmres-step1.csv')
            except (OSError,ValueError) as e:errors.append(str(e))
        if not errors:
            try:
                u0=read(folder/'potential-state-step0.csv')['u']
                if mode=='direct':
                    initial_reference=u0
                    reference={step:read(folder/f'potential-state-step{step}.csv')['u']
                               for step in range(1,args.steps+1)}
                    if any(len(u)!=len(u0) for u in reference.values()):
                        raise ValueError('Direct snapshot DOF counts changed')
                else:
                    if u0!=initial_reference:raise ValueError('Initial state differs from direct')
                    row['per_step_errors']=[]
                    for step in range(1,args.steps+1):
                        u=read(folder/f'potential-state-step{step}.csv')['u']
                        if len(u)!=metrics['dofs']:raise ValueError('Full-precision snapshot DOF count mismatch')
                        err=error_metrics(u,reference[step])
                        if any(v is None or not math.isfinite(v) for v in err.values()):
                            raise ValueError(f'Nonfinite solution comparison at step {step}')
                        row['per_step_errors'].append(dict(step=step,
                            iterations=metrics['iterations_per_step'][step-1],
                            true_relative_residual=metrics['relative_residuals'][step-1],**err))
                    row['full_precision_vs_direct']=err
                    row.update(inspect_accounting(text,mode,metrics))
                    if 'fgmres' in mode:row['fgmres']=inspect_fgmres(folder,text,metrics)
                    if mode=='full-cg16' and args.case=='10cells':
                        row['source_first_step_parity']=first_residual(text)==first_residual(previous)
                        if not row['source_first_step_parity']:raise ValueError('CG16 first-step parity failed')
                if mode in old_controls and args.case=='10cells':
                    if saved_control(folder,mode,args.steps)!=old_controls[mode]:
                        raise ValueError('Saved one-step control parity failed')
                    row['saved_control_parity']=True
            except (OSError,ValueError) as e:errors.append(str(e))
        row['status']='FAIL' if errors else 'PASS';records[mode]=row;save(root,records,args.steps,args.case)
        (folder/'metrics.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
        print('  '+row['status'],flush=True)
        if errors and mode not in candidate_modes:
            print('Control failed; stopping before candidate.',flush=True);break
    if digest(model)!=sha:raise RuntimeError('Executable changed during experiment')
    print('\n'+(root/'summary.txt').read_text());print('Results:',root)
    return int(len(records)!=len(modes) or any(r['status']!='PASS' for r in records.values()))

if __name__=='__main__':raise SystemExit(main())
