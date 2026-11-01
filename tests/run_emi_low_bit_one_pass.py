#!/usr/bin/env python3
"""One-pass low-bit restriction ladder, reusing validated controls from the SAME binary."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from collect_emi_scaling_smoke import collect
from run_emi_mixed_precision_validation import digest,correction_metrics,logical_bytes,objective
from run_emi_quantization_diagnostic import comparison
from analyze_emi_potential_diagnostics import read,error_metrics


def value(cmd,flag):
    if cmd.count(flag)!=1: raise ValueError('Expected exactly one '+flag)
    return cmd[cmd.index(flag)+1]


def candidate_command(saved,bits,out,prolongation_bits=16):
    cmd=list(saved)
    required={'--mpi':'1','--sdc':'0','--bddc':'1','--direct':'0',
              '--bddcCompression':'1','--bddcUseCg':'1','--bddcCgOriginalResidual':'1',
              '--bddcRelativeTolerance':'1e-7','--bddcRestrictionPasses':'1',
              '--bddcRestrictGroupSize':'0','--bddcRestrictZeroPreserving':'1',
              '--bddcProlongateCompressionBits':'16','--bddcConstantCorrection':'1',
              '--potentialDiagnostics':'1','--bddcHuffman':'1'}
    for flag,wanted in required.items():
        if value(cmd,flag)!=wanted: raise ValueError('Saved experiment must use '+flag+' '+wanted)
    for flag,new in [('--bddcRestrictCompressionBits',str(bits)),
                     ('--bddcProlongateCompressionBits',str(prolongation_bits)),('--dir',str(out))]:
        value(cmd,flag);cmd[cmd.index(flag)+1]=new
    return cmd


def save(root,records):
    (root/'summary.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    lines=['One restriction pass; corrected rho; true tolerance 1e-7.',
           'case          Rbits Pbits solve/checks  iterations  bytes saved vs16  full-precision relative L2']
    for r in records:
        m=r.get('metrics',{}); iters=sum(m.get('iterations_per_step',[]))
        saving=r.get('saved_vs16_percent')
        err=r.get('full_precision_vs_direct',{}).get('relative_l2')
        lines.append(f"{r['case']:13s} {r['bits']:5d} {r['prolongation_bits']:5d} {r['status']:12s} {iters:10d} "
                     f"{f'{saving:.2f}%' if saving is not None else '--':>17s} "
                     f"{f'{err:.4e}' if err is not None else '--':>26s}")
        for error in r['errors']: lines.append('  '+error)
        if 'saved_prolongation_control' in r:
            paired=r['saved_prolongation_control']
            lines.append(f"  Compared with saved R{r['bits']}/P{paired['prolongation_bits']}: "
                         f"{paired['saved_percent']:.2f}% fewer logical bytes; "
                         f"iterations {paired['iterations']} -> {iters}.")
    lines+=['','PASS means every step converged and finite comparisons completed; it is not a physical-accuracy certificate.',
            'Raw and centered errors are reported, not required to be exactly no worse than 16-bit.',
            'Communication savings are a separate result and include the post-CG correction broadcasts.',
            'Directional bit depths are listed explicitly; equal depths do not change the quantizer policies or float32 metadata.',
            'No additional restriction pass, precision fallback or tolerance relaxation is allowed.',
            'The post-CG scalar correction remains enabled, identically to the saved controls.',
            'Byte counters include intra-rank transfers; timings include diagnostics. Controls were reused, not retimed.']
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--reference-campaign',type=Path,required=True)
    p.add_argument('--cases',nargs='+',default=['10cells'])
    p.add_argument('--bits',nargs='+',type=int,default=[12,10,8])
    p.add_argument('--prolongation-bits',type=int,default=16)
    p.add_argument('--comparison-campaign',type=Path,
                   help='Saved low-bit ladder to compare at the same restriction precision; only prolongation may differ')
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    if len(set(args.bits))!=len(args.bits) or any(b<2 or b>15 for b in args.bits):
        p.error('Choose distinct restriction precisions in 2..15')
    if len(set(args.cases))!=len(args.cases): p.error('Cases must be distinct')
    if not 2<=args.prolongation_bits<=32: p.error('Prolongation precision must be in 2..32')
    repo=args.repo.resolve()
    ref=args.reference_campaign.resolve()
    envinfo=json.loads((ref/'environment.json').read_text())
    previous=json.loads((ref/'summary.json').read_text())
    model=repo/'emiModel'
    if digest(model)!=envinfo['executable_sha256']:
        p.error('Binary differs from the saved controls. Do not rebuild for this ladder.')
    steps=envinfo['steps'];ranks=envinfo['bddc_ranks']
    seed_mode=f"full-zero-r{envinfo['candidate_restriction_bits']}-p16"
    seeds={};baselines={}
    for case in args.cases:
        if case not in previous: p.error('Missing reference case '+case)
        runs=previous[case]['runs']
        for mode in ['direct','uncompressed-cg','full-r16-p16',seed_mode]:
            if runs[mode]['status']!='PASS': p.error('Reference solver checks failed: '+case+'/'+mode)
        seed=shlex.split((ref/case/seed_mode/'command.sh').read_text())
        candidate_command(seed,args.bits[0],Path('validation-only'))
        if str(model) not in seed: p.error('Saved command uses another executable')
        if int(value(seed,'--maximumNumberOfTimeSteps'))!=steps or int(value(seed,'-np'))!=ranks:
            p.error('Saved step/rank metadata mismatch')
        baseline=shlex.split((ref/case/'full-r16-p16/command.sh').read_text())
        for flag in ['--bddcCgOriginalResidual','--bddcConstantCorrection']:
            if value(baseline,flag)!='1': p.error('Baseline must use the same corrected solver')
        if value(baseline,'--bddcRestrictionPasses')!='1':
            p.error('Baseline must use one restriction pass')
        for mode,filename in [('direct','emiLast.vtu'),('uncompressed-cg','emiBDDCLast.vtu')]:
            if not (ref/case/mode/filename).is_file(): p.error('Missing saved reference solution')
        # Read only the new algebraic CSV reference; do not parse meshes here.
        read(ref/case/'direct'/f'potential-state-step{steps}.csv')
        seeds[case]=seed;baselines[case]=runs['full-r16-p16']
    paired_controls={}
    if args.comparison_campaign:
        paired_root=args.comparison_campaign.resolve()
        paired_env=json.loads((paired_root/'environment.json').read_text())
        if paired_env['reference_environment']!=envinfo or Path(paired_env['reference_campaign']).resolve()!=ref:
            p.error('Comparison campaign must use the same saved controls and executable')
        paired_rows=json.loads((paired_root/'summary.json').read_text())
        for case in args.cases:
            for bits in args.bits:
                matches=[r for r in paired_rows if r['case']==case and r['bits']==bits]
                if len(matches)!=1 or matches[0]['status']!='PASS':
                    p.error(f'Missing successful unique comparison R{bits} for {case}')
                control=matches[0]
                pb=control.get('prolongation_bits',paired_env['prolongation_bits'])
                old_cmd=shlex.split((paired_root/case/f'full-zero-r{bits}-p{pb}'/'command.sh').read_text())
                wanted=candidate_command(seeds[case],bits,Path('comparison-only'),args.prolongation_bits)
                for flag in ['--dir','--bddcProlongateCompressionBits']:
                    value(old_cmd,flag)
                    old_cmd[old_cmd.index(flag)+1]=value(wanted,flag)
                if old_cmd!=wanted: p.error('Comparison differs in settings other than prolongation and output')
                if logical_bytes(control)<=0: p.error('Comparison byte count must be positive')
                paired_controls[case,bits]=(control,pb)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=(args.output or repo/'output'/f'emi-low-bit-one-pass-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True,exist_ok=False)
    (root/'environment.json').write_text(json.dumps(dict(
        reference_campaign=str(ref),reference_environment=envinfo,
        cases=args.cases,restriction_bits=args.bits,prolongation_bits=args.prolongation_bits,
        comparison_campaign=str(args.comparison_campaign.resolve()) if args.comparison_campaign else None,
        restriction_passes=1,tolerance=1e-7,controls_reused=True,
        acceptance='Convergence at every step and finite comparisons; bytes and accuracy reported separately',
        preserved_reference_summary=previous),indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE',None)
    records=[]
    print('Low-bit one-pass output:',root,flush=True)
    for case in args.cases:
        reference=read(ref/case/'direct'/f'potential-state-step{steps}.csv')['u']
        for bits in args.bits:
            if digest(model)!=envinfo['executable_sha256']: raise RuntimeError('Binary changed during ladder')
            folder=root/case/f'full-zero-r{bits}-p{args.prolongation_bits}';folder.mkdir(parents=True)
            cmd=candidate_command(seeds[case],bits,folder,args.prolongation_bits)
            (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
            print(f'Running {case}: R{bits}/P{args.prolongation_bits}, ONE restriction pass ...',flush=True)
            with (folder/'run.log').open('w') as stream:
                rc=subprocess.run(cmd,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT).returncode
            text=(folder/'run.log').read_text(errors='replace')
            m=collect(text,steps,ranks,rc,1e-7)
            corrections,extra,errors=correction_metrics(text,steps,ranks)
            row=dict(case=case,bits=bits,prolongation_bits=args.prolongation_bits,status='FAIL',metrics=m,
                     errors=list(m['errors'])+errors,constant_corrections=corrections)
            m['constant_correction_payload_bytes']=extra if not errors else 0
            if not row['errors']:
                try:
                    row['comparison_to_direct']=comparison(repo,folder/'emiBDDCLastCompression.vtu',
                        ref/case/'direct/emiLast.vtu',folder/'compare-direct.log')
                    row['comparison_to_uncompressed']=comparison(repo,folder/'emiBDDCLastCompression.vtu',
                        ref/case/'uncompressed-cg/emiBDDCLast.vtu',folder/'compare-uncompressed.log')
                    u=read(folder/f'potential-state-step{steps}.csv')['u']
                    row['full_precision_vs_direct']=error_metrics(u,reference)
                except (RuntimeError,OSError,ValueError) as error:
                    row['errors'].append(str(error))
            row['status']='FAIL' if row['errors'] else 'PASS'
            if row['status']=='PASS':
                row['saved_vs16_percent']=100*(1-logical_bytes(row)/logical_bytes(baselines[case]))
                row['communication_better_than16']=logical_bytes(row)<logical_bytes(baselines[case])
                row['legacy_strict_objective_failures']=objective(baselines[case],row)
                if (case,bits) in paired_controls:
                    control,pb=paired_controls[case,bits]
                    row['saved_prolongation_control']=dict(
                        prolongation_bits=pb,logical_bytes=logical_bytes(control),
                        iterations=sum(control['metrics']['iterations_per_step']),
                        full_precision_vs_direct=control.get('full_precision_vs_direct'),
                        saved_percent=100*(1-logical_bytes(row)/logical_bytes(control)))
            records.append(row)
            (folder/'metrics.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
            save(root,records)
            print('  '+row['status']+'; logs retained; continuing.',flush=True)
    if digest(model)!=envinfo['executable_sha256']: raise RuntimeError('Binary changed during ladder')
    print('\n'+(root/'summary.txt').read_text())
    print('Results:',root)
    return int(any(r['status']!='PASS' for r in records))


if __name__=='__main__': sys.exit(main())
