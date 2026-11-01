#!/usr/bin/env python3
"""Validate rebuilt standard-CG controls, then test one-pass flexible CG at 16/16 and 12/12."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

from collect_emi_scaling_smoke import collect
from run_emi_low_bit_one_pass import candidate_command,value
from run_emi_mixed_precision_validation import digest,history_parity,correction_metrics,logical_bytes
from run_emi_quantization_diagnostic import comparison
from analyze_emi_potential_diagnostics import read,error_metrics

MARKER='BDDC flexible CG: one-term Polak-Ribiere; original residual; existing two-scalar reduction; no extra transfer pass'
MODES=['uncompressed-cg','full-r16-p16','flex-r16-p16','flex-zero-r12-p12']


def make_command(saved,mode,out):
    """The 12/12 candidate differs from the old 14/16 seed only in bits, solver and output."""
    cmd=list(saved)
    if '--bddcFlexibleCg' in cmd: raise ValueError('Reference must predate the flexible-CG experiment')
    if mode=='flex-zero-r12-p12':
        cmd=candidate_command(cmd,12,out,12)
    else:
        value(cmd,'--dir');cmd[cmd.index('--dir')+1]=str(out)
    if mode.startswith('flex-'):
        if value(cmd,'--bddcCgOriginalResidual')!='1': raise ValueError('Flexible CG needs original rho')
    cmd+=['--bddcFlexibleCg',str(int(mode.startswith('flex-')))]
    return cmd


def save(root,records):
    lines=['One-pass flexible-CG experiment; true relative tolerance 1e-7.',
           'case          mode                 status  iterations  logical bytes   saved vs CG16  saved vs FCG16  final relative L2']
    for case,runs in records.items():
        for mode,r in runs.items():
            success=r['status']=='PASS'
            iterations=sum(r['metrics'].get('iterations_per_step',[]))
            count=logical_bytes(r) if success else None
            savings=[]
            for base in ['full-r16-p16','flex-r16-p16']:
                b=runs.get(base,{})
                savings.append(100*(1-count/logical_bytes(b)) if success and b.get('status')=='PASS' else None)
            r['saved_vs_standard16_percent'],r['saved_vs_flexible16_percent']=savings
            err=r.get('full_precision_vs_direct',{}).get('relative_l2')
            fmt=lambda x: '--' if x is None else f'{x:.2f}%'
            lines.append(f"{case:13s} {mode:20s} {r['status']:6s} {iterations:10d} "
                         f"{str(count) if count is not None else '--':>14s} {fmt(savings[0]):>15s} "
                         f"{fmt(savings[1]):>15s} {f'{err:.4e}' if err is not None else '--':>18s}")
            for error in r['errors']: lines.append('  '+error)
    lines+=['','PASS requires all steps to meet 1e-7 and finite final comparisons; no hidden accuracy threshold.',
            'Standard controls must reproduce saved convergence CSVs; the rebuilt binary is used for all new runs.',
            'Direct solutions are reused from the saved campaign. Raw/centered final errors remain in JSON.',
            'Flexible CG uses one prior direction (Polak-Ribiere), not full-history orthogonalization.',
            'No precision fallback, extra restriction pass or tolerance relaxation. Float32 metadata is unchanged.',
            'Logical byte totals include post-CG correction deliveries and intra-rank transfers.',
            'Flexible CG reuses the existing two-double reduction; no additional solver message or payload.',
            'Timing includes diagnostics and output, so this is not a performance benchmark.']
    (root/'summary.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--reference-campaign',type=Path,required=True)
    p.add_argument('--cases',nargs='+',default=['10cells','strong_h16'])
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    repo=args.repo.resolve();ref=args.reference_campaign.resolve();model=repo/'emiModel'
    oldenv=json.loads((ref/'environment.json').read_text())
    old=json.loads((ref/'summary.json').read_text())
    if len(set(args.cases))!=len(args.cases):p.error('Cases must be distinct')
    steps=oldenv['steps'];ranks=oldenv['bddc_ranks'];sha=digest(model)
    seed_mode=f"full-zero-r{oldenv['candidate_restriction_bits']}-p16"
    commands={}
    for case in args.cases:
        if case not in old:p.error('Missing reference case '+case)
        for mode in ['direct','uncompressed-cg','full-r16-p16',seed_mode]:
            if old[case]['runs'][mode]['status']!='PASS':p.error('Reference solver failed: '+case+'/'+mode)
        for mode in MODES:
            source='full-r16-p16' if mode=='flex-r16-p16' else seed_mode if mode==MODES[-1] else mode
            saved=shlex.split((ref/case/source/'command.sh').read_text())
            cmd=make_command(saved,mode,Path('validation-only'))
            if str(model) not in cmd:p.error('Reference command uses another executable path')
            for flag,want in [('--maximumNumberOfTimeSteps',str(steps)),('-np',str(ranks)),
                              ('--bddcRestrictionPasses','1'),('--bddcRelativeTolerance','1e-7'),
                              ('--bddcConditionEstimate','0'),('--potentialDiagnostics','1')]:
                if value(cmd,flag)!=want:p.error('Incompatible saved flag '+flag)
            if mode!='uncompressed-cg':
                for flag in ['--bddcCgOriginalResidual','--bddcConstantCorrection','--bddcCgDiagnostics']:
                    if value(cmd,flag)!='1':p.error('Expected '+flag+' 1')
            commands[case,mode]=saved
        if not (ref/case/'direct/emiLast.vtu').is_file():p.error('Missing direct solution')
        read(ref/case/'direct'/f'potential-state-step{steps}.csv')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=(args.output or repo/'output'/f'emi-flexible-cg-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True,exist_ok=False)
    (root/'environment.json').write_text(json.dumps(dict(executable_sha256=sha,
        reference_campaign=str(ref),reference_environment=oldenv,preserved_reference_summary=old,
        cases=args.cases,steps=steps,ranks=ranks,tolerance=1e-7,restriction_passes=1,
        solver='one-term truncated flexible CG, Polak-Ribiere; standard-CG controls rerun'),indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE',None)
    records={}
    print('Flexible-CG validation output:',root,flush=True)
    for case in args.cases:
        records[case]={}
        direct=read(ref/case/'direct'/f'potential-state-step{steps}.csv')['u']
        for mode in MODES:
            if digest(model)!=sha:raise RuntimeError('Executable changed during validation')
            folder=root/case/mode;folder.mkdir(parents=True)
            cmd=make_command(commands[case,mode],mode,folder)
            (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
            print(f'Running {case}: {mode}, {steps} steps ...',flush=True)
            with (folder/'run.log').open('w') as stream:
                rc=subprocess.run(cmd,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT).returncode
            text=(folder/'run.log').read_text(errors='replace')
            m=collect(text,steps,ranks,rc,1e-7)
            r=dict(status='FAIL',return_code=rc,metrics=m,errors=list(m['errors']))
            r['errors'][:0]=list(dict.fromkeys(re.findall(r'what\(\):\s*([^\r\n]+)',text)))
            if mode.startswith('flex-') and text.count(MARKER)!=1:
                r['errors'].append('Missing/duplicate flexible-CG mode marker; rebuild the executable')
            if mode!='uncompressed-cg':
                corrections,extra,errors=correction_metrics(text,steps,ranks)
                r['constant_corrections']=corrections;r['errors']+=errors
                m['constant_correction_payload_bytes']=extra if not errors else 0
            if mode in MODES[:2]:
                r['saved_convergence_parity']=history_parity(ref/case/mode,folder)
                if not r['saved_convergence_parity']:r['errors'].append('Standard-CG convergence history changed')
                if not r['errors'] and logical_bytes(r)!=logical_bytes(old[case]['runs'][mode]):
                    r['errors'].append('Standard-CG logical byte count changed')
            if not r['errors']:
                try:
                    filename='emiBDDCLast.vtu' if mode=='uncompressed-cg' else 'emiBDDCLastCompression.vtu'
                    r['comparison_to_direct']=comparison(repo,folder/filename,ref/case/'direct/emiLast.vtu',folder/'compare-direct.log')
                    r['full_precision_vs_direct']=error_metrics(read(folder/f'potential-state-step{steps}.csv')['u'],direct)
                except (ValueError,RuntimeError,OSError) as error:r['errors'].append(str(error))
            r['status']='FAIL' if r['errors'] else 'PASS'
            records[case][mode]=r
            save(root,records)
            (folder/'metrics.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
            print('  '+r['status'],flush=True)
            if r['status']=='FAIL' and mode in MODES[:2]:
                print('  Standard control failed validation; skipping experiments for this case.',flush=True)
                break
    if digest(model)!=sha:raise RuntimeError('Executable changed during validation')
    print('\n'+(root/'summary.txt').read_text())
    print('Results:',root)
    return int(any(len(r)!=len(MODES) or any(v['status']!='PASS' for v in r.values()) for r in records.values()))

if __name__=='__main__':raise SystemExit(main())
