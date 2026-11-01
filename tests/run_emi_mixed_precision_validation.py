#!/usr/bin/env python3
"""Validate one-pass mixed restriction precision/P16 against R16/P16."""
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
import sys

from collect_emi_scaling_smoke import collect, NUMBER
from run_emi_quantization_diagnostic import comparison

CASES = ('10cells','weak_8_h16','strong_h16')
MODES = ('direct','uncompressed-cg','full-r16-p16','full-zero-r12-p16')


def modes_for(bits):
    return MODES[:3]+(f'full-zero-r{bits}-p16',)


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def logical_bytes(row):
    metrics=row.get('metrics',{})
    payload=metrics.get('transmitted_sent_bytes')
    return None if payload is None else payload+metrics.get('constant_correction_payload_bytes',0)


def objective(baseline,candidate):
    failures=[]
    if baseline.get('status')!='PASS' or candidate.get('status')!='PASS':
        return ['Baseline or candidate did not pass solver/comparison checks']
    b=logical_bytes(baseline)
    c=logical_bytes(candidate)
    if b is None or c is None or not 0<c<b:
        failures.append('Candidate total logical bytes do not beat this matched 16/16 baseline')
    for key in ['relative L2 error','max absolute error']:
        old=baseline.get('comparison_to_direct',{}).get(key)
        new=candidate.get('comparison_to_direct',{}).get(key)
        if old is None or new is None or not math.isfinite(old) or not math.isfinite(new) or not 0<=new<=old:
            failures.append('Candidate final raw '+key+' is missing or worse than baseline')
    return failures


def command(model,launcher,base,refine,steps,ranks,mode,out,
            candidate_original_rho=False,candidate_cg_diagnostics=False,
            matched_original_rho=False,matched_cg_diagnostics=False,candidate_bits=12,
            potential_diagnostics=False,matched_constant_correction=False):
    modes=modes_for(candidate_bits)
    if mode not in modes: raise ValueError('Unknown mode for candidate precision: '+mode)
    direct=mode=='direct'
    compressed=mode in modes[2:]
    candidate=mode==modes[3]
    cmd=launcher+['--allow-run-as-root','-np',str(1 if direct else ranks),
                  '-x','OPENBLAS_NUM_THREADS','-x','OMP_NUM_THREADS',str(model),
                  '--mpi','1','--sdc','0','--refine',str(refine)]
    for flag,suffix in [('--input','_mesh.vtu'),('--extra_set','_list_extracellular.txt'),
                        ('--intra_set','_list_intracellular.txt'),('--excited','_early_excited.txt')]:
        cmd += [flag,str(base)+suffix]
    if direct:
        cmd += ['--directBackend','mumps-spd']
    cmd += ['--sigma_i','3','--sigma_e','20','--C_m','1','--R','0.1',
            '--dt','0.01','--finalTime',format(steps*.01,'.17g'),
            '--maximumNumberOfTimeSteps',str(steps),'--nThreads','1',
            '--quiet','0','--vtk','1','--diagnostics','0','--profile','1',
            '--bddc',str(int(not direct)),'--direct',str(int(direct)),
            '--bddcCompression',str(int(compressed)),
            '--bddcCompressionBits','16','--bddcRestrictCompressionBits',str(candidate_bits) if candidate else '16',
            '--bddcProlongateCompressionBits','16',
            '--bddcRestrictZeroPreserving',str(int(candidate)),
            '--bddcRestrictGroupSize','0','--bddcRestrictionPasses','1',
            '--bddcUseCg','1','--bddcRelativeTolerance','0' if direct else '1e-7','--bddcIterations','300',
            '--masterWorker','0','--bddcSplitInterfaces','0','--bddcInterfaceTypes','7',
            '--bddcMaterialAssembly','0' if direct else '1','--bddcInterfaceAssembly','0' if direct else '2',
            '--bddcExcludeCornersFromCompression','1',
            '--bddcCompressionMinInterfaceSize','0',
            '--bddcCompressionThresholdInterfaceTypes','4',
            '--bddcVerbose','0','--bddcConvergenceReport',str(int(not direct)),
            '--bddcCompressionReport',str(int(not direct)),
            '--bddcConditionEstimate','0','--bddcHistogram','0','--bddcCgDiagnostics',
            str(int((candidate and candidate_cg_diagnostics) or
                    (compressed and matched_cg_diagnostics)))]
    if candidate_original_rho or matched_original_rho:
        cmd += ['--bddcCgOriginalResidual',
                str(int((candidate and candidate_original_rho) or
                        (compressed and matched_original_rho)))]
    for flag in ['--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                 '--bddcZigzag','--bddcBitlength','--bddcHuffman']:
        cmd += [flag,str(int(compressed))]
    if potential_diagnostics: cmd += ['--potentialDiagnostics','1']
    if matched_constant_correction:
        cmd += ['--bddcConstantCorrection',str(int(compressed))]
    return cmd+['--dir',str(out)]


def matched_command(old,new,allow_restriction_bits=False):
    """Allow only explicitly selected experiment differences."""
    def normalized(cmd):
        cmd=list(cmd)
        ignored=['--dir','--bddcCgDiagnostics','--bddcCgOriginalResidual','--potentialDiagnostics',
                 '--bddcConstantCorrection']
        if allow_restriction_bits: ignored.append('--bddcRestrictCompressionBits')
        for flag in ignored:
            if cmd.count(flag)>1: raise ValueError('Duplicate '+flag)
            if flag in cmd:
                i=cmd.index(flag); del cmd[i:i+2]
        return cmd
    return normalized(old)==normalized(new)


def history_parity(previous,current):
    old={p.name:p.read_bytes() for p in previous.glob('bddc-convergence-step*.csv')}
    new={p.name:p.read_bytes() for p in current.glob('bddc-convergence-step*.csv')}
    return bool(old) and old==new


def original_rho(cmd):
    flag='--bddcCgOriginalResidual'
    return flag in cmd and cmd[cmd.index(flag)+1]=='1'


def constant_correction(cmd):
    flag='--bddcConstantCorrection'
    return flag in cmd and cmd[cmd.index(flag)+1]=='1'


def correction_metrics(text,steps,ranks):
    rows=re.findall(rf'^BDDC constant correction step (\d+): delta=({NUMBER}) sum_r=({NUMBER}) '
                    rf'one_A_one=({NUMBER}) pre_rel=({NUMBER}) post_rel=({NUMBER})\s*$',text,re.M)
    payload=re.findall(r'^BDDC constant correction additional broadcast payload bytes: (\d+)\s*$',text,re.M)
    errors=[]
    if [int(r[0]) for r in rows]!=list(range(1,steps+1)):
        errors.append('Missing or out-of-order constant corrections')
    if any(not all(math.isfinite(float(v)) for v in r[1:]) or float(r[3])<=0
           or not 0<=float(r[4])<=1e-7 or not 0<=float(r[5])<=1e-7 for r in rows):
        errors.append('Invalid constant correction or failed residual recheck')
    expected=32*(ranks-1)*steps
    if len(payload)!=1 or int(payload[0])!=expected:
        errors.append('Missing or incorrect constant-correction broadcast payload count')
    return rows,expected,errors


def direct_errors(text,steps,rc):
    rows=re.findall(rf'^Direct global residual step (\d+): abs=({NUMBER}) rel=({NUMBER})\s*$',text,re.M)
    errors=[]
    if rc or not re.search(r'^EMI model completed\s*$',text,re.M):
        errors.append('Direct run did not complete normally')
    if [int(r[0]) for r in rows]!=list(range(1,steps+1)):
        errors.append('Missing, duplicate or out-of-order direct residual summaries')
    if any(not math.isfinite(float(v)) or float(v)<0 for r in rows for v in r[1:]):
        errors.append('Invalid direct residual')
    if any(float(r[2])>1e-7 for r in rows):
        errors.append('Direct residual exceeds 1e-7')
    return errors


def save(root,results,candidate_bits=12):
    modes=modes_for(candidate_bits)
    (root/'summary.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    lines=[f'Full-compression R{candidate_bits}/P16 validation; one pass, tolerance 1e-7.',
           'case          solver checks   objective   baseline bytes  candidate bytes  saved vs16']
    for name,r in results.items():
        b=logical_bytes(r['runs'].get(MODES[2],{}))
        c=logical_bytes(r['runs'].get(modes[3],{}))
        saving=f'{100*(1-c/b):.2f}%' if b and c is not None else '--'
        solver='PASS' if len(r['runs'])==4 and all(x['status']=='PASS' for x in r['runs'].values()) else 'FAIL'
        lines.append(f'{name:13s} {solver:14s} {"PASS" if r["objective_passed"] else "NOT MET":10s} {str(b):>14s} {str(c):>16s} {saving:>11s}')
        for failure in r['failures']: lines.append('  '+failure)
        candidate=r['runs'].get(modes[3],{})
        baseline=r['runs'].get(MODES[2],{})
        if baseline.get('original_residual_rho') and candidate.get('original_residual_rho'):
            lines.append('  Both compressed configurations use original-local-residual CG rho.')
        elif candidate.get('original_residual_rho'):
            lines.append('  Candidate CG rho uses original local residual; control rho remains legacy.')
        if candidate.get('constant_correction'):
            lines.append('  Both compressed runs use post-CG constant correction; bytes include added broadcast deliveries.')
    lines += ['', 'Every time step must meet the true residual tolerance.',
              'Accuracy compares FINAL exported solutions with direct; intermediate-step solution errors are not assessed.',
              'Objective requires fewer total logical bytes AND no worse raw final L2/max error than matched 16/16.',
              'Counters include intra-rank exchanges; timings include setup, validation and output. This is not a scaling benchmark.']
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--cases',nargs='+',choices=CASES,default=list(CASES))
    p.add_argument('--steps',type=int,default=10)
    p.add_argument('--ranks',type=int,default=8)
    p.add_argument('--mesh-dir',type=Path,default=Path('../input_emi_mesh_blocks/weak_strong_scaling'))
    p.add_argument('--output',type=Path)
    p.add_argument('--potential-diagnostics',action='store_true',
                   help='record read-only full-precision DOF states and analyze mean error per step')
    p.add_argument('--matched-constant-correction',action='store_true',
                   help='post-CG constant correction and true-residual recheck in both compressed runs')
    p.add_argument('--candidate-restriction-bits',type=int,default=12,
                   help='zero-preserving candidate restriction bits, 2..32; prolongation remains 16')
    p.add_argument('--candidate-original-rho',action='store_true',
                   help='opt in to original local residual for candidate CG rho only')
    p.add_argument('--candidate-cg-diagnostics',action='store_true',
                   help='observe the candidate; timings include observer overhead')
    p.add_argument('--matched-original-rho',action='store_true',
                   help='use original-local-residual rho for BOTH compressed configurations')
    p.add_argument('--matched-cg-diagnostics',action='store_true',
                   help='observe BOTH compressed configurations; timings include observer overhead')
    p.add_argument('--reference-campaign',type=Path,
                   help='check saved commands and exact histories where scalar mode is unchanged')
    args=p.parse_args()
    if args.matched_constant_correction and not args.matched_original_rho:
        p.error('The controlled constant-correction experiment requires --matched-original-rho')
    if args.steps<1 or args.ranks<1 or len(set(args.cases))!=len(args.cases):
        p.error('Steps/ranks must be positive and cases distinct')
    if not 2<=args.candidate_restriction_bits<=32:
        p.error('Candidate restriction bits must be in 2..32')
    bits=args.candidate_restriction_bits
    modes=modes_for(bits)
    repo=args.repo.resolve(); model=repo/'emiModel'
    if not model.is_file() or not os.access(model,os.X_OK): p.error('Missing executable '+str(model))
    mesh_dir=args.mesh_dir if args.mesh_dir.is_absolute() else repo/args.mesh_dir
    inputs={}
    for name in args.cases:
        base=repo/'input/10Cells3d_10extra' if name=='10cells' else mesh_dir/(name+'_kaskade')
        for suffix in ['_mesh.vtu','_list_extracellular.txt','_list_intracellular.txt','_early_excited.txt']:
            if not Path(str(base)+suffix).is_file(): p.error('Missing input '+str(base)+suffix)
        inputs[name]=base
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE',None)
    launcher=shlex.split(env.get('MPIEXEC','mpirun'))
    if not launcher: p.error('Empty MPIEXEC')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=args.output or repo/'output'/f'emi-mixed-validation-{stamp}-{os.getpid()}'
    if not root.is_absolute(): root=repo/root
    root=root.resolve(); root.mkdir(parents=True,exist_ok=False)
    sha=digest(model)
    metadata=dict(executable_sha256=sha,hostname=os.uname().nodename,steps=args.steps,
                  bddc_ranks=args.ranks,direct_ranks=1,tolerance=1e-7,restriction_passes=1,
                  group_size=0,cg_diagnostics=args.candidate_cg_diagnostics or args.matched_cg_diagnostics,
                  cg_diagnostic_scope=('both compressed' if args.matched_cg_diagnostics else
                                       'candidate only' if args.candidate_cg_diagnostics else 'off'),
                  candidate_original_rho=args.candidate_original_rho,
                  candidate_cg_diagnostics=args.candidate_cg_diagnostics,
                  matched_original_rho=args.matched_original_rho,
                  matched_cg_diagnostics=args.matched_cg_diagnostics,
                  candidate_restriction_bits=bits,
                  potential_diagnostics=args.potential_diagnostics,
                  matched_constant_correction=args.matched_constant_correction,
                  reference_campaign=str(args.reference_campaign.resolve()) if args.reference_campaign else None,
                  cases={n:dict(mesh_base=str(inputs[n]),refine=2 if n=='10cells' else 0) for n in args.cases},
                  physics=dict(sigma_i=3,sigma_e=20,C_m=1,R=.1,dt=.01),
                  baseline='Original R16/P16, full coding',
                  candidate=f'Zero-preserving R{bits} / original P16, full coding',
                  scope='Every-step residuals and final-state raw accuracy; logical bytes, not hardware traffic')
    (root/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    if args.reference_campaign:
        reference=args.reference_campaign.resolve()
        reference_bits=json.loads((reference/'environment.json').read_text()).get('candidate_restriction_bits',12)
        reference_modes=modes_for(reference_bits)
        for name,base in inputs.items():
            for mode,old_mode in zip(modes,reference_modes):
                old=shlex.split((reference/name/old_mode/'command.sh').read_text())
                new=command(model,launcher,base,2 if name=='10cells' else 0,
                            args.steps,args.ranks,mode,root/name/mode,
                            args.candidate_original_rho,args.candidate_cg_diagnostics,
                            args.matched_original_rho,args.matched_cg_diagnostics,bits,args.potential_diagnostics,
                            args.matched_constant_correction)
                if not matched_command(old,new,allow_restriction_bits=mode==modes[3]):
                    p.error('Saved parameters differ for '+name+'/'+mode)
    print('Validation output:',root,flush=True)
    results={}
    for name,base in inputs.items():
        record=dict(runs={},objective_passed=False,failures=[])
        results[name]=record
        for mode in modes:
            if digest(model)!=sha: raise RuntimeError('Executable changed during validation')
            folder=root/name/mode; folder.mkdir(parents=True)
            cmd=command(model,launcher,base,2 if name=='10cells' else 0,
                        args.steps,args.ranks,mode,folder,
                        args.candidate_original_rho,args.candidate_cg_diagnostics,
                        args.matched_original_rho,args.matched_cg_diagnostics,bits,args.potential_diagnostics,
                        args.matched_constant_correction)
            (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
            print(f'Running {name}: {mode}, {args.steps} steps ...',flush=True)
            with (folder/'run.log').open('w') as f:
                rc=subprocess.run(cmd,cwd=repo,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
            text=(folder/'run.log').read_text(errors='replace')
            row=dict(status='FAIL',return_code=rc,errors=[],
                     original_residual_rho=original_rho(cmd),constant_correction=constant_correction(cmd))
            if mode=='direct':
                row['errors']=direct_errors(text,args.steps,rc)
                if not (folder/'emiLast.vtu').is_file(): row['errors'].append('Missing direct final solution')
            else:
                metrics=collect(text,args.steps,args.ranks,rc,1e-7)
                row['metrics']=metrics; row['errors']=list(metrics['errors'])
                if constant_correction(cmd):
                    corrections,extra,errors=correction_metrics(text,args.steps,args.ranks)
                    row['constant_corrections']=corrections
                    row['errors']+=errors
                    metrics['constant_correction_payload_bytes']=extra if not errors else 0
                    metrics['total_logical_bytes_with_correction']=logical_bytes(row)
                if args.reference_campaign:
                    old_mode=reference_modes[modes.index(mode)]
                    saved=shlex.split((reference/name/old_mode/'command.sh').read_text())
                    same_bits=(saved[saved.index('--bddcRestrictCompressionBits')+1]
                               ==cmd[cmd.index('--bddcRestrictCompressionBits')+1])
                    same_scalar=original_rho(saved)==original_rho(cmd)
                    same_correction=constant_correction(saved)==constant_correction(cmd)
                    if same_scalar and same_bits and same_correction:
                        row['saved_history_parity']=history_parity(reference/name/old_mode,folder)
                        if mode in ('uncompressed-cg','full-r16-p16') and not original_rho(cmd):
                            row['legacy_history_parity']=row['saved_history_parity']
                        if not row['saved_history_parity']:
                            row['errors'].append('Saved convergence history changed with the same scalar mode and precision')
                    else:
                        row['history_parity_note']='Full history parity not required: scalar, precision or post-CG correction changed intentionally'
                        if same_scalar and same_bits and not same_correction:
                            old_first=reference/name/old_mode/'bddc-convergence-step1.csv'
                            new_first=folder/'bddc-convergence-step1.csv'
                            row['first_step_cg_parity']=(old_first.is_file() and new_first.is_file()
                                                       and old_first.read_bytes()==new_first.read_bytes())
                            if not row['first_step_cg_parity']:
                                row['errors'].append('First CG trajectory changed before the post-solve correction')
                if cmd[cmd.index('--bddcCgDiagnostics')+1]=='1':
                    diagnostics=sorted(folder.glob('bddc-cg-diagnostics-step*.csv'))
                    row['cg_diagnostic_files']=[str(f) for f in diagnostics]
                    if not diagnostics:
                        row['errors'].append('Requested CG diagnostic CSVs are missing')
                (folder/'metrics.json').write_text(json.dumps(metrics,indent=2,allow_nan=False)+'\n')
                if not row['errors']:
                    solution=folder/('emiBDDCLast.vtu' if mode=='uncompressed-cg' else 'emiBDDCLastCompression.vtu')
                    try:
                        row['comparison_to_direct']=comparison(repo,solution,root/name/'direct/emiLast.vtu',folder/'compare-direct.log')
                        if mode!='uncompressed-cg':
                            row['comparison_to_uncompressed']=comparison(repo,solution,root/name/'uncompressed-cg/emiBDDCLast.vtu',folder/'compare-uncompressed.log')
                    except (RuntimeError,OSError,ValueError) as error: row['errors'].append(str(error))
            if rc:
                exceptions=list(dict.fromkeys(
                    message.strip() for message in re.findall(r'what\(\):\s*([^\r\n]+)',text)))
                row['errors'][:0]=['Application exception: '+message for message in exceptions]
            row['status']='FAIL' if row['errors'] else 'PASS'
            record['runs'][mode]=row
            print('  '+row['status']+(' '+ '; '.join(row['errors']) if row['errors'] else ''),flush=True)
            save(root,results,bits)
            if row['status']!='PASS':
                record['failures'].append(mode+' failed; inspect '+str(folder/'run.log'))
                break
            if row.get('saved_history_parity'):
                print('  Saved convergence CSV parity (same scalar mode): PASS',flush=True)
        if not record['failures']:
            record['failures']=objective(record['runs'][modes[2]],record['runs'][modes[3]])
        record['objective_passed']=not record['failures']
        save(root,results,bits)
    if digest(model)!=sha: raise RuntimeError('Executable changed during validation')
    potential_complete=True
    if args.potential_diagnostics:
        from analyze_emi_potential_diagnostics import analyze
        potential_complete=analyze(root)['complete']
    print('\n'+(root/'summary.txt').read_text())
    print('Results:',root)
    return int(not potential_complete or any(not r['objective_passed'] for r in results.values()))


if __name__=='__main__':
    sys.exit(main())
