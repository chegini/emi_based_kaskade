#!/usr/bin/env python3
"""Replay saved commands with shadow diagnostics off/on, preserving numerical options."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess

from collect_emi_scaling_smoke import NUMBER
from run_emi_low_bit_one_pass import value
from run_emi_mixed_precision_validation import digest


def replay_command(saved, folder, enabled, constraint_call=0, grouping_call=0):
    cmd=list(saved)
    if value(cmd,'--bddcRelativeTolerance')!='1e-7' or value(cmd,'--bddcRestrictionPasses')!='1':
        raise ValueError('Replay requires tolerance 1e-7 and one restriction pass')
    if value(cmd,'--potentialDiagnostics')!='1':
        raise ValueError('Replay needs potential snapshots for parity')
    value(cmd,'--dir');cmd[cmd.index('--dir')+1]=str(folder)
    flag='--bddcProlongationDiagnostics'
    if flag in cmd:
        if value(cmd,flag)!='0':raise ValueError('Source must have shadow diagnostics disabled')
        cmd[cmd.index(flag)+1]=str(int(enabled))
    else:cmd += [flag,str(int(enabled))]
    flag='--bddcProlongationConstraintCall'
    if flag in cmd:
        if value(cmd,flag)!='0':raise ValueError('Source must have constraint probe disabled')
        cmd[cmd.index(flag)+1]=str(constraint_call if enabled else 0)
    elif constraint_call:
        cmd += [flag,str(constraint_call if enabled else 0)]
    flag='--bddcProlongationGroupingCall'
    if flag in cmd:
        if value(cmd,flag)!='0':raise ValueError('Source must have grouping probe disabled')
        cmd[cmd.index(flag)+1]=str(grouping_call if enabled else 0)
    elif grouping_call:
        cmd += [flag,str(grouping_call if enabled else 0)]
    if grouping_call and (value(cmd,'--bddcCompression')=='0' or
                         ( '--bddcProlongateZeroPreserving' in cmd and
                           value(cmd,'--bddcProlongateZeroPreserving')!='0')):
        raise ValueError('Grouping replay requires compressed ORIGINAL affine prolongation')
    return cmd


def signature(text):
    residuals=re.findall(rf'BDDC global residual step (\d+): abs=({NUMBER}) rel=({NUMBER}) '
                         rf'target=({NUMBER}) status=(\S+) iterations=(\d+)',text)
    failures=re.findall(r'BDDC step (\d+) failed after (\d+) solver calls: '
                        r'(solver breakdown|iteration limit|solver stagnation)',text)
    return dict(residuals=residuals, failures=failures,
                completed=bool(re.search(r'^EMI model completed\s*$',text,re.MULTILINE)))


def snapshots(root):
    return {p.name:digest(p) for p in sorted(root.glob('potential-*.csv'))}


def inspect(root):
    results=[]
    for path in sorted(root.glob('bddc-prolongation-diagnostics-step*.csv')):
        rows=list(csv.DictReader(path.open()))
        if not rows:raise ValueError('Empty diagnostic: '+str(path))
        for row in rows:
            for key in ['rho_before_p','rho_after_p','rho_change_p','p_error_norm','weight_partition_defect']:
                if not math.isfinite(float(row[key])):raise ValueError('Nonfinite '+key)
        last=rows[-1]
        before=float(last['rho_before_p']);after=float(last['rho_after_p'])
        if before>0 and after<=0:
            finding='Prolongation error changes positive true r.z to nonpositive.'
        elif before<=0:
            finding='Exact-prolongation shadow is already nonpositive; inspect restriction/local/coarse response.'
        else:
            finding='Both true r.z values are positive at this observation.'
        result=dict(file=path.name,rows=len(rows),last=last,finding=finding,
                    max_p_error_rel=max(float(r['p_error_rel']) for r in rows if r['p_error_rel']),
                    max_weight_partition_defect=max(float(r['weight_partition_defect']) for r in rows))
        results.append(result)
    if not results:raise ValueError('No prolongation diagnostic CSVs were produced')
    return results


def inspect_constraints(root, call):
    path=root/'bddc-prolongation-constraints-summary.csv'
    rows=list(csv.DictReader(path.open()))
    if len(rows)!=1 or int(rows[0]['solver_call'])!=call:
        raise ValueError('Missing or incorrect target constraint observation')
    row={k:float(v) for k,v in rows[0].items()}
    if not all(math.isfinite(v) for v in row.values()):
        raise ValueError('Nonfinite constraint observation')
    if row['groups']<=0:raise ValueError('No coarse constraints measured')
    groups=list(csv.DictReader((root/'bddc-prolongation-constraints-groups.csv').open()))
    if len(groups)!=int(row['groups']):raise ValueError('Incomplete constraint group table')
    for group in groups:
        if not all(math.isfinite(float(v)) for v in group.values()):
            raise ValueError('Nonfinite group observation')
    worst=sorted(groups,key=lambda g:float(g['rho_error_with_interior']))[:5]
    print('  Constraint shadow at call',call)
    print('    true rho: actual={rho_after_p:.7e}, mean-corrected={rho_mean_corrected:.7e}, '
          'float32-mean-corrected={rho_float32_mean_corrected:.7e}'.format(**row))
    print('    Unattributed rho error:',row['unattributed_rho_error'])
    print('    Remaining mean defect (exact / float32):',row['max_remaining_mean_error'],
          row['max_remaining_float32_mean_error'])
    print('    Estimated extra scalar payload/call (32 / 64 bits):',
          int(row['estimated_payload32_per_call']),int(row['estimated_payload64_per_call']))
    for group in worst:
        print('    group',group['group'],'dofs',group['dofs'],
              'rho error including interior:',group['rho_error_with_interior'],
              'mean component:',group['rho_mean_error_with_interior'])
    print('    Shadow only: globally known mean rounded to float32; not an implemented distributed codec.')
    print('    Payload estimate excludes framing, IDs, reductions and any additional messages.')
    return dict(summary=row,worst_groups=worst)


def inspect_grouping(root, call):
    with (root/'bddc-prolongation-grouping.csv').open() as stream:
        raw=list(csv.DictReader(stream))
    rows=[{k:float(v) for k,v in r.items()} for r in raw]
    if not rows or any(int(r['solver_call'])!=call for r in rows):
        raise ValueError('Missing grouping target observation')
    if any(not math.isfinite(v) for r in rows for v in r.values()):
        raise ValueError('Nonfinite grouping observation')
    keys=[(int(r['bits']),int(r['group_size'])) for r in rows]
    required={(16,0),(8,0),*((8,g) for g in [8,16,32,64,128,256,512])}
    if len(keys)!=len(set(keys)) or not required.issubset(keys):
        raise ValueError('Incomplete or duplicate grouping candidates')
    current=[r for r in rows if r['current_precision']==1]
    if len(current)!=1 or current[0]['relative_difference_vs_actual']>1e-8:
        raise ValueError('Whole-buffer producer reconstruction failed parity (>1e-8 relative)')
    baseline=next(r for r in rows if r['bits']==16 and r['group_size']==0)
    if baseline['codec_valid']!=1 or baseline['modeled_p_bytes']<=0:
        raise ValueError('Missing valid same-state P16 byte model')
    print('  Grouped prolongation screen at call',call)
    print('    bits group   true r.z       metadata   modeled P bytes   P saving vs16')
    for r in rows:
        valid=r['codec_valid']==1
        saving=100*(1-r['modeled_p_bytes']/baseline['modeled_p_bytes']) if valid else None
        r['sampled_p_saved_vs16_percent']=saving
        r['screen_eligible']=bool(valid and r['rho_candidate']>0 and saving>0)
        saving_text=f'{saving:.2f}%' if valid else 'invalid codec'
        print(f"    {int(r['bits']):4} {int(r['group_size']):5} {r['rho_candidate']: .6e} "
              f"{int(r['metadata_bytes']):10} {int(r['modeled_p_bytes']):17} {saving_text:>15}")
    print('    Same frozen correction and current Huffman table; P16 here is a shadow, not a solver run.')
    print('    Group 0 reproduces the original grid; other rows add run counts and float32 scale/min pairs.')
    print('    Modeled bytes include padding/exact corners; exclude restriction, MPI framing and codebook training.')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--case',default='10cells')
    p.add_argument('--modes',nargs='+')
    p.add_argument('--output',type=Path)
    p.add_argument('--constraint-call',type=int,default=0,
                   help='probe coarse means at this solver call of step 1; 0 disables')
    p.add_argument('--grouping-call',type=int,default=0,
                   help='screen grouped affine prolongation at this call of step 1; 0 disables')
    args=p.parse_args()
    if min(args.constraint_call,args.grouping_call)<0:p.error('Probe calls must be nonnegative')
    args.modes=args.modes or (['full-cg16','full-fcg8','full-sd8'] if args.grouping_call else
                            ['uncompressed-cg','full-cg16','full-fcg8','full-sd8'])
    repo=args.repo.resolve();source=args.campaign.resolve();model=repo/'emiModel';sha=digest(model)
    commands={};previous={}
    for mode in args.modes:
        if Path(mode).name!=mode or mode in commands:p.error('Use distinct mode folder names')
        folder=source/args.case/mode
        cmd=shlex.split((folder/'command.sh').read_text())
        if str(model) not in cmd:p.error('Saved command uses a different executable path')
        replay_command(cmd,Path('validation-only'),True,args.constraint_call,args.grouping_call)
        commands[mode]=cmd
        previous[mode]=signature((folder/'run.log').read_text(errors='replace'))
        if not previous[mode]['residuals']:p.error('Source has no usable true-residual summary: '+mode)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=(args.output or repo/'output'/f'emi-prolongation-replay-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True,exist_ok=False)
    (root/'environment.json').write_text(json.dumps(dict(executable_sha256=sha,
        source_campaign=str(source),case=args.case,modes=args.modes,constraint_call=args.constraint_call,
        grouping_call=args.grouping_call,
        purpose='read-only shadow, no algorithm changes; off/on parity'),indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE',None)
    records={};bad=False
    print('Prolongation replay:',root,flush=True)
    for mode,saved in commands.items():
        record={};records[mode]=record
        for label,enabled in [('off',False),('on',True)]:
            if digest(model)!=sha:raise RuntimeError('Executable changed during replay')
            folder=root/mode/label;folder.mkdir(parents=True)
            cmd=replay_command(saved,folder,enabled,args.constraint_call,args.grouping_call)
            (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
            print(f'Running {mode}: diagnostic {label} ...',flush=True)
            with (folder/'run.log').open('w') as f:
                rc=subprocess.run(cmd,cwd=repo,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
            text=(folder/'run.log').read_text(errors='replace')
            record[label]=dict(return_code=rc,signature=signature(text),snapshots=snapshots(folder))
            print('  application exit:',rc,flush=True)
        record['source_outcome_parity']=record['off']['signature']==previous[mode]
        record['observer_parity']=(record['off']['signature']==record['on']['signature'] and
            bool(record['off']['signature']['residuals']) and bool(record['off']['snapshots']) and
            record['off']['snapshots']==record['on']['snapshots'] and
            (record['off']['return_code']==0)==(record['on']['return_code']==0))
        record['errors']=[]
        if not record['observer_parity']:record['errors'].append('Observer outcome/snapshot parity failed')
        if not record['source_outcome_parity']:record['errors'].append('Rebuilt off run differs from saved outcome')
        try:
            record['diagnostics']=inspect(root/mode/'on')
            if args.constraint_call:
                record['constraints']=inspect_constraints(root/mode/'on',args.constraint_call)
            if args.grouping_call:
                record['grouping']=inspect_grouping(root/mode/'on',args.grouping_call)
            if any(r['max_weight_partition_defect']>1e-10 for r in record['diagnostics']):
                record['errors'].append('Transfer weights do not form the expected partition')
            if mode.startswith('uncompressed') and any(r['max_p_error_rel']>1e-8 for r in record['diagnostics']):
                record['errors'].append('Uncompressed shadow differs materially from actual correction (>1e-8 relative)')
        except (ValueError,OSError) as error:record['errors'].append(str(error))
        print('  Observer/source parity:',record['observer_parity'],record['source_outcome_parity'])
        for result in record.get('diagnostics',[]):
            last=result['last']
            print(' ',result['file'],'call',last['solver_call'],
                  'rho before/after:',last['rho_before_p'],last['rho_after_p'],
                  'internal rho:',last.get('rho_internal','--'),'relative P error:',last['p_error_rel'])
            print('   ',result['finding'])
        for error in record['errors']:print('  INVALID DIAGNOSTIC:',error)
        bad=bad or bool(record['errors'])
        (root/'summary.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
        if record['errors']:
            print('Stopping: resolve diagnostic validity before interpreting later runs.')
            break
    if digest(model)!=sha:raise RuntimeError('Executable changed during replay')
    if args.grouping_call and not bad:
        selected=[records.get(m,{}).get('grouping',[]) for m in ['full-fcg8','full-sd8']]
        if all(selected):
            eligible=[{int(r['group_size']) for r in rows if r['bits']==8 and r['screen_eligible']}
                      for rows in selected]
            print('Common 8-bit groups with positive rho AND fewer modeled P bytes than same-state P16:',
                  sorted(eligible[0]&eligible[1]))
        print('Snapshot screening only. No candidate was used to advance the solver; convergence remains untested.')
    print('Expected solver failures are retained. Diagnostic parity is not solver convergence.')
    print('Extra shadow solves/reductions are NOT compression traffic; do not compare these timings.')
    print('Results:',root)
    return int(bad)


if __name__=='__main__':raise SystemExit(main())
