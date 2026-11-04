#!/usr/bin/env python3
"""Read-only first-packet audit; replay saved failing Pepe CG16 with observer off/on."""
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

STAGES=('wire','metadata','symbols','coarse','fine','offsets','neighbors','weights','nodal','values')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def read_vector(path, stage):
    lines=path.read_text().splitlines()
    if not lines:raise ValueError('Empty audit file: '+str(path))
    convert=float if stage in ('metadata','values','weights') else int
    values=[convert(x) for x in lines[1:]]
    if len(values)!=int(lines[0]):raise ValueError('Truncated audit file: '+str(path))
    if stage in ('metadata','values','weights') and not all(math.isfinite(x) for x in values):
        raise ValueError('Nonfinite audit values: '+str(path))
    return values

def inspect_packets(folder):
    result={'phases':{},'mismatches':[],'incomplete':[],'scope':'First packet on each directed interface; encoded nodal payload, with full packet byte comparison including trailers.'}
    send_files=list(folder.glob('*_send_*_wire.txt'))
    if not send_files:raise ValueError('No audit packets; verify rebuilt executable and MPI environment forwarding')
    for path in sorted(send_files):
        match=re.fullmatch(r'(restrict|prolongate)_send_(\d+)_(\d+)_wire.txt',path.name)
        if not match:raise ValueError('Unknown audit filename: '+path.name)
        phase,producer,consumer=match.groups()
        counts=result['phases'].setdefault(phase,{'packets':0,'checked':{s:0 for s in STAGES},'mismatched':{s:0 for s in STAGES},'sender_inverse_checked':0,'sender_inverse_mismatched':0})
        counts['packets']+=1
        for stage in STAGES:
            sender=folder/f'{phase}_send_{producer}_{consumer}_{stage}.txt'
            receiver=folder/f'{phase}_recv_{consumer}_{producer}_{stage}.txt'
            if not sender.exists() or not receiver.exists():
                result['incomplete'].append({'phase':phase,'producer':producer,'consumer':consumer,'stage':stage})
                continue
            a=read_vector(sender,stage);b=read_vector(receiver,stage)
            counts['checked'][stage]+=1
            if a!=b:
                counts['mismatched'][stage]+=1
                first=next((i for i,(x,y) in enumerate(zip(a,b)) if x!=y),min(len(a),len(b)))
                result['mismatches'].append({'phase':phase,'producer':producer,'consumer':consumer,'stage':stage,
                    'sender_count':len(a),'receiver_count':len(b),'first_index':first,
                    'sender':a[first] if first<len(a) else None,'receiver':b[first] if first<len(b) else None})
        own_inverse=folder/f'{phase}_send_{producer}_{consumer}_selfnodal.txt'
        own_nodal=folder/f'{phase}_send_{producer}_{consumer}_nodal.txt'
        if own_inverse.exists() and own_nodal.exists():
            counts['sender_inverse_checked']+=1
            a=read_vector(own_nodal,'nodal');b=read_vector(own_inverse,'nodal')
            if a!=b:
                counts['sender_inverse_mismatched']+=1
                first=next((i for i,(x,y) in enumerate(zip(a,b)) if x!=y),min(len(a),len(b)))
                result['mismatches'].append({'phase':phase,'producer':producer,'consumer':consumer,
                    'stage':'sender_local_inverse','first_index':first,
                    'sender':a[first] if first<len(a) else None,'receiver':b[first] if first<len(b) else None})
        else:result['incomplete'].append({'phase':phase,'producer':producer,'consumer':consumer,'stage':'sender_local_inverse'})
    for path in folder.glob('*_recv_*_wire.txt'):
        phase,_,consumer,producer,_=path.stem.split('_')
        if not (folder/f'{phase}_send_{producer}_{consumer}_wire.txt').exists():
            result['incomplete'].append({'phase':phase,'producer':producer,'consumer':consumer,'stage':'missing_sender_wire'})
    result['status']='MISMATCH' if result['mismatches'] else ('INCOMPLETE' if result['incomplete'] or set(result['phases'])!={'restrict','prolongate'} else 'PASS')
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run',type=Path,default=Path('output/emi-cg-split-check-pepe-20261104T155620Z-3959/cg16'))
    p.add_argument('--reference-run',type=Path,help='Validate a repair against a successful matched quantization-only run')
    p.add_argument('--kaskade',type=Path,default=Path('/kaskade7'))
    args=p.parse_args();repo=Path.cwd();source=args.source_run.resolve()
    saved=shlex.split((source/'command.sh').read_text())
    def value(flag):
        if saved.count(flag)!=1:raise ValueError('Expected one '+flag)
        return saved[saved.index(flag)+1]
    required={'--bddcUseCg':1,'--bddcFlexibleCg':0,'--bddcFgmres':0,'--bddcCgOriginalResidual':1,
              '--maximumNumberOfTimeSteps':1,'--bddcRestrictionPasses':1,'--refine':0,'-np':8,
              '--bddcCompression':1,'--bddcCompressionBits':16,'--bddcRestrictCompressionBits':16,
              '--bddcProlongateCompressionBits':16,'--bddcGraphLifting':1,'--bddcGraphLiftingAfterQuantization':1,
              '--bddcZigzag':1,'--bddcBitlength':1,'--bddcHuffman':1,'--bddcProlongateConstantSplit':0,
              '--bddcRestrictGroupSize':0,'--bddcProlongateGroupSize':0,'--bddcRestrictZeroPreserving':0,
              '--bddcProlongateZeroPreserving':0,'--bddcRelativeTolerance':1e-7,'--bddcIterations':300,
              '--bddcHistogram':0,'--bddcCgDiagnostics':0,'--bddcProlongationDiagnostics':0}
    for flag,wanted in required.items():
        if float(value(flag))!=wanted:raise ValueError('Unexpected source option: '+flag)
    if Path(value('--input')).name!='pepe_combi_domi.vtu':raise ValueError('Expected Pepe source')
    reference=args.reference_run.resolve() if args.reference_run else source
    if args.reference_run:
        reference_command=shlex.split((reference/'command.sh').read_text())
        expected=list(saved)
        for flag in ('--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                     '--bddcZigzag','--bddcBitlength','--bddcHuffman'):
            expected[expected.index(flag)+1]='0'
        expected[expected.index('--dir')+1]=reference_command[reference_command.index('--dir')+1]
        if reference_command!=expected:
            raise ValueError('Reference must differ only in the five reversible coding flags and output directory')
    model=repo/'emiModel'
    if saved.count(str(model))!=1:raise ValueError('Source executable differs from this repository')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=repo/'output'/f'emi-transfer-audit-pepe-{stamp}-{os.getpid()}';root.mkdir(exist_ok=False)
    snapshots=[repo/'EMI_model.cpp',repo/'EMI_bddc.hh',args.kaskade/'mg/bddcSpaceTransfer.hpp',args.kaskade/'mg/bddcTransferAudit.hh',args.kaskade/'mg/bddcLiftingPlan.hh',args.kaskade/'mg/bddcSpaceTransfer.hh',args.kaskade/'mg/bddc.hpp',Path(__file__)]
    manifest={'source_run':str(source),'source_command':saved,'binary_sha256':digest(model),
              'source_sha256':{str(x):digest(x) for x in snapshots},'reference_run':str(reference),'reference_trace_sha256':digest(reference/'bddc-convergence-step1.csv'),
              'purpose':'Producer-plan repair check' if args.reference_run else 'Read-only first-packet audit'}
    (root/'environment.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Transfer audit:',root,flush=True)
    records={}
    for mode in ('off','on'):
        if digest(model)!=manifest['binary_sha256']:raise RuntimeError('Executable changed')
        out=root/mode;out.mkdir();cmd=list(saved);cmd[cmd.index('--dir')+1]=str(out)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        env.pop('KASKADE_RESTRICTION_CAPTURE',None);env.pop('KASKADE_TRANSFER_AUDIT_DIR',None)
        if mode=='on':
            env['KASKADE_TRANSFER_AUDIT_DIR']=str(out/'packets')
            i=cmd.index(str(model));cmd[i:i]=['-x','KASKADE_TRANSFER_AUDIT_DIR']
        exports=['OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1']
        if mode=='on':exports+=['KASKADE_TRANSFER_AUDIT_DIR='+env['KASKADE_TRANSFER_AUDIT_DIR']]
        (out/'command.sh').write_text(shlex.join(['env','-u','KASKADE_RESTRICTION_CAPTURE','-u','KASKADE_TRANSFER_AUDIT_DIR']+exports+cmd)+'\n')
        print('Running saved CG16, audit '+mode+' ...',flush=True)
        with (out/'run.log').open('w') as f:rc=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
        records[mode]={'exit_code':rc}
        print('  application exit:',rc,flush=True)
    name='bddc-convergence-step1.csv'
    def trace(folder):
        with (folder/name).open() as f:rows=list(csv.DictReader(f))
        if not rows:raise ValueError('Missing convergence trace')
        return rows
    a=trace(root/'off');b=trace(root/'on');original=trace(source)
    parity=a==b and records['off']['exit_code']==records['on']['exit_code']
    source_parity=a==original
    reference_trace=trace(reference)
    reference_parity=a==reference_trace
    reference_converged=(math.isfinite(float(reference_trace[-1]['relative_residual'])) and
                         float(reference_trace[-1]['relative_residual'])<=1e-7)
    comparison_ok=reference_parity and (reference_converged if args.reference_run else True)
    summary={'runs':records,'observer_parity':parity,'saved_source_trace_parity':source_parity,
             'reference_trace_parity':reference_parity,'reference_converged':reference_converged}
    print('Observer parity:',parity,flush=True)
    if args.reference_run:
        print('Old failing trace reproduced:',source_parity,'(False expected after repair)',flush=True)
        print('Iterations:',a[-1]['iteration'],'; true relative residual:',a[-1]['relative_residual'],flush=True)
    else:
        print('Saved source trace parity:',source_parity,flush=True)
    print('Reference trace parity:',reference_parity,flush=True)
    if parity:
        summary['audit']=inspect_packets(root/'on/packets')
        audit=summary['audit']
        for phase,counts in audit['phases'].items():
            print(phase,': packets=',counts['packets'],'; stage mismatches=',counts['mismatched'],flush=True)
            print('  Sender inverse:',counts['sender_inverse_mismatched'],'mismatches /',counts['sender_inverse_checked'],'checked',flush=True)
        print('Packet audit:',audit['status'],'; incomplete stage records:',len(audit['incomplete']),flush=True)
        for item in audit['mismatches'][:6]:print('  First differences:',item,flush=True)
    else:print('STOP: observer parity failed; do not interpret packet comparisons.',flush=True)
    if args.reference_run:
        summary['repair_passed']=(parity and comparison_ok and
            all(r['exit_code']==0 for r in records.values()) and
            summary.get('audit',{}).get('status')=='PASS' and
            all(summary['audit']['phases'][phase]['packets']==204
                for phase in ('restrict','prolongate')))
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Results:',root,flush=True)
    print('Diagnostic I/O is not compression traffic or benchmark timing.')
    if args.reference_run:
        passed=summary['repair_passed']
        print('Producer-plan repair validation:', 'PASS' if passed else 'FAIL',flush=True)
        return int(not passed)
    # A detected mismatch is the diagnostic finding; failed observation is an error.
    return int(not parity or not comparison_ok or summary.get('audit',{}).get('status')=='INCOMPLETE')

if __name__=='__main__':raise SystemExit(main())
