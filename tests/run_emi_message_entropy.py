#!/usr/bin/env python3
"""Read-only message entropy pass, exact-checked against untimed validation runs."""
import argparse,csv,json,os,re,shlex
from pathlib import Path
from types import SimpleNamespace
import run_emi_hb_time_study as h
from summarize_emi_message_entropy import summarize


def write_json(path,data):
 temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,indent=2,allow_nan=False,default=str)+'\n');temp.replace(path)


def command(a,case,folder,original_kaskade):
 cfg,cmd=h.command(a,case,folder,validation=True,reuse=True)
 image=cmd.index(str(a.slurm_container))
 # Overlay symlinks point to the original absolute Kaskade path.
 cmd[image:image]=['--bind',f'{original_kaskade}:{original_kaskade}']
 return cfg,cmd


def environment(case,folder):
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
 for key in list(env):
  if key.startswith(('KASKADE_FROZEN','KASKADE_RESTRICTION_CAPTURE','KASKADE_TRANSFER_AUDIT','APPTAINERENV_KASKADE_')):env.pop(key)
 values=dict(KASKADE_LIFTING_OPTIMIZATION=str(case.get('optimization',0)),KASKADE_LIFTING_PROFILE_PREFIX='',
  KASKADE_TIMESTEP_BYTES='1',KASKADE_TRANSFER_TIMING='0',KASKADE_PACKET_TIMING='0',
  EMI_MESSAGE_ENTROPY_PREFIX=str(folder/'counts/entropy'),EMI_ENTROPY_MAX_BINS='500000')
 for k,v in values.items():env[k]=v;env['APPTAINERENV_'+k]=v
 return env


def aggregate(case,rows):
 out=[]
 for phase in ('restrict','prolongate'):
  for scope,routes in [('inter_rank',{'same_node','inter_node'}),('inter_node',{'inter_node'})]:
   q=[r for r in rows if r['phase']==phase and r['route'] in routes and r['stage']=='encoder_coefficients']
   wire=[r for r in rows if r['phase']==phase and r['route'] in routes and r['stage']=='routed_bytes']
   original=[r for r in rows if r['phase']==phase and r['route'] in routes and r['stage']=='quantized_input']
   lengths=[r for r in rows if r['phase']==phase and r['route'] in routes and r['stage']=='huffman_bitlength']
   n=sum(r['symbols'] for r in q);n_input=sum(r['symbols'] for r in original)
   assert n==n_input,(case['name'],phase,scope,n,n_input)
   weighted=lambda rr:sum(r['ideal_message_entropy_bits'] for r in rr)/sum(r['symbols'] for r in rr) if sum(r['symbols'] for r in rr) else None
   b=sum(r['symbols'] for r in wire)
   out.append(dict(configuration=case['name'],method=case['variant'],bits=case['bits'],depth=case['depth'],split=case.get('split','off'),phase=phase,scope=scope,
    quantized_coefficients=n_input,encoder_coefficients=n,coefficient_messages=sum(r['messages'] for r in q),
    quantized_input_H0_bits_per_coefficient=weighted(original),encoder_H0_bits_per_coefficient=weighted(q),
    huffman_length_H0_bits_per_length_symbol=weighted(lengths),
    routed_bytes=b,serialized_byte_H0_bits_per_byte=weighted(wire),
    effective_routed_bits_per_quantized_coefficient=8*b/n if n else None))
 return out


def save_csv(path,rows):
 if not rows:return
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main(argv=None):
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--performance-results',type=Path,required=True);p.add_argument('--executable',type=Path,required=True)
 p.add_argument('--overlay',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--plan-only',action='store_true')
 args=p.parse_args(argv)
 manifest=json.loads((args.performance_results/'manifest.json').read_text())
 a=SimpleNamespace(**manifest);original_kaskade=Path(a.kaskade)
 a.slurm_container=Path(a.slurm_container);a.output=args.output.resolve();a.kaskade=args.overlay.resolve();a.executable=args.executable.resolve()
 baseline_runs=json.loads((args.performance_results/'runs.json').read_text())
 refs={r['configuration']:r for r in baseline_runs if r.get('validation') and r.get('reuse') and r['status']=='PASS' and (r['configuration']=='raw' or r.get('exact_parity'))}
 if args.plan_only:
  for case in manifest['cases']:
   print(case['name'], 'DIAGNOSE' if case['name'] in refs else 'SKIP: baseline invalid')
   print(shlex.join(command(a,case,a.output/case['name'],original_kaskade)[1]))
  return 0
 h.require(int(os.environ.get('SLURM_JOB_NUM_NODES',0))>=a.nodes,'Insufficient allocation')
 a.output.mkdir(parents=True,exist_ok=False)
 executable_hash=h.digest(a.executable)
 # The input snapshot must be the same as the performance experiment.
 inputs=h.mesh_inputs(Path(h.__file__).resolve().parents[1],a.case)
 for key,path in inputs.items():h.require(h.digest(Path(path))==manifest['input_sha256'][key],'Mesh changed: '+key)
 write_json(a.output/'manifest.json',dict(performance_results=args.performance_results,performance_executable_sha256=manifest['executable_sha256'],
  observer_executable_sha256=executable_hash,overlay=json.loads((a.kaskade/'entropy-overlay.json').read_text()),
  input_sha256=manifest['input_sha256'],nodes=a.nodes,ranks=a.ranks,steps=a.steps,
  note='Entropy diagnostic only. Timing comparisons must use original performance results, with entropy disabled.'))
 records=[];all_rows=[];summary=[];failed=False;passed={}
 for case in manifest['cases']:
  name=case['name'];record=dict(configuration=name,status='SKIP',reason='No complete exact-validated baseline')
  records.append(record)
  if name not in refs:
   write_json(a.output/'runs.json',records);print('Entropy skip: '+name,flush=True);continue
  baseline=refs[name];folder=a.output/name;folder.mkdir()
  cfg,cmd=command(a,case,folder,original_kaskade)
  (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
  print('Entropy diagnostic: '+name,flush=True)
  try:
   rc,timeout,wall=h.study.execute(cmd,folder,cfg,environment(case,folder))
   log=(folder/'run.log').read_text(errors='replace');metrics=h.collect(log,a.steps,a.ranks,rc,1e-7)
   h.require(not timeout and not metrics['errors'],'; '.join(metrics['errors'])+(' timeout' if timeout else ''))
   rank_rows,_,_=h.profiles(folder,a.ranks,a.steps)
   packed=re.findall(r'^BDDC packet transport bytes: original=(\d+) actual=(\d+) packing=([01])$',log,re.M)
   h.require(len(packed)==1,'Missing actual packet bytes');actual=int(packed[0][1])
   h.require(actual==baseline['packed_bytes']==sum(r['packed_bytes'] for r in rank_rows),'Entropy observer changed byte count')
   h.require(not any(r['packed_bytes'] for r in rank_rows if r['step']==0),'Setup R/P messages require separate attribution')
   residuals=[r['relative_residual'] for r in rank_rows if r['rank']==0 and r['step']>0]
   h.require(metrics['iterations_per_step']==baseline['iterations_per_step'] and residuals==baseline['residual_history'],'Entropy observer changed convergence')
   h.require(h.exact_trajectory(Path(baseline['folder']),folder,a.steps),'Entropy observer changed state or full residual trace')
   rows,ranks=summarize(folder/'counts',folder/'summary')
   h.require(ranks==list(range(a.ranks)),'Missing rank entropy files')
   h.require(sum(r['symbols'] for r in rows if r['stage']=='routed_bytes')==actual,'Entropy byte count differs from actual routed bytes')
   if case['bits'] is not None:
    h.require(all(any(r['phase']==phase and r['stage']=='encoder_coefficients' and r['symbols']>0 for r in rows) for phase in ['restrict','prolongate']),'Missing actual coefficient population')
   h.require(not any(r['stage'] in ['quantized_input','encoder_coefficients'] and not r['encoded'] and r['symbols'] for r in rows),'Unencoded candidate coefficients present; do not interpret them as sent symbols')
   aggregated=aggregate(case,rows)
   if case['variant']=='B-recenter':
    other=passed.get(name.replace('B-recenter','B'));h.require(other is not None,'Missing B baseline for recentering')
    for current,old in zip(aggregated,other):
     for field in ['quantized_coefficients','coefficient_messages']:
      h.require(current[field]==old[field],'Recentering changed population')
     for field in ['quantized_input_H0_bits_per_coefficient','encoder_H0_bits_per_coefficient']:
      x,y=current[field],old[field]
      h.require(x is None and y is None or x is not None and y is not None and abs(x-y)<1e-8,'Recentering changed per-message H0')
   record.update(status='PASS',reason=None,exact_parity='EXACT_STATES_FULL_RESIDUAL_TRACES_ITERATIONS_AND_PACKET_BYTES',packed_bytes=actual,iterations_per_step=metrics['iterations_per_step'],baseline_folder=baseline['folder'])
   passed[name]=aggregated
   all_rows.extend(dict(configuration=name,**r) for r in rows);summary.extend(aggregated)
   print('  PASS: entropy + exact numerical/packet parity',flush=True)
  except Exception as e:
   failed=True;record.update(status='FAIL',reason=str(e));print('  FAIL: '+str(e),flush=True)
  write_json(a.output/'runs.json',records)
  save_csv(a.output/'entropy-comparison.csv',all_rows);save_csv(a.output/'message-entropy-summary.csv',summary)
 h.require(h.digest(a.executable)==executable_hash,'Observer executable changed')
 write_json(a.output/'status.json',dict(status='FAIL' if failed else 'PASS',passed=len(passed),skipped=sum(r['status']=='SKIP' for r in records),failed=sum(r['status']=='FAIL' for r in records),
  definitions=['H0=-sum p*log2(p); primary message-weighted H0=sum N_message*H0_message / sum N_message.',
  'Coefficient populations exclude exact corners, split constants, headers and synthetic training seeds.',
  'Routed-byte H0 includes actual R/P serialization. This is bits/byte, not bits/coefficient.',
  'Raw has no quantized coefficient entropy; its byte entropy is reported separately.',
  'effective routed bits/coefficient includes metadata/corners/constants; it is not a pure code length.',
  'Entropy instrumentation is excluded from all performance measurements.']))
 print('Message entropy results: '+str(a.output),flush=True)
 return int(failed)
if __name__=='__main__':raise SystemExit(main())
