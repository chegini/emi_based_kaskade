#!/usr/bin/env python3
"""Twelve controlled ownership configurations; never submits a SLURM job.

--plan-only is safe outside an allocation. Execution requires an existing
8-node/16-rank allocation and launches only the explicitly selected pilot.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import statistics

import run_emi_hb_time_study as hb
from collect_emi_scaling_smoke import collect
from analyze_emi_potential_diagnostics import read, error_metrics


def cases():
    result=[]
    for policy,mode in [('P2','baseline'),('P1','coupled'),('P3','separated')]:
        for variant in ['raw','A','B','B-recenter']:
            result.append(dict(name=policy+'-'+variant,variant=variant,policy=policy,
                partition_mode=mode,bits=None if variant=='raw' else 16,depth=1,
                split='off' if variant=='raw' else 'on',
                optimization=5 if variant=='A' else 4 if variant!='raw' else 0))
    return result


def readcsv(path):
    with path.open() as f:return list(csv.DictReader(f))


def writecsv(path, rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def serialized_entropy(folder, transport):
    """Read pooled byte entropy, verifying counts against actual routed bytes."""
    rows=readcsv(folder/'serialized-message-byte-entropy.csv')
    expected={(p,s) for p in ['restriction','prolongation'] for s in ['same_node','inter_node','inter_rank']}
    if len(rows)!=6 or {(r['phase'],r['scope']) for r in rows}!=expected:
        raise ValueError('Incomplete/duplicate serialized byte entropy summaries')
    values={}
    for row in rows:
        phase,scope=row['phase'],row['scope']
        samples=[r for r in transport if r['phase']==phase and r['scope']==scope]
        size=sum(int(r['sent_bytes']) for r in samples)
        messages=sum(int(r['mpi_messages']) for r in samples)
        if (int(row['total_message_bytes']),int(row['message_count']))!=(size,messages):
            raise ValueError('Serialized entropy population differs from actual MPI bytes/messages')
        if row['entropy_aggregation']!='pooled_all_messages_all_ranks' or row['run_status']!='PASS':
            raise ValueError('Invalid byte entropy aggregation/run status')
        value=row['serialized_byte_entropy_bits_per_byte']
        if size:
            value=float(value)
            if not math.isfinite(value) or not 0<=value<=8+1e-12:
                raise ValueError('Invalid empirical byte entropy')
        elif value!='':raise ValueError('Empty stream must have undefined entropy')
        else:value=None
        values[phase+'_'+scope+'_byte_entropy_bits_per_byte']=value
    return values


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def arguments(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan-only',action='store_true')
    p.add_argument('--byte-entropy',action='store_true',help='Opt-in pooled final-buffer entropy, raw and compressed')
    p.add_argument('--case',choices=['pepe'],default='pepe')
    p.add_argument('--nodes',type=int,default=8);p.add_argument('--ranks',type=int,default=16)
    p.add_argument('--steps',type=int,default=3)
    p.add_argument('--executable',type=Path,default=Path('emiModel'))
    p.add_argument('--kaskade',type=Path,default=Path('/data/numerik/people/fchegini/kaskade7_new'))
    p.add_argument('--slurm-container',type=Path,default=Path('/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif'))
    p.add_argument('--output',type=Path,default=Path('partition-pilot-results'))
    p.add_argument('--repeats',type=int,default=3);p.add_argument('--detail-repeats',type=int,default=3)
    p.add_argument('--warmups',type=int,default=1)
    p.add_argument('--case-minutes',type=int,default=30)
    p.add_argument('--intra-nodes',default='');p.add_argument('--extra-nodes',default='')
    p.add_argument('--balance-factor',type=float,default=1)
    p.add_argument('--seed',type=int,default=20261115)
    a=p.parse_args(argv)
    if (a.nodes,a.ranks,a.steps)!=(8,16,3):p.error('This pilot is fixed to Pepe, 8 nodes, 16 ranks, 3 steps')
    if min(a.repeats,a.detail_repeats)<3 or a.warmups<1:p.error('Require >=3 repeats and >=1 warmup')
    if not math.isfinite(a.balance_factor) or a.balance_factor<1:p.error('Balance factor must be >=1')
    return a


def command(a, case, folder, validation):
    cfg,cmd=hb.command(a,case,folder,validation=validation,reuse=True)
    # The fixed-topology HB application guard forbids the heavyweight CG
    # observer. Keep the supported true-residual history and state snapshots.
    cmd[cmd.index('--bddcCgDiagnostics')+1]='0'
    if getattr(a,'byte_entropy',False):cmd+=['--bddcSerializedByteEntropy','1']
    cmd += ['--bddcPartitionMode',case['partition_mode'],'--bddcPartitionReport','1',
            '--bddcPartitionBalanceFactor',str(a.balance_factor)]
    if case['policy']=='P3':
        if a.intra_nodes:cmd += ['--bddcIntracellularNodes',a.intra_nodes]
        if a.extra_nodes:cmd += ['--bddcExtracellularNodes',a.extra_nodes]
    if cmd[0]=='srun':cmd.insert(1,'--distribution=block:block')
    return cfg,cmd


def observations(folder, ranks, steps):
    transport=[];timing=[];completion=[]
    for rank in range(ranks):
        for stem,target in [('transport',transport),('timing',timing),('steps',completion)]:
            target.extend(dict(rank=rank,**r) for r in readcsv(folder/f'partition-{stem}-rank-{rank}.csv'))
    assert len(transport)==ranks*steps*2*3
    assert len(timing)==ranks*steps*2 and len(completion)==ranks*steps
    for rank in range(ranks):
        for step in range(1,steps+1):
            for phase in ['restriction','prolongation']:
                rr={r['scope']:r for r in transport if r['rank']==rank and int(r['step'])==step and r['phase']==phase}
                assert set(rr)=={'inter_rank','inter_node','same_node'}
                for field in ['raw_reference_bytes','prepacking_bytes','sent_bytes','mpi_messages']:
                    assert int(rr['inter_rank'][field])==int(rr['inter_node'][field])+int(rr['same_node'][field])
    return transport,timing,completion


def validate_topology(folder, transport):
    """Reconcile actual raw-reference bytes/message buckets with actual endpoints."""
    edges=readcsv(folder/'partition-interfaces.csv')
    predicted={}
    for e in edges:
        a,b=int(e['rank_a']),int(e['rank_b'])
        if a==b:continue
        scope='inter_node' if int(e['route'])==2 else 'same_node'
        for source,target in [(a,b),(b,a)]:
            key=(source,scope)
            entry=predicted.setdefault(key,dict(raw=0,peers=set()))
            entry['raw']+=16+8*int(e['interface_dofs']);entry['peers'].add(target)
    for r in transport:
        if r['scope']=='inter_rank':continue
        expected=predicted.get((r['rank'],r['scope']),dict(raw=0,peers=set()))
        applications=int(r['phase_applications'])
        assert int(r['raw_reference_bytes'])==expected['raw']*applications,'Actual R/P endpoints differ from setup graph'
        assert int(r['mpi_messages'])==len(expected['peers'])*applications,'Aggregated MPI message count differs from endpoints'


def main(argv=None):
    a=arguments(argv);configurations=cases()
    if a.plan_only:
        print(json.dumps(dict(status='PLAN_ONLY_NOT_EXECUTED',configurations=configurations,
            primary_trajectories=12*(a.repeats+a.warmups),
            diagnostic_trajectories=12*(a.detail_repeats+a.warmups),validation_trajectories=12,
            warning='No network-only elapsed time; no full depth/precision sweep.'),indent=2))
        return 0
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Pilot execution requires an explicitly approved existing SLURM allocation')
    if int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=8 or int(os.environ.get('SLURM_NTASKS','0'))!=16:
        raise RuntimeError('Pilot allocation must have exactly 8 nodes and 16 ranks')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    manifest=dict(vars(a),cases=configurations,executable_sha256=digest(a.executable),
        input_sha256={k:digest(Path(v)) for k,v in hb.mesh_inputs(Path(__file__).resolve().parents[1],a.case).items()},
        allocation={k:os.environ.get(k) for k in ['SLURM_JOB_ID','SLURM_JOB_NODELIST']},
        notes=['Only ownership varies across policies. State errors use each policy own raw reference.',
               'Inter-rank=inter-node+same-node. MPI times are mixed-scope collectives.',
               'Primary tier 0; kernel diagnostic tier 2. Parent/child timers are nested.',
               '28.47% is a prior graph-derived logical-payload bound, not a speedup.',
               'Map/CPU placement and matrix fingerprints must stay fixed across repeats.'])
    (root/'manifest.json').write_text(json.dumps(manifest,default=str,indent=2)+'\n')
    rows=[];validated={};state_errors=[];placement=None;problem=None;coarse=None
    all_steps=[];all_transport=[];all_timing=[];phase_runs=[]
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    for key in list(env):
        if key.startswith(('KASKADE_FROZEN','KASKADE_RESTRICTION_CAPTURE','KASKADE_TRANSFER_AUDIT','EMI_MESSAGE_ENTROPY_PREFIX','APPTAINERENV_KASKADE_')):
            env.pop(key)

    def save():
        (root/'runs.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
        writecsv(root/'validation-state-errors.csv',state_errors)
        writecsv(root/'steps.csv',all_steps);writecsv(root/'transport.csv',all_transport)
        writecsv(root/'timing.csv',all_timing)
        writecsv(root/'phase-runs.csv',phase_runs)

    def run(case,tier,repeat,validation):
        nonlocal placement,problem,coarse
        name=('validate' if validation else f't{tier}-r{repeat}')+'-'+case['name']
        folder=root/name;folder.mkdir()
        cfg,cmd=command(a,case,folder,validation)
        entry=dict(case,configuration=case['name'],tier=tier,repeat=repeat,validation=validation,
                   warmup=not validation and repeat<0,status='FAILED',folder=str(folder))
        rows.append(entry);print('Running '+name,flush=True)
        runenv=dict(env,KASKADE_LIFTING_OPTIMIZATION=str(case['optimization']),
            APPTAINERENV_KASKADE_LIFTING_OPTIMIZATION=str(case['optimization']),
            KASKADE_TIMESTEP_BYTES='1',APPTAINERENV_KASKADE_TIMESTEP_BYTES='1',
            KASKADE_TRANSFER_TIMING=str(tier),APPTAINERENV_KASKADE_TRANSFER_TIMING=str(tier),
            KASKADE_PACKET_TIMING=str(int(tier>0)),APPTAINERENV_KASKADE_PACKET_TIMING=str(int(tier>0)))
        (folder/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        try:
            rc,timeout,launcher=hb.study.execute(cmd,folder,cfg,runenv)
            log=(folder/'run.log').read_text(errors='replace')
            metrics=collect(log,a.steps,a.ranks,rc,1e-7)
            if timeout or metrics['errors']:raise RuntimeError('; '.join(metrics['errors'])+(' TIMEOUT' if timeout else ''))
            step_rows,stages,totals=hb.profiles(folder,a.ranks,a.steps)
            transport,timing,completion=observations(folder,a.ranks,a.steps)
            validate_topology(folder,transport)
            if a.byte_entropy:entry.update(serialized_entropy(folder,transport))
            current_placement=(folder/'partition-ranks.csv').read_bytes()
            placement_rows=readcsv(folder/'partition-ranks.csv')
            assert len(placement_rows)==a.ranks and len({r['hostname'] for r in placement_rows})==a.nodes,'Actual physical node count differs from requested pilot'
            node_hosts={}
            for r in placement_rows:node_hosts.setdefault(r['node_group'],set()).add(r['hostname'])
            assert len(node_hosts)==a.nodes and all(len(v)==1 for v in node_hosts.values()),'Shared-memory groups and physical hosts disagree'
            fingerprint=readcsv(folder/'partition-problem.csv')[0]
            invariant={k:v for k,v in fingerprint.items() if k not in ['mode','planning_s']}
            if placement is None:placement=current_placement
            if problem is None:problem=invariant
            assert placement==current_placement,'Rank-to-host or CPU binding changed within the pilot'
            assert problem==invariant,'Ownership changed the mesh/operator/overlap invariants'
            coarse_rows=(folder/'partition-coarse.csv').read_bytes()
            if coarse is None:coarse=coarse_rows
            assert coarse==coarse_rows,'Ownership changed the actual coarse constraint rows'
            if case['name'] in validated:
                previous=Path(validated[case['name']]['folder'])
                assert (previous/'partition-ownership.csv').read_bytes()==(folder/'partition-ownership.csv').read_bytes(),'Ownership changed between repeats'
            for step in range(1,a.steps+1):
                done=[r for r in completion if int(r['step'])==step]
                assert all(r['status']=='PASS' and int(r['accepted_updates'])==metrics['iterations_per_step'][step-1] for r in done)
                for phase in ['restriction','prolongation']:
                    assert all(int(r['phase_applications'])==metrics['iterations_per_step'][step-1] for r in transport if int(r['step'])==step and r['phase']==phase)
            actual=sum(int(r['sent_bytes']) for r in transport if r['scope']=='inter_rank')
            assert actual==sum(r['packed_bytes'] for r in step_rows),'New and existing byte counters disagree'
            if case['variant']=='raw':
                assert all(r['raw_reference_bytes']==r['sent_bytes'] for r in transport),'Raw reference serialization differs from actual raw messages'
            entry.update(metrics,status='PASS',launcher_s=launcher,
                simulation_s=max(t['simulation_s'] for t in totals),
                application_compute_s=max(t['application_compute_s'] for t in totals),
                fcg_call_s=max(sum(r['fcg_s'] for r in step_rows if r['rank']==rank) for rank in range(a.ranks)),
                retained_bddc_setup_s=max(r['wall_s'] for r in step_rows if not r['step']),
                application_setup_s=max(t['application_setup_s'] for t in totals),
                packed_bytes=actual,inter_node_bytes=sum(int(r['sent_bytes']) for r in transport if r['scope']=='inter_node'),
                inter_node_messages=sum(int(r['mpi_messages']) for r in transport if r['scope']=='inter_node'),
                mpi_payload_collective_s=max(sum(float(r['mpi_payload_collective_s']) for r in timing if r['rank']==rank) for rank in range(a.ranks)),
                process_peak_rss_bytes=max(t['process_peak_rss_bytes'] for t in totals),
                partition_ownership_sha256=digest(folder/'partition-ownership.csv'))
            work=readcsv(folder/'partition-work.csv')
            for scope in ['rank','node']:
                for field in ['local_dofs','local_matrix_nnz']:
                    values=[int(r[field]) for r in work if r['scope']==scope]
                    entry[scope+'_'+field+'_max_over_mean']=max(values)/statistics.mean(values)
            for kind,data in [('steps',step_rows),('transport',transport),('timing',timing)]:
                target={'steps':all_steps,'transport':all_transport,'timing':all_timing}[kind]
                target.extend(dict(run=name,configuration=case['name'],policy=case['policy'],tier=tier,repeat=repeat,validation=validation,warmup=entry['warmup'],**r) for r in data)
            for phase in ['restriction','prolongation']:
                tr=[r for r in transport if r['phase']==phase]
                tt=[r for r in timing if r['phase']==phase]
                applications=sum(int(r['phase_applications']) for r in tr if r['rank']==0 and r['scope']=='inter_rank')
                for scope in ['inter_rank','inter_node','same_node']:
                    samples=[r for r in tr if r['scope']==scope]
                    sent=sum(int(r['sent_bytes']) for r in samples)
                    count=sum(int(r['mpi_messages']) for r in samples)
                    pr=dict(run=name,configuration=case['name'],policy=case['policy'],tier=tier,
                        repeat=repeat,validation=validation,warmup=entry['warmup'],phase=phase,
                        scope=scope,phase_applications=applications,
                        accepted_updates=sum(metrics['iterations_per_step']),sent_bytes=sent,mpi_messages=count,
                        mean_sent_bytes=sent/count if count else None,
                        max_sent_bytes=max((int(r['max_sent_bytes']) for r in samples if r['max_sent_bytes']),default=None),
                        raw_reference_bytes=sum(int(r['raw_reference_bytes']) for r in samples),
                        bytes_per_application=sent/applications if applications else None,
                        bytes_per_accepted_update=sent/sum(metrics['iterations_per_step']))
                    phase_runs.append(pr)
                # Per-rank accumulated phase durations, then rank maximum.
                # Only inter-rank rows carry durations; network-only rows do not.
                for field in ['compression_s','decompression_s','phase_wall_s','serialization_s',
                              'deserialization_s','mpi_payload_collective_s','mpi_all_calls_s',
                              'float_forward_s','integer_forward_s','float_inverse_s','integer_inverse_s',
                              'quantize_s','dequantize_s','huffman_encode_s','huffman_decode_s','recenter_s','uncenter_s']:
                    value=max(sum(float(r[field]) for r in tt if r['rank']==rank) for rank in range(a.ranks)) if all(r[field] for r in tt) else None
                    phase_runs[-3][field+'_max_rank']=value
            if validation:
                reference=validated.get(case['policy']+'-raw')
                if reference:
                    for step in range(a.steps+1):
                        v=read(folder/f'potential-state-step{step}.csv')
                        raw=read(Path(reference['folder'])/f'potential-state-step{step}.csv')
                        if step==0:assert v['u']==raw['u'],'Initial states differ'
                        state_errors.append(dict(configuration=case['name'],policy=case['policy'],step=step,reference=reference['configuration'],**error_metrics(v['u'],raw['u'])))
                if case['variant']=='raw' and case['policy']!='P2':
                    reference=validated['P2-raw']
                    for step in range(a.steps+1):
                        v=read(folder/f'potential-state-step{step}.csv');raw=read(Path(reference['folder'])/f'potential-state-step{step}.csv')
                        if step==0:assert v['u']==raw['u'],'Raw placement changed the initial state'
                        state_errors.append(dict(configuration=case['name'],policy=case['policy'],step=step,reference='P2-raw',**error_metrics(v['u'],raw['u'])))
                if case['variant']=='B-recenter':
                    assert hb.exact_trajectory(Path(validated[case['policy']+'-B']['folder']),folder,a.steps),'Recentering changed reconstructed states/residual traces'
            print('  PASS',flush=True)
        except Exception as error:
            entry['status']='FAILED';entry['error']=str(error);print('  FAILED: '+str(error),flush=True)
            # Retain completed/partial failure records without normalizing by 3.
            for rank in range(a.ranks):
                path=folder/f'partition-steps-rank-{rank}.csv'
                if path.is_file():
                    entry.setdefault('failure_step_records',[]).extend(dict(rank=rank,**r) for r in readcsv(path))
        save();return entry

    # Gate on all three raw mappings before drawing any compression conclusions.
    for case in [c for c in configurations if c['variant']=='raw']:
        row=run(case,2,-1,True)
        if row['status']!='PASS':raise RuntimeError('Raw ownership correctness gate failed; inspect '+row['folder'])
        validated[case['name']]=row
    for case in [c for c in configurations if c['variant']!='raw']:
        if case['variant']=='B-recenter' and case['policy']+'-B' not in validated:
            rows.append(dict(case,configuration=case['name'],status='SKIPPED',reason='Corresponding B validation failed; B+R was not independently attempted',validation=True));save();continue
        row=run(case,2,-1,True)
        if row['status']=='PASS':validated[case['name']]=row
    active=[c for c in configurations if c['name'] in validated]
    rng=random.Random(a.seed)
    for tier,repeats in [(0,a.repeats),(2,a.detail_repeats)]:
        for repeat in range(-a.warmups,repeats):
            order=list(active);rng.shuffle(order)
            for case in order:run(case,tier,repeat,False)
    assert digest(a.executable)==manifest['executable_sha256'],'Executable changed during pilot'
    summarize(root,rows,configurations)
    phase_summary=[]
    for case in configurations:
        for tier in [0,2]:
            for phase in ['restriction','prolongation']:
                for scope in ['inter_rank','inter_node','same_node']:
                    samples=[r for r in phase_runs if r['configuration']==case['name'] and r['tier']==tier and r['phase']==phase and r['scope']==scope and not r['validation'] and not r['warmup']]
                    out=dict(configuration=case['name'],tier=tier,phase=phase,scope=scope,successful_repeats=len(samples))
                    fields=['phase_applications','accepted_updates','sent_bytes','mpi_messages','mean_sent_bytes','max_sent_bytes','raw_reference_bytes','bytes_per_application','bytes_per_accepted_update']
                    fields+=sorted({k for r in samples for k in r if k.endswith('_max_rank')})
                    for field in fields:
                        values=[r[field] for r in samples if r.get(field) is not None]
                        out[field+'_mean']=statistics.mean(values) if values else None
                        out[field+'_sample_sd']=statistics.stdev(values) if len(values)>1 else None
                    phase_summary.append(out)
    writecsv(root/'phase-comparison.csv',phase_summary)
    print('Results: '+str(root),flush=True)
    return int(any(r['status']!='PASS' for r in rows))


def summarize(root, runs, configurations):
    summaries=[];paired=[]
    for case in configurations:
        measured=[r for r in runs if r['configuration']==case['name'] and r.get('tier')==0 and not r.get('validation') and not r.get('warmup') and r['status']=='PASS']
        attempted=[r for r in runs if r['configuration']==case['name']]
        status='FAILED' if any(r['status']=='FAILED' for r in attempted) else 'SKIPPED' if any(r['status']=='SKIPPED' for r in attempted) else 'PASS' if len(measured)>=3 else 'INCOMPLETE'
        out=dict(configuration=case['name'],policy=case['policy'],method=case['variant'],status=status,successful_primary_repeats=len(measured))
        for field in ['simulation_s','application_compute_s','fcg_call_s','retained_bddc_setup_s','application_setup_s','packed_bytes','inter_node_bytes','inter_node_messages','mpi_payload_collective_s','process_peak_rss_bytes',
                      'rank_local_dofs_max_over_mean','node_local_dofs_max_over_mean',
                      'rank_local_matrix_nnz_max_over_mean','node_local_matrix_nnz_max_over_mean']:
            values=[r[field] for r in measured]
            out[field+'_mean']=statistics.mean(values) if values else None
            out[field+'_sample_sd']=statistics.stdev(values) if len(values)>1 else None
        out['iterations_per_step_observations']=json.dumps([r['iterations_per_step'] for r in measured])
        summaries.append(out)
        for r in measured:
            raw=next((x for x in runs if x['configuration']==case['policy']+'-raw' and x.get('tier')==0 and x.get('repeat')==r['repeat'] and not x.get('validation') and x['status']=='PASS'),None)
            if raw:
                n=sum(r['iterations_per_step']);nr=sum(raw['iterations_per_step'])
                paired.append(dict(configuration=case['name'],policy=case['policy'],repeat=r['repeat'],
                    network_byte_saving_percent=100*(1-r['inter_node_bytes']/raw['inter_node_bytes']),
                    application_speedup=raw['application_compute_s']/r['application_compute_s'],
                    fcg_call_speedup=raw['fcg_call_s']/r['fcg_call_s'],
                    mixed_payload_collective_time_reduction_percent=100*(1-r['mpi_payload_collective_s']/raw['mpi_payload_collective_s']) if raw['mpi_payload_collective_s'] else None,
                    inter_node_bytes_per_accepted_update=r['inter_node_bytes']/n,
                    additional_accepted_updates=n-nr,application_difference_s=r['application_compute_s']-raw['application_compute_s']))
    writecsv(root/'comparison.csv',summaries);writecsv(root/'paired-comparison.csv',paired)
    (root/'report.txt').write_text('Controlled EMI partition pilot; only existing approved allocation used.\n'
        'Scope: 8 nodes, 16 ranks, 3 timesteps, 16 bits, depth 1, split ON.\n'
        'MPI timing is mixed-scope collective elapsed time, not inter-node-only time.\n'
        'Parent/child kernel counters and independent rank maxima must not be added.\n'
        'Primary runtime tier 0; detailed codec timings in tier 2.\n'
        'Raw-reference bytes are counterfactual all-raw serialization of the same active endpoints.\n'
        'Use matched raw trajectories for full-run savings and speedup.\n'
        'The 28.47% bound is a prior graph-derived logical-byte estimate, not a measured speedup.\n')


if __name__=='__main__':raise SystemExit(main())
