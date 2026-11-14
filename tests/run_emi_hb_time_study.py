#!/usr/bin/env python3
"""Repeated fixed-topology EMI trajectories; mode 0 is the runtime endpoint."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import random
import re
import shlex

import run_emi_solver_study as study
from collect_emi_scaling_smoke import collect
from run_emi_mpi_routing_check import digest, routing_counts
from run_emi_histogram_start_check import mesh_inputs, require
from analyze_emi_potential_diagnostics import read, error_metrics


def arguments(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--kaskade',type=Path,required=True)
    p.add_argument('--slurm-container',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case',default='pepe')
    p.add_argument('--scaling',choices=['pilot','strong','weak'],default='pilot')
    p.add_argument('--nodes',type=int,default=2)
    p.add_argument('--ranks',type=int,default=4)
    p.add_argument('--steps',type=int,default=10)
    p.add_argument('--bits',type=int,nargs='+',choices=[8,12,16,24,32],default=[16])
    p.add_argument('--lifting-optimizations',type=int,nargs='+',choices=range(8),default=[0])
    p.add_argument('--hierarchy-profile',action='store_true',help='Separate diagnostic run; setup CSV I/O perturbs timings')
    p.add_argument('--depth-split-study',action='store_true',help='A split OFF/ON and all four B split/recenter combinations')
    p.add_argument('--interface-inventory',action='store_true',help='Setup graph inventory during validation only')
    p.add_argument('--depths',type=int,nargs='+',choices=range(65),default=[1,4])
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--detail-repeats',type=int,default=3)
    p.add_argument('--tiers',type=int,nargs='+',choices=[0,1,2],default=[0,1,2])
    p.add_argument('--warmups',type=int,default=1)
    p.add_argument('--case-minutes',type=int,default=60)
    p.add_argument('--seed',type=int,default=6111327)
    p.add_argument('--plan-only',action='store_true')
    a=p.parse_args(argv)
    if a.nodes<2 or a.ranks<a.nodes or a.ranks%a.nodes:p.error('Use >=2 nodes with equal ranks per node')
    if min(a.repeats,a.detail_repeats)<3 or a.warmups<1:p.error('Require >=3 repeats and >=1 warmup')
    if a.steps<2:p.error('Use >=2 timesteps to distinguish startup from steady work')
    if 0 not in a.tiers or len(set(a.tiers))!=len(a.tiers):p.error('Unique tiers including primary mode 0 required')
    if len(set(a.bits))!=len(a.bits) or len(set(a.depths))!=len(a.depths):p.error('Duplicate configuration')
    if 0 not in a.lifting_optimizations or len(set(a.lifting_optimizations))!=len(a.lifting_optimizations):p.error('Unique optimizations including original 0 required')
    if a.lifting_optimizations[0]!=0:p.error('Original optimization 0 must come first for exact validation')
    if 1 not in a.depths:p.error('Depth 1 baseline is required')
    if a.depth_split_study and a.depths[0]!=1:p.error('Depth/split study requires depth 1 first for method-A exact validation')
    if a.depth_split_study and a.lifting_optimizations!=[0]:p.error('Depth/split mode selects validated masks A=5 and B=4 itself')
    if a.scaling=='weak' and a.case!=f'weak_{2*a.ranks}_h16':
        p.error('Weak family requires weak_<2*ranks>_h16; two cells/rank, fixed h16')
    return a


def configurations(a):
    if getattr(a,'depth_split_study',False):
        # A already has reversible integer lifting/zigzag; recentering is a B representation.
        return [dict(name='raw',variant='raw',bits=None,depth=1,optimization=0,split='off')]+[
            dict(name=f'{v}-b{b}-d{d}-split-{split}',baseline_name=None,
                 variant=v,bits=b,depth=d,split=split,optimization=5 if v=='A' else 4)
            for b in a.bits for d in a.depths for split in ['off','on']
            for v in ['A','B','B-recenter']]
    return [dict(name='raw',variant='raw',bits=None,depth=1,optimization=0)]+[
        dict(name=f'{v}-b{b}-d{d}'+(f'-o{o}' if o else ''),
             baseline_name=f'{v}-b{b}-d{d}',variant=v,bits=b,depth=d,optimization=o)
        for b in a.bits for d in a.depths for o in getattr(a,'lifting_optimizations',[0])
        for v in ['A','B','B-recenter']]


def command(a,case,folder,*,validation=False,reuse=True):
    repo=Path(__file__).resolve().parents[1]
    cfg=study.arguments(['--repo',str(repo),'--case',a.case,'--solvers','fcg',
        '--split',case.get('split','off'),'--bits',str(case['bits'] or 16),
        '--ranks',str(a.ranks),'--expected-nodes',str(a.nodes),
        '--slurm-container',str(a.slurm_container),'--kaskade',str(a.kaskade),
        '--executable',str(a.executable),'--case-minutes',str(a.case_minutes)])
    cfg.steps=a.steps;cfg.lifting_levels=case['depth'];cfg.lossless_packing=1
    cfg.transform_order='after' if case['variant']=='A' else 'before'
    mode=dict(case,solver='fcg',direct=False);mode['split']=case.get('split','off')
    result=study.make_command(cfg,mode,folder,['mpirun'])
    for flag,value in [('potentialDiagnostics',int(validation)),('bddcConvergenceReport',int(validation))]:
        result[result.index('--'+flag)+1]=str(value)
    result+=['--bddcTimeStudy','1','--bddcReuseLiftingPlans',str(int(reuse)),
             '--bddcTransformFirstRecentering',str(int(case['variant']=='B-recenter'))]
    # The build and results share the job directory, which can live on scratch.
    # Bind it explicitly rather than relying on image/site default mounts.
    image=result.index(str(cfg.slurm_container))
    job_root=a.output.resolve().parent
    result[image:image]=['--bind',f'{job_root}:{job_root}']
    return cfg,result


def profiles(folder,ranks,steps):
    rank_rows=[];stage_rows=[];totals=[]
    for rank in range(ranks):
        with (folder/f'hb-time-rank-{rank}.csv').open() as f:rows=list(csv.DictReader(f))
        require([int(r['step']) for r in rows]==list(range(steps+1)),'Incomplete timestep records')
        for row in rows:
            row={k:(int(v) if k in ['step','calls'] or k.endswith('_bytes') else float(v)) for k,v in row.items()}
            require(all(math.isfinite(v) and v>=0 for v in row.values()),'Invalid timestep sample')
            rank_rows.append(dict(rank=rank,**row))
        with (folder/f'hb-time-rank-{rank}-stages.csv').open() as f:
            for row in csv.DictReader(f):
                stage_rows.append(dict(rank=rank,**{k:(float(v) if k=='seconds' else int(v) if k in ['step','call','level','calls'] else v) for k,v in row.items()}))
        with (folder/f'hb-time-rank-{rank}-total.csv').open() as f:
            values=list(csv.DictReader(f))
        require(len(values)==1,'Missing simulation timer')
        totals.append({k:float(v) for k,v in values[0].items()})
    require(all(math.isfinite(t) and t>=0 for row in totals for t in row.values()),'Invalid simulation time')
    return rank_rows,stage_rows,totals


def exact_trajectory(first,second,steps):
    files=[f'potential-state-step{s}.csv' for s in range(steps+1)]
    files += [f'bddc-convergence-step{s}.csv' for s in range(1,steps+1)]
    return all(digest(first/f)==digest(second/f) for f in files)


def main(argv=None):
    a=arguments(argv);cases=configurations(a)
    if a.plan_only:
        print(json.dumps(dict(configuration=vars(a),cases=cases,
            notes=['Generic graph hierarchy; no path-only assumption on EMI interfaces.',
                   'No adaptive recentering: A, B/off and B/on remain explicit.',
                   'Warmups are entire trajectories. Setup is retained within each trajectory.']),default=str,indent=2))
        for case in cases:print(shlex.join(command(a,case,a.output/case['name'])[1]))
        return 0
    require(int(os.environ.get('SLURM_JOB_NUM_NODES','0'))>=a.nodes,'Insufficient Slurm allocation')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    repo=Path(__file__).resolve().parents[1]
    executable_hash=digest(a.executable)
    (root/'manifest.json').write_text(json.dumps(dict(vars(a),cases=cases,
        executable_sha256=executable_hash,
        input_sha256={k:digest(Path(v)) for k,v in mesh_inputs(repo,a.case).items()},
        allocation={k:os.environ.get(k) for k in ['SLURM_JOB_ID','SLURM_JOB_NODELIST']},
        notes=['Generic hierarchy, no coarse update. Explicit split field; original study defaults OFF.',
               'Actual routed bytes include same-node inter-rank MPI; not all network wire bytes.',
               'Live A/B/depth trajectories can differ; compare state error and actual FCG counts.',
               'Validation snapshots/traces enabled only for validation runs.',
               'Primary runtime uses mode 0; detailed kernel observations have timer overhead.',
               'MPI waits overlap other ranks work; never add separate rank maxima.']),default=str,indent=2)+'\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    for k in list(env):
        if k.startswith(('KASKADE_FROZEN','KASKADE_RESTRICTION_CAPTURE','KASKADE_TRANSFER_AUDIT','APPTAINERENV_KASKADE_')):env.pop(k)
    records=[];steps_out=[];stages_out=[];errors=[];validated={};failed=set()

    def save():
        for name,rows in [('runs',records),('errors',errors)]:
            target=root/(name+'.json');temporary=root/(name+'.json.tmp')
            temporary.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n');temporary.replace(target)

    def run(case,tier,repeat,validation=False,reuse=True,suffix=''):
        name=('validate' if validation else f't{tier}-r{repeat}')+'-'+case['name']+suffix
        folder=root/name;folder.mkdir()
        cfg,cmd=command(a,case,folder,validation=validation,reuse=reuse)
        optimization=str(case.get('optimization',0))
        prefix=str(folder/'lifting-setup') if a.hierarchy_profile or (a.interface_inventory and validation and reuse) else ''
        runenv=dict(env,KASKADE_LIFTING_OPTIMIZATION=optimization,
                   APPTAINERENV_KASKADE_LIFTING_OPTIMIZATION=optimization,
                   KASKADE_LIFTING_PROFILE_PREFIX=prefix,
                   APPTAINERENV_KASKADE_LIFTING_PROFILE_PREFIX=prefix,KASKADE_TIMESTEP_BYTES='1',APPTAINERENV_KASKADE_TIMESTEP_BYTES='1',KASKADE_TRANSFER_TIMING=str(tier),KASKADE_PACKET_TIMING=str(int(tier>0)),
                   APPTAINERENV_KASKADE_TRANSFER_TIMING=str(tier),APPTAINERENV_KASKADE_PACKET_TIMING=str(int(tier>0)))
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        row=dict(case,configuration=case['name'],tier=tier,repeat=repeat,
                 validation=validation,warmup=not validation and repeat<0,reuse=reuse,folder=str(folder),status='FAIL')
        records.append(row);print('Running '+name,flush=True)
        try:
            rc,timeout,wall=study.execute(cmd,folder,cfg,runenv)
            log=(folder/'run.log').read_text(errors='replace')
            metrics=collect(log,a.steps,a.ranks,rc,1e-7)
            require(not timeout and not metrics['errors'],'; '.join(metrics['errors'])+(' timeout' if timeout else ''))
            rows,stages,total=profiles(folder,a.ranks,a.steps)
            row.update(metrics,launcher_s=wall,simulation_s=max(t['simulation_s'] for t in total),
                       application_compute_s=max(t['application_compute_s'] for t in total),
                       process_peak_rss_bytes=max((t['process_peak_rss_bytes'] for t in total if 'process_peak_rss_bytes' in t),default=None),
                       application_setup_s=max(t['application_setup_s'] for t in total),validation_only=validation,status='PASS')
            work=[t.get('owned_local_dofs',0) for t in total]
            row['owned_dofs_rank_mean']=sum(work)/len(work)
            row['owned_dofs_rank_min']=min(work);row['owned_dofs_rank_max']=max(work)
            row['owned_subdomains_total']=sum(t.get('owned_subdomains',0) for t in total)
            row['routes']=routing_counts(log,a.nodes)
            packet=re.findall(r'^BDDC packet transport bytes: original=(\d+) actual=(\d+) packing=([01])$',log,re.M)
            require(len(packet)==1,'Missing transport bytes')
            row['packed_bytes']=int(packet[0][1])
            require(sum(r['packed_bytes'] for r in rows)==row['packed_bytes'],'Per-step transport accounting differs')
            if a.depth_split_study:
                require(all(r['restrict_packed_bytes']+r['prolongate_packed_bytes']==r['packed_bytes'] for r in rows),
                        'Per-phase transport accounting differs')
            require(int(packet[0][2])==int(case['bits'] is not None),'Wrong packing mode')
            if reuse and tier>0 and case['bits'] is not None:
                require(not any(r['stage']=='plan_setup' and r['step']>0 for r in stages),'Hierarchy rebuilt during timesteps')
                require(any(r['stage']=='plan_setup' and r['step']==0 for r in stages),'Hierarchy setup was not observed')
            for r in rows:steps_out.append(dict(run=name,configuration=case['name'],tier=tier,repeat=repeat,validation=validation,warmup=row['warmup'],**r))
            for r in stages:stages_out.append(dict(run=name,configuration=case['name'],tier=tier,repeat=repeat,validation=validation,warmup=row['warmup'],**r))
            row['residual_history']=[r['relative_residual'] for r in rows if r['rank']==0 and r['step']>0]
            if not validation:
                ref=validated[case['name']]
                require(row['packed_bytes']==ref['packed_bytes'],'Repeat/timer changed packed byte count')
                require(row['iterations_per_step']==ref['iterations_per_step'] and row['residual_history']==ref['residual_history'],
                        'Profiling/repeat changed convergence trajectory')
            row['hierarchy_setup_payload']=[dict(rank=int(rank),local=int(local),exported=int(remote)) for rank,local,remote in
                re.findall(r'^BDDC integer lifting setup rank=(\d+) local_plan_bytes=(\d+) exported_remote_bytes=(\d+)',log,re.M)]
            print('  PASS',flush=True)
        except Exception as error:
            row['status']='FAIL';row['error']=str(error);print('  FAIL: '+str(error),flush=True)
        save();return row

    raw=run(cases[0],2,-1,validation=True);require(raw['status']=='PASS','Raw trajectory failed; inspect '+raw['folder'])
    validated['raw']=raw
    for case in cases[1:]:
        if case['variant']=='B-recenter' and case['name'].replace('B-recenter','B') in failed:
            failed.add(case['name']);continue
        if not a.depth_split_study and case.get('optimization',0) and case['baseline_name'] not in validated:
            failed.add(case['name']);continue
        original=None
        if case['variant']!='B-recenter' and (a.depth_split_study or not case.get('optimization',0)):original=run(case,0,-1,validation=True,reuse=False,suffix='-uncached')
        current=run(case,2,-1,validation=True)
        if current['status']!='PASS' or original and original['status']!='PASS':
            failed.add(case['name']);continue
        baseline=validated[case['baseline_name']] if case.get('optimization',0) and not a.depth_split_study else original or validated[case['name'].replace('B-recenter','B')]
        if case.get('optimization',0) and not a.depth_split_study:
            require(current['packed_bytes']==baseline['packed_bytes'],'Plan optimization changed packed byte count')
            require(current['iterations_per_step']==baseline['iterations_per_step'],'Plan optimization changed FCG iterations')
        require(exact_trajectory(Path(baseline['folder']),Path(current['folder']),a.steps),
                'CACHE/RECENTERING CHANGED SOLUTION OR RESIDUAL TRACE: '+case['name'])
        if a.depth_split_study and case['variant']=='A' and case['depth']!=1:
            first=validated[case['name'].replace('-d'+str(case['depth'])+'-','-d1-')]
            require(exact_trajectory(Path(first['folder']),Path(current['folder']),a.steps),
                    'INTEGER LIFTING DEPTH CHANGED METHOD A RECONSTRUCTION: '+case['name'])
            current['exact_depth_parity']='EXACT_METHOD_A_DEPTH1_STATES_AND_RESIDUAL_TRACES'
        current['exact_parity']='EXACT_ALL_TIMESTEP_STATES_AND_RESIDUAL_TRACES'
        validated[case['name']]=current
        for step in range(a.steps+1):
            reference=read(Path(raw['folder'])/f'potential-state-step{step}.csv')
            value=read(Path(current['folder'])/f'potential-state-step{step}.csv')
            if step==0:require(value['u']==reference['u'],'Initial states differ')
            errors.append(dict(configuration=case['name'],step=step,**error_metrics(value['u'],reference['u'])))
        save()
    if a.hierarchy_profile:
        save()
        print('Setup diagnostics only; no runtime comparison claimed: '+str(root),flush=True)
        return 1 if failed else 0
    active=[c for c in cases if c['name'] in validated]
    rng=random.Random(a.seed)
    for tier in a.tiers:
        repeats=a.detail_repeats if tier==2 else a.repeats
        for repeat in range(-a.warmups,repeats):
            order=list(active);rng.shuffle(order)
            for case in order:run(case,tier,repeat)
    require(digest(a.executable)==executable_hash,'Executable changed during study')
    save()
    for name,rows in [('steps',steps_out),('stages',stages_out)]:
        (root/(name+'.json')).write_text(json.dumps(rows,allow_nan=False)+'\n')
    import summarize_emi_hb_time_study as summarize
    summarize.summarize(root)
    if a.depth_split_study:
        import summarize_emi_hb_depth_split
        summarize_emi_hb_depth_split.summarize(root)
    complete=all(r['status']=='PASS' for r in records) and not failed
    (root/('PASS' if complete else 'FAIL')).write_text('All requested configurations completed.\n' if complete else 'Inspect retained failed configurations; no speed claims for these.\n')
    print('Results: '+str(root),flush=True)
    return 0 if complete else 1

if __name__=='__main__':raise SystemExit(main())
