#!/usr/bin/env python3
"""Three Pepe correctness trajectories, P2/P1/P3, B+R at 16 bits/depth 1/split ON.

One timestep per placement. Diagnostic timings are not a repeated performance
benchmark. Requires an existing two-node/four-rank SLURM allocation; no submission.
"""
import argparse
import json
import os
from pathlib import Path

import run_emi_partition_pilot as pilot
from analyze_emi_potential_diagnostics import read, error_metrics
from collect_emi_scaling_smoke import collect


def arguments(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-only',action='store_true')
    parser.add_argument('--uncompressed',action='store_true',help='Three raw FCG references; no compression or splitting')
    parser.add_argument('--byte-entropy',action='store_true',help='Observe pooled final serialized MPI bytes; diagnostic overhead')
    parser.add_argument('--executable',type=Path,default=Path('emiModel'))
    parser.add_argument('--kaskade',type=Path,default=Path('/data/numerik/people/fchegini/kaskade7_new'))
    parser.add_argument('--slurm-container',type=Path,default=Path('/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif'))
    parser.add_argument('--output',type=Path,default=Path('partition-smoke-results'))
    a=parser.parse_args(argv)
    # Hold all scientific settings fixed; only placement changes in these runs.
    a.case='pepe';a.nodes=2;a.ranks=4;a.steps=1;a.case_minutes=20
    a.intra_nodes='0';a.extra_nodes='1';a.balance_factor=1.0
    return a


def cases(uncompressed=False):
    variant='raw' if uncompressed else 'B-recenter'
    return [c for c in pilot.cases() if c['variant']==variant]


def command(a,case,folder):
    return pilot.command(a,case,folder,validation=True)


def check_placement(folder,case,a):
    ranks=pilot.readcsv(folder/'partition-ranks.csv')
    owners=pilot.readcsv(folder/'partition-ownership.csv')
    if len(ranks)!=a.ranks or len({r['hostname'] for r in ranks})!=a.nodes:
        raise ValueError('Actual physical placement differs from two nodes/four ranks')
    rank_map={int(r['rank']):r for r in ranks}
    if set(rank_map)!=set(range(a.ranks)) or len({r['node_group'] for r in ranks})!=a.nodes:
        raise ValueError('Invalid rank/shared-node mapping')
    host_groups={}
    for r in ranks:host_groups.setdefault(r['node_group'],set()).add(r['hostname'])
    if any(len(h)!=1 for h in host_groups.values()):
        raise ValueError('Shared-memory groups and physical hosts disagree')
    ids=[int(r['subdomain_id']) for r in owners]
    if sorted(ids)!=list(range(len(owners))):raise ValueError('Missing/duplicate subdomain owner')
    for r in owners:
        if r['region'] not in ['I','E']:raise ValueError('Unknown I/E material classification')
        rank=int(r['owner_rank']);node=str(r['node_group'])
        if rank not in rank_map or rank_map[rank]['node_group']!=node:
            raise ValueError('Owner/node metadata disagree')
        if case['policy']=='P2' and rank!=int(r['subdomain_id'])%a.ranks:
            raise ValueError('P2 does not reproduce cyclic ownership')
        if case['policy']=='P3' and node!=('0' if r['region']=='I' else '1'):
            raise ValueError('P3 does not separate I/E onto the requested node pools')
    return ranks,owners


def main(argv=None):
    a=arguments(argv);configurations=cases(a.uncompressed)
    description='Uncompressed FCG, no splitting' if a.uncompressed else '16-bit B+R, depth 1, split ON'
    reference_name=configurations[0]['name']
    if a.plan_only:
        print(json.dumps(dict(status='PLAN_ONLY_NOT_EXECUTED',nodes=a.nodes,ranks=a.ranks,
            steps=a.steps,trajectories=len(configurations),configurations=configurations,
            P3_node_pools=dict(intracellular=[0],extracellular=[1]),
            note='Three one-step raw references' if a.uncompressed else 'Correctness smoke test; no raw reference, timing repetitions or full pilot.'),indent=2))
        return 0
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Smoke execution requires an existing SLURM allocation')
    if (int(os.environ.get('SLURM_JOB_NUM_NODES','0')),int(os.environ.get('SLURM_NTASKS','0')))!=(2,4):
        raise RuntimeError('Smoke allocation must have exactly two nodes and four ranks')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    executable_hash=pilot.digest(a.executable)
    (root/'manifest.json').write_text(json.dumps(dict(vars(a),configurations=configurations,
        executable_sha256=executable_hash,diagnostic_timing_mode=2,
        notes=['Only three one-step trajectories; no performance-speedup claim.',
               'Cross-placement solutions compared with '+reference_name,
               'MPI elapsed time mixes inter-node and same-node buffers.']),default=str,indent=2)+'\n')
    rows=[];reference=None;reference_placement=None;reference_owners=None;state_errors=[]
    for case in configurations:
        folder=root/case['name'];folder.mkdir()
        cfg,cmd=command(a,case,folder)
        row=dict(configuration=case['name'],policy=case['policy'],method='Raw' if a.uncompressed else 'B+R',
                 bits=case['bits'],lifting_depth=None if a.uncompressed else 1,
                 splitting=case['split'].upper(),steps=1,status='FAILED',folder=str(folder))
        rows.append(row)
        (folder/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        for key in list(env):
            if key.startswith(('KASKADE_','APPTAINERENV_KASKADE_','EMI_MESSAGE_ENTROPY_PREFIX')):
                env.pop(key)
        for key,value in [('KASKADE_LIFTING_OPTIMIZATION',str(case['optimization'])),('KASKADE_TRANSFER_TIMING','2'),
                          ('KASKADE_PACKET_TIMING','1'),('KASKADE_TIMESTEP_BYTES','1')]:
            env[key]=value;env['APPTAINERENV_'+key]=value
        print('Running '+case['name']+' (Pepe, FCG, '+description+')',flush=True)
        try:
            rc,timeout,launcher=pilot.hb.study.execute(cmd,folder,cfg,env)
            log=(folder/'run.log').read_text(errors='replace')
            if timeout or rc:
                row['exit_code']=rc
                row['application_log_excerpt']='\n'.join(log.splitlines()[:16])
                raise RuntimeError(('TIMEOUT' if timeout else 'Application exit code '+str(rc))+
                    '; see '+str(folder/'run.log')+'\n'+row['application_log_excerpt'])
            metrics=collect(log,1,a.ranks,rc,1e-7)
            row.update(metrics)
            if timeout or metrics['errors']:
                raise RuntimeError('; '.join(metrics['errors'])+(' TIMEOUT' if timeout else ''))
            row['solver_status']='PASS'
            steps,stages,totals=pilot.hb.profiles(folder,a.ranks,1)
            transport,timing,completed=pilot.observations(folder,a.ranks,1)
            pilot.validate_topology(folder,transport)
            if a.byte_entropy:row.update(pilot.serialized_entropy(folder,transport))
            placement,owners=check_placement(folder,case,a)
            if len(owners)!=metrics['subdomains']:raise ValueError('Owner map does not cover the reported subdomains')
            if any(r['status']!='PASS' or int(r['accepted_updates'])!=metrics['iterations_per_step'][0] for r in completed):
                raise ValueError('Accepted update/completion records disagree with solver summaries')
            sent=sum(int(r['sent_bytes']) for r in transport if r['scope']=='inter_rank')
            if sent!=sum(r['packed_bytes'] for r in steps):raise ValueError('Old/new packet byte counters disagree')
            for scope in ['inter_node','same_node','inter_rank']:
                scoped=[r for r in transport if r['scope']==scope]
                row[scope+'_bytes']=sum(int(r['sent_bytes']) for r in scoped)
                row[scope+'_messages']=sum(int(r['mpi_messages']) for r in scoped)
            for phase in ['restriction','prolongation']:
                scoped=[r for r in transport if r['phase']==phase and r['scope']=='inter_rank']
                applications={int(r['phase_applications']) for r in scoped}
                if len(applications)!=1:raise ValueError('Ranks disagree on phase application counts')
                row[phase+'_applications']=applications.pop()
                row[phase+'_bytes']=sum(int(r['sent_bytes']) for r in scoped)
            row.update(accepted_iterations=metrics['iterations_per_step'][0],
                true_relative_residual=metrics['relative_residuals'][0],
                diagnostic_application_s=max(t['application_compute_s'] for t in totals),
                diagnostic_fcg_s=max(sum(r['fcg_s'] for r in steps if r['rank']==rank) for rank in range(a.ranks)),
                diagnostic_mixed_payload_collective_s=max(sum(float(r['mpi_payload_collective_s']) for r in timing if r['rank']==rank) for rank in range(a.ranks)),
                owner_map_sha256=pilot.digest(folder/'partition-ownership.csv'))
            invariants=pilot.readcsv(folder/'partition-problem.csv')[0]
            invariants={k:v for k,v in invariants.items() if k not in ['mode','planning_s']}
            invariants['coarse_rows']=(folder/'partition-coarse.csv').read_text()
            if case['policy']=='P2':
                reference=invariants;reference_placement=placement
                reference_owners=[(r['subdomain_id'],r['owner_rank']) for r in owners]
            elif reference is not None:
                if placement!=reference_placement:raise ValueError('Rank/host/CPU placement changed across runs')
                for step in [0,1]:
                    values=read(folder/f'potential-state-step{step}.csv')['u']
                    baseline=read(root/reference_name/f'potential-state-step{step}.csv')['u']
                    if step==0 and values!=baseline:raise ValueError('Initial states differ')
                    state_errors.append(dict(configuration=case['name'],step=step,
                        reference=reference_name,**error_metrics(values,baseline)))
                if invariants!=reference:
                    row['status']='REVIEW_REQUIRED'
                    row['error']='Operator/metadata fingerprint differs from P2; inspect roundoff/reduction order before interpreting'
                else:row['status']='PASS'
            else:
                row['status']='REVIEW_REQUIRED';row['error']='P2 reference unavailable; independent solve converged'
            if case['policy']=='P2':row['status']='PASS'
            if reference_owners is not None:
                row['same_owners_as_P2']=[(r['subdomain_id'],r['owner_rank']) for r in owners]==reference_owners
            print('  '+row['status']+': iterations='+str(row['accepted_iterations'])+
                  ', true relative residual='+str(row['true_relative_residual']),flush=True)
        except Exception as error:
            row['status']='FAILED';row['error']=str(error)
            print('  FAILED: '+str(error),flush=True)
        (root/'runs.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
        pilot.writecsv(root/'comparison.csv',rows)
        pilot.writecsv(root/'state-errors-vs-P2.csv',state_errors)
    if pilot.digest(a.executable)!=executable_hash:raise RuntimeError('Executable changed during check')
    (root/'report.txt').write_text('Pepe partition correctness check: 2 nodes, 4 ranks, 1 timestep.\n'
        'FCG; '+description+'.\n'
        'P2 cyclic; P1 coupled; P3 intracellular node group 0 / extracellular group 1.\n'
        'Exactly three trajectories. Cross-placement errors use '+reference_name+'.\n'
        'No repeated timings, no speedup/accuracy-equivalence claim.\n'
        'Elapsed MPI collectives mix inter-node and same-node communication.\n'
        'Inspect comparison.csv, state-errors-vs-P2.csv and per-run residual/rho/topology records.\n'+
        '\n'.join(r['configuration']+': '+r['status'] for r in rows)+'\n')
    print('Results: '+str(root),flush=True)
    return int(any(r['status']!='PASS' for r in rows))


if __name__=='__main__':raise SystemExit(main())
