#!/usr/bin/env python3
"""Depth/split comparisons. Primary runtime is uninstrumented mode 0.
Per-level profiles are setup diagnostics, not geometric face counts or kernel times.
"""
from collections import defaultdict
import csv
import json
from pathlib import Path
import statistics as stats
from summarize_emi_hb_time_study import mean_ci,write_csv


def inventory(root,runs,required=False):
    rows=[];plans=[]
    for run in runs:
        if not run['validation'] or not run.get('reuse',False) or run['status']!='PASS' or run['variant']=='raw':continue
        files=list(Path(run['folder']).glob('lifting-setup*.csv'))
        if required and not files:raise RuntimeError('Missing interface inventory: '+run['folder'])
        for file in files:
            groups=defaultdict(list)
            with file.open() as f:
                for r in csv.DictReader(f):groups[r['plan_id']].append(r)
            for pid,events in groups.items():
                final=next((r for r in events if r['kind']=='graph_plan'),None)
                first=next((r for r in events if r['kind']=='predictor_split'),None)
                if not final or not first:continue
                stages=[first]+sorted((r for r in events if r['kind']=='elimination' and int(r['detail'])>0),key=lambda r:int(r['level']))
                n=int(first['nodes_before']);requested=int(first['requested_depth']);actual=int(final['level'])
                common=dict(configuration=run['configuration'],source=file.name,plan_id=pid,
                    producer=int(first['producer']),consumer=int(first['consumer']),role=first['role'],phase=first['phase'],
                    requested_depth=requested,actual_depth=actual,interface_dofs=int(first['interface_dofs']),
                    exact_dofs=int(first['exact_dofs']),group_graph_dofs=n)
                # Every successful prediction must partition its active coefficients.
                previous=n
                for r in stages:
                    c,f=int(r['coarse']),int(r['detail']);active=c+f
                    if active!=previous:raise RuntimeError('Inconsistent hierarchy partition '+str(common))
                    rows.append(dict(common,level=int(r['level']),active_dofs=active,coarse_dofs=c,detail_dofs=f,
                        coarse_fraction_of_active=c/active if active else None,
                        coarse_fraction_of_original=c/n if n else None,
                        edges=int(r['edges_after']),max_degree=int(r['max_degree']),
                        setup_component_s=float(r['seconds']),process_peak_rss_bytes=int(r['process_peak_rss_bytes'])))
                    previous=c
                if len(stages)!=actual:raise RuntimeError('Missing actual hierarchy stages '+str(common))
                plans.append(dict(common,final_coarse_dofs=previous,final_coarse_fraction=previous/n if n else None,
                    reached_requested_depth=actual>=requested,construction_s=float(final['seconds']),
                    plan_storage_bytes=int(final['plan_storage_bytes'])))
    write_csv(root/'interface-level-shrinkage.csv',rows)
    write_csv(root/'interface-plan-inventory.csv',plans)
    groups=defaultdict(list)
    for r in plans:groups[r['configuration'],r['role']].append(r)
    summary=[]
    for (name,role),ps in groups.items():
        summary.append(dict(configuration=name,role=role,constructed_group_plans=len(ps),
            requested_depth=ps[0]['requested_depth'],reached_requested=sum(p['reached_requested_depth'] for p in ps),
            actual_depth_min=min(p['actual_depth'] for p in ps),actual_depth_max=max(p['actual_depth'] for p in ps),
            group_dofs_min=min(p['group_graph_dofs'] for p in ps),group_dofs_max=max(p['group_graph_dofs'] for p in ps),
            final_coarse_fraction_mean=stats.mean(p['final_coarse_fraction'] for p in ps if p['final_coarse_fraction'] is not None)))
    write_csv(root/'interface-depth-summary.csv',summary)
    return plans


def summarize(root):
    root=Path(root);m=json.loads((root/'manifest.json').read_text());runs=json.loads((root/'runs.json').read_text())
    steps=json.loads((root/'steps.json').read_text());stages=json.loads((root/'stages.json').read_text())
    with (root/'combined-comparison.csv').open() as f:combined={r['configuration']:r for r in csv.DictReader(f)}
    measured=[r for r in runs if r['status']=='PASS' and not r['validation'] and not r['warmup']]
    valid={r['configuration']:r for r in runs if r['validation'] and r.get('exact_parity')}
    controls={(r['configuration'],r['tier'],r['repeat']):r for r in measured}
    byrun=defaultdict(list);sg=defaultdict(list)
    for r in steps:byrun[r['run']].append(r)
    for r in stages:
        if r['call']==0 and r['step']>0:sg[r['run']].append(r)
    groups=defaultdict(list)
    for r in measured:groups[r['configuration'],r['tier']].append(r)
    def maxsum(rows,key):
        byrank=defaultdict(float)
        for r in rows:byrank[r['rank']]+=r[key]
        return max(byrank.values(),default=0.)
    def samples(out,key,values):
        for s,v in mean_ci(values).items():out[key+'_'+s]=v
    table=[];phase=[];paired=[]
    for (name,tier),selected in sorted(groups.items()):
        if tier!=0:continue
        c=selected[0];generic=combined[name]
        row=dict(configuration=name,method=c['variant'],bits=c['bits'],depth_cap=c['depth'],split=c.get('split','off'),
            recenter=c['variant']=='B-recenter',optimization=c.get('optimization',0),
            exact_implementation_parity=valid.get(name,{}).get('exact_parity','raw control'),
            iterations_per_timestep=json.dumps(c['iterations_per_step']),
            exported_remote_plan_bytes=sum(v['exported'] for v in valid.get(name,{}).get('hierarchy_setup_payload',[])),
            max_state_relative_l2_vs_raw=generic['max_trajectory_relative_l2_vs_raw'])
        for key in ['packed_bytes','saving_percent','iterations','trajectory_s','application_s','fcg_s','setup_s','steady_step_s','process_peak_rss_bytes','trajectory_difference_s','application_difference_s']:
            for suffix in ['mean','sd','low','high']:row[key+'_'+suffix]=generic.get(key+'_'+suffix)
        samples(row,'bytes_per_iteration',[r['packed_bytes']/sum(r['iterations_per_step']) for r in selected])
        # Mode 1 timing is separate; never add independent rank maxima to primary wall.
        detail=groups.get((name,1),[])
        for key in ['hierarchy_setup_s','compression_s','decompression_s','mpi_s','packet_exchange_s','exchange_s','recenter_s','uncenter_s']:
            for suffix in ['mean','sd']:row['mode1_'+key+'_'+suffix]=generic.get('aggregate_'+key+'_'+suffix)
        table.append(row)
        for direction in ['restrict','prolongate']:
            entry=dict(configuration=name,phase=direction)
            phase_values=[];raw_values=[]
            for r in selected:
                rr=byrun[Path(r['folder']).name]
                raw=controls['raw',0,r['repeat']]
                rawrr=byrun[Path(raw['folder']).name]
                phase_values.append(sum(x[direction+'_packed_bytes'] for x in rr))
                raw_values.append(sum(x[direction+'_packed_bytes'] for x in rawrr))
            samples(entry,'packed_bytes',phase_values)
            samples(entry,'saving_percent',[100*(1-v/b) for v,b in zip(phase_values,raw_values) if b])
            samples(entry,'bytes_per_iteration',[v/sum(r['iterations_per_step']) for v,r in zip(phase_values,selected)])
            for stage in ['compression','decompression','mpi','packet_exchange']:
                values=[maxsum([s for s in sg[Path(r['folder']).name] if s['phase']==direction and s['stage']==stage and s['level']==0],'seconds') for r in detail]
                if values:samples(entry,'mode1_'+stage+'_s',values)
            phase.append(entry)
        for kind,reference in [('raw','raw'),('recenter',name.replace('B-recenter','B')),('split',name.replace('-split-on','-split-off')),('depth1',name.replace('-d'+str(c['depth'])+'-','-d1-'))]:
            if reference==name:continue
            if kind=='recenter' and c['variant']!='B-recenter':continue
            if kind=='split' and c.get('split')!='on':continue
            if kind=='depth1' and c['depth']==1:continue
            current=[];baseline=[]
            for r in selected:
                b=controls.get((reference,0,r['repeat']))
                if b:current.append(r);baseline.append(b)
            if not current:continue
            out=dict(configuration=name,reference=reference,comparison=kind)
            for key in ['simulation_s','application_compute_s','packed_bytes']:
                samples(out,key+'_difference',[r[key]-b[key] for r,b in zip(current,baseline)])
            samples(out,'iteration_difference',[sum(r['iterations_per_step'])-sum(b['iterations_per_step']) for r,b in zip(current,baseline)])
            paired.append(out)
    write_csv(root/'depth-split-comparison.csv',table)
    write_csv(root/'phase-comparison.csv',phase)
    write_csv(root/'depth-split-paired.csv',paired)
    plans=inventory(root,runs,m.get('interface_inventory',False))
    failed=[r for r in runs if r['status']!='PASS']
    lines=[f"EMI depth/splitting study: {m['case']}, {m['nodes']} nodes/{m['ranks']} ranks, {m['steps']} timesteps, bits={m['bits']}.",
        'Depth cap applies to BOTH R/P; splitting applies only to prolongation. A has no recentering switch.',
        'Runtime: mode 0 mean +/- sample SD; exploratory paired 95% t intervals in depth-split-paired.csv.',
        'State error is against raw FCG; exact quantized-code/reconstruction parity is a separate native test and live B/recenter validation.',
        'Bytes include routed application packets (also same-node inter-rank MPI), metadata; exported prediction-plan setup bytes are separate; not NIC wire bytes.',
        'phase-comparison.csv splits routed R/P bytes and nested mode-1 codec/MPI times. Never add separate rank maxima.',
        'Interface inventory identifies constructed owner-group graphs, which can repeat across union/peer cache entries; these are not unique geometric faces.',
        'Inventory runs have setup CSV I/O and are excluded from runtime statistics. Actual successful levels exclude unsuccessful terminal splits.',
        f'Failed retained runs: {len(failed)}; graph constructions inventoried: {len(plans)}.']
    for r in table:
        lines.append(f"{r['configuration']}: savings={float(r['saving_percent_mean']):.3f}%, bytes={float(r['packed_bytes_mean']):.0f}, iterations={r['iterations_per_timestep']}, bytes/iteration={r['bytes_per_iteration_mean']:.1f}, application={float(r['application_s_mean']):.3f} +/- {float(r['application_s_sd'] or 0):.3f} s, max state error={r['max_state_relative_l2_vs_raw']}")
    compressed=[r for r in table if r['method']!='raw']
    if compressed:
        fastest=min(compressed,key=lambda r:float(r['application_s_mean']))
        smallest=min(compressed,key=lambda r:float(r['packed_bytes_mean']))
        lines+=['Fastest observed compressed configuration: '+fastest['configuration']+'.',
                'Fewest full-trajectory bytes: '+smallest['configuration']+'.',
                'Select a larger-mesh candidate only after reviewing state-error requirements, exact parity, paired runtime intervals and interface shrinkage. No automatic accuracy threshold or scaling claim.']
    (root/'depth-split-report.txt').write_text('\n'.join(lines)+'\n')
    return table

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('results',type=Path);summarize(p.parse_args().results)
