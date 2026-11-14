#!/usr/bin/env python3
"""Observed paired trajectory comparisons and explicitly conditional break-even models."""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import statistics as stats


def mean_ci(values):
    if not values:return dict(n=0,mean=None,sd=None,low=None,high=None)
    mean=stats.mean(values);n=len(values)
    if n<2:return dict(n=n,mean=mean,sd=None,low=None,high=None)
    # Two-sided 95% Student t, conservative tabulated bounds for unlisted df.
    bounds=[(1,12.7063),(2,4.3027),(3,3.1825),(4,2.7765),(5,2.5706),
            (6,2.4469),(7,2.3646),(8,2.3061),(9,2.2622),(10,2.2282),
            (15,2.1315),(20,2.0860),(30,2.0423),(60,2.0003)]
    df=n-1;critical=next((v for d,v in reversed(bounds) if d<=df),12.7063)
    sd=stats.stdev(values);margin=critical*sd/math.sqrt(n)
    return dict(n=n,mean=mean,sd=sd,low=mean-margin,high=mean+margin)


def projected_break_even(startup_difference,steady_gain):
    """First persistent amortized benefit under the stationary cost model."""
    if steady_gain<0:return None
    if steady_gain==0:return 1 if startup_difference<0 else None
    return max(1,math.floor(startup_difference/steady_gain)+1)


def write_csv(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def summarize(root):
    root=Path(root);manifest=json.loads((root/'manifest.json').read_text())
    runs=json.loads((root/'runs.json').read_text())
    if (root/'steps.json').is_file() and (root/'stages.json').is_file():
        steps=json.loads((root/'steps.json').read_text())
        stages=json.loads((root/'stages.json').read_text())
    else:
        # Recover completed runs without waiting for the final aggregation.
        # Raw per-rank files are written before a run is marked PASS.
        from run_emi_hb_time_study import profiles
        steps=[];stages=[]
        for run in runs:
            if run['status']!='PASS':continue
            rows,detail,_=profiles(Path(run['folder']),manifest['ranks'],manifest['steps'])
            common={k:run[k] for k in ['configuration','tier','repeat','validation','warmup']}
            common['run']=Path(run['folder']).name
            steps.extend(dict(common,**r) for r in rows)
            stages.extend(dict(common,**r) for r in detail)
    errors=json.loads((root/'errors.json').read_text())
    count=manifest['steps'];repeats=manifest['repeats']
    measured=[r for r in runs if not r['validation'] and not r['warmup'] and r['status']=='PASS']
    byrun=defaultdict(list)
    for row in steps:byrun[row['run']].append(row)
    stagebyrun=defaultdict(list)
    for row in stages:
        if row['call']==0:stagebyrun[row['run']].append(row)

    def run_name(row):return Path(row['folder']).name
    def step_rows(row,limit=None):
        return [r for r in byrun[run_name(row)] if r['step']>0 and (limit is None or r['step']<=limit)]
    def maxsum(rows,key):
        groups=defaultdict(float)
        for r in rows:groups[r['rank']]+=r[key]
        return max(groups.values(),default=0.)
    def prefix(row,n):
        return maxsum([r for r in byrun[run_name(row)] if r['step']<=n],'wall_s')
    def stage_time(row,stage,step_filter=lambda s:s>0):
        return maxsum([r for r in stagebyrun[run_name(row)] if r['stage']==stage and r['level']==0 and step_filter(r['step'])],'seconds')
    def steady(row):return maxsum([r for r in step_rows(row) if r['step']>1],'wall_s')/(count-1)

    # Paired original-versus-optimized trajectories within the same campaign.
    # IDs/order are randomized together, so repeats share an allocation.
    controls={(r['configuration'],r['tier'],r['repeat']):r for r in measured}
    optimized=[]
    for row in measured:
        if not row.get('optimization',0) or not row.get('baseline_name'):continue
        baseline=controls.get((row['baseline_name'],row['tier'],row['repeat']))
        if baseline is None:continue
        optimized.append(dict(configuration=row['configuration'],tier=row['tier'],repeat=row['repeat'],
            optimization=row['optimization'],
            simulation_difference_s=row['simulation_s']-baseline['simulation_s'],
            application_difference_s=row['application_compute_s']-baseline['application_compute_s'],
            setup_difference_s=prefix(row,0)-prefix(baseline,0),
            peak_rss_difference_bytes=(row['process_peak_rss_bytes']-baseline['process_peak_rss_bytes']
                if row.get('process_peak_rss_bytes') is not None and baseline.get('process_peak_rss_bytes') is not None else None),
            packed_bytes_equal=row['packed_bytes']==baseline['packed_bytes'],
            iterations_equal=row['iterations_per_step']==baseline['iterations_per_step']))
    write_csv(root/'optimization-paired-runs.csv',optimized)
    optgroups=defaultdict(list)
    for row in optimized:optgroups[row['configuration'],row['tier']].append(row)
    comparisons=[]
    for (name,tier),rows in sorted(optgroups.items()):
        out=dict(configuration=name,tier=tier,n=len(rows),optimization=rows[0]['optimization'])
        for key in ['simulation_difference_s','application_difference_s','setup_difference_s','peak_rss_difference_bytes']:
            for statistic,value in mean_ci([r[key] for r in rows if r[key] is not None]).items():out[key+'_'+statistic]=value
        out['packed_bytes_equal']=all(r['packed_bytes_equal'] for r in rows)
        out['iterations_equal']=all(r['iterations_equal'] for r in rows)
        comparisons.append(out)
    write_csv(root/'optimization-comparison.csv',comparisons)
    summaries=[];prefixes=[];components=[]
    groups=defaultdict(list)
    for row in measured:groups[row['configuration'],row['tier']].append(row)
    raw={(r['tier'],r['repeat']):r for r in measured if r['configuration']=='raw'}
    for (name,tier),selected in sorted(groups.items()):
        valid=[r for r in selected if (tier,r['repeat']) in raw]
        header=dict(configuration=name,tier=tier,n=len(valid),expected_n=manifest['detail_repeats'] if tier==2 else repeats,
                    bits=selected[0]['bits'],depth=selected[0]['depth'],variant=selected[0]['variant'],
                    split=selected[0].get('split','off'),optimization=selected[0].get('optimization',0))
        header['complete']=header['n']==header['expected_n']
        metrics={
            'trajectory_s':[r['simulation_s'] for r in valid],
            'application_s':[r['application_compute_s'] for r in valid],
            'reported_application_wall_s':[r['wall_seconds'] for r in valid],
            'process_peak_rss_bytes':[r['process_peak_rss_bytes'] for r in valid if r.get('process_peak_rss_bytes') is not None],
            'model_setup_s':[r['application_setup_s'] for r in valid],
            'fcg_s':[maxsum(step_rows(r),'fcg_s') for r in valid],
            'setup_s':[prefix(r,0) for r in valid],
            'per_timestep_difference_s':[(maxsum(step_rows(r),'wall_s')-maxsum(step_rows(raw[tier,r['repeat']]),'wall_s'))/count for r in valid],
            'first_step_s':[maxsum([s for s in step_rows(r) if s['step']==1],'wall_s') for r in valid],
            'steady_step_s':[steady(r) for r in valid],
            'iterations':[sum(r['iterations_per_step']) for r in valid],
            'packed_bytes':[r['packed_bytes'] for r in valid],
            'saving_percent':[100*(1-r['packed_bytes']/raw[tier,r['repeat']]['packed_bytes']) for r in valid],
            'inter_node_bytes':[sum(s['inter_node_packed_bytes'] for s in step_rows(r)) for r in valid],
            'inter_node_saving_percent':[100*(1-sum(s['inter_node_packed_bytes'] for s in step_rows(r))/sum(s['inter_node_packed_bytes'] for s in step_rows(raw[tier,r['repeat']]))) for r in valid
                if sum(s['inter_node_packed_bytes'] for s in step_rows(raw[tier,r['repeat']]))],
            'trajectory_difference_s':[r['simulation_s']-raw[tier,r['repeat']]['simulation_s'] for r in valid],
            'application_difference_s':[r['application_compute_s']-raw[tier,r['repeat']]['application_compute_s'] for r in valid],
            'fcg_per_call_s':[maxsum(step_rows(r),'fcg_s')/sum(r['iterations_per_step']) for r in valid if sum(r['iterations_per_step'])],
            'iteration_ratio':[sum(r['iterations_per_step'])/sum(raw[tier,r['repeat']]['iterations_per_step']) for r in valid if sum(raw[tier,r['repeat']]['iterations_per_step'])],
        }
        if tier>0:
            for stage in ['compression','decompression','mpi','plan_setup','codebook_setup',
                          'float_forward','float_inverse','integer_forward','integer_inverse',
                          'quantize','dequantize','huffman_encode','huffman_decode','recenter','uncenter','packet_exchange']:
                metrics[stage+'_s']=[stage_time(r,stage) for r in valid]
            for key in ['pack_s','unpack_s','exchange_s','verify_s']:
                metrics[key]=[maxsum(step_rows(r),key) for r in valid]
            metrics['hierarchy_setup_s']=[stage_time(r,'plan_setup',lambda s:s==0) for r in valid]
        for key,values in metrics.items():
            for suffix,value in mean_ci(values).items():header[key+'_'+suffix]=value
        header['max_trajectory_relative_l2_vs_raw']=max((e['relative_l2'] for e in errors if e['configuration']==name and e['relative_l2'] is not None),default=None)
        header['final_trajectory_relative_l2_vs_raw']=next((e['relative_l2'] for e in errors if e['configuration']==name and e['step']==count),None)
        summaries.append(header)

        if tier==0 and name!='raw' and header['complete']:
            startup=[];gain=[]
            for r in valid:
                b=raw[tier,r['repeat']]
                tail=steady(r)-steady(b)
                startup.append(prefix(r,1)-prefix(b,1)-tail)
                gain.append(-tail)
            gain_stats=mean_ci(gain)
            header['conditional_break_even_steps']=projected_break_even(stats.mean(startup),gain_stats['mean'])
            header['steady_gain_s_mean']=gain_stats['mean'];header['steady_gain_s_ci_low']=gain_stats['low'];header['steady_gain_s_ci_high']=gain_stats['high']
            header['steady_gain_resolved']=gain_stats['low'] is not None and gain_stats['low']>0
            core_gains=[];core_startup=[]
            for r in valid:
                b=raw[tier,r['repeat']]
                ci=sum(r['iterations_per_step'][1:]);ri=sum(b['iterations_per_step'][1:])
                if not ci or not ri:continue
                ct=maxsum([s for s in step_rows(r) if s['step']>1],'fcg_s')/ci
                rt=maxsum([s for s in step_rows(b) if s['step']>1],'fcg_s')/ri
                gain_per_raw_call=rt-(ci/ri)*ct
                first_delta=prefix(r,0)-prefix(b,0)+maxsum([s for s in step_rows(r) if s['step']==1],'fcg_s')-maxsum([s for s in step_rows(b) if s['step']==1],'fcg_s')
                core_gains.append(gain_per_raw_call)
                core_startup.append(first_delta+b['iterations_per_step'][0]*gain_per_raw_call)
            core=mean_ci(core_gains)
            header['fcg_gain_per_raw_call_s_mean']=core['mean']
            header['fcg_gain_per_raw_call_s_ci_low']=core['low']
            header['conditional_fcg_raw_calls_break_even']=projected_break_even(stats.mean(core_startup),core['mean']) if core_startup else None

            header['effective_startup_difference_s_mean']=stats.mean(startup)
            observed=[]
            for n in range(1,count+1):
                delta=mean_ci([prefix(r,n)-prefix(raw[tier,r['repeat']],n) for r in valid])
                entry=dict(configuration=name,steps=n,**delta)
                entry['compressed_fcg_calls_mean']=stats.mean(sum(r['iterations_per_step'][:n]) for r in valid)
                entry['raw_fcg_calls_mean']=stats.mean(sum(raw[tier,r['repeat']]['iterations_per_step'][:n]) for r in valid)
                prefixes.append(entry);observed.append(entry)
            header['observed_mean_benefit_from_step']=next((r['steps'] for i,r in enumerate(observed) if all(x['mean']<0 for x in observed[i:])),None)
            header['observed_ci_benefit_from_step']=next((r['steps'] for i,r in enumerate(observed) if all(x['high'] is not None and x['high']<0 for x in observed[i:])),None)

    # Exchange averages are rank-local per phase and FCG call. Aggregate stage
    # counters count nested scalar/packet operations, NOT network messages.
    stage_groups=defaultdict(list)
    stage_index={}
    for run in measured:
        if run['tier']==0:continue
        for r in stagebyrun[run_name(run)]:
            stage_index[(run_name(run),r['step'],r['rank'],r['phase'],r['stage'],r['level'])]=r
            stage_groups[(run['configuration'],run['tier'],r['step'],r['rank'],r['phase'],r['stage'],r['level'])].append(r)
    for key,rows in stage_groups.items():
        name,tier,step,rank,phase,stage,level=key
        entry=dict(configuration=name,tier=tier,step=step,rank=rank,phase=phase,stage=stage,level=level)
        timings=[];calls=[];per_exchange=[]
        for run in groups[name,tier]:
            key=(run_name(run),step,rank,phase,stage,level)
            current=stage_index.get(key,{})
            timings.append(current.get('seconds',0.));calls.append(current.get('calls',0))
            exchange=stage_index.get((run_name(run),step,rank,phase,'packet_exchange',0),{})
            if exchange.get('calls'):per_exchange.append(current.get('seconds',0.)/exchange['calls'])
        for suffix,value in mean_ci(timings).items():entry['seconds_'+suffix]=value
        entry['operation_calls_mean']=stats.mean(calls)
        entry['seconds_per_phase_exchange_mean']=stats.mean(per_exchange) if per_exchange else None
        components.append(entry)
    combined=[]
    for primary in summaries:
        if primary['tier']!=0:continue
        row=dict(primary)
        for tier,prefix_label in [(1,'aggregate_'),(2,'detailed_')]:
            match=next((r for r in summaries if r['configuration']==primary['configuration'] and r['tier']==tier),None)
            if match:
                row.update({prefix_label+k:v for k,v in match.items() if k not in ['configuration','variant','bits','depth']})
        combined.append(row)
    timestep=[]
    for (name,tier),selected in sorted(groups.items()):
        for step in range(1,count+1):
            valid=[r for r in selected if (tier,r['repeat']) in raw]
            row=dict(configuration=name,tier=tier,step=step)
            values={
                'wall_s':[maxsum([s for s in step_rows(r) if s['step']==step],'wall_s') for r in valid],
                'fcg_s':[maxsum([s for s in step_rows(r) if s['step']==step],'fcg_s') for r in valid],
                'iterations':[r['iterations_per_step'][step-1] for r in valid],
                'wall_difference_s':[maxsum([s for s in step_rows(r) if s['step']==step],'wall_s')-maxsum([s for s in step_rows(raw[tier,r['repeat']]) if s['step']==step],'wall_s') for r in valid],
                'packed_bytes':[sum(s['packed_bytes'] for s in step_rows(r) if s['step']==step) for r in valid],
            }
            if tier>0:
                for stage in ['compression','decompression','mpi','packet_exchange','codebook_setup']:
                    values[stage+'_s']=[stage_time(r,stage,lambda s:s==step) for r in valid]
            for key,data in values.items():
                for suffix,value in mean_ci(data).items():row[key+'_'+suffix]=value
            timestep.append(row)
    write_csv(root/'combined-comparison.csv',combined)
    write_csv(root/'timestep-comparison.csv',timestep)
    write_csv(root/'comparison.csv',summaries);write_csv(root/'observed-prefixes.csv',prefixes)
    write_csv(root/'components-by-step-level-rank.csv',components)
    lines=['Fixed-topology EMI time study',
        f"Case {manifest['case']}; {manifest['nodes']} nodes / {manifest['ranks']} ranks; {count} timesteps.",
        'Negative paired differences mean compression is faster. Mode 0 is primary.',
        'All intervals are exploratory paired 95% Student-t intervals; no multiple-comparison adjustment.',
        'Detailed scalar timers perturb kernels. Operation/call and step timings are nested.',
        'Trajectory timer covers retained-BDDC setup and timesteps; application compute timer includes mesh/matrix setup and excludes profiler-file I/O and final reporting.',
        'Reported application/launcher walls are retained separately.',
        'Projected break-even assumes the measured post-first-step cost stays stationary; only observed prefixes are measured.',
        'For gain <= 0, more timesteps cannot amortize a positive setup penalty.',
        'Hierarchy setup is one-time; first codebook training and packet-plan setup remain charged to the first timestep.',
        'Per-step wall includes required true-residual assembly/broadcast; solve-call wall excludes caller checks.',
        'Same-node MPI bytes are excluded only in inter_node_bytes, not in packed_bytes.',
        'No bandwidth-only wire-speed claim is inferred from routed application bytes.',
        'Reconstruction column is trajectory state error vs uncompressed FCG, not transform-only error.',
        'Recommended configuration requires complete exact-parity validation, acceptable state error and a resolved total-runtime improvement.',
        'Adaptive per-packet recentering is not part of this experiment.']
    primary=[r for r in summaries if r['tier']==0 and r['configuration']!='raw' and r['complete']]
    winners=[r for r in primary if r['trajectory_difference_s_high']<0 and r['application_difference_s_high']<0]
    if winners:
        best=min(winners,key=lambda r:r['trajectory_s_mean'])
        lines.append('Fastest resolved trajectory/application benefit: '+best['configuration']+'; review reported trajectory error before acceptance.')
    else:lines.append('No configuration demonstrates a resolved improvement in both trajectory and application runtime.')
    for r in primary:
        lines.append(f"{r['configuration']}: bytes={r['packed_bytes_mean']:.0f}, savings={r['saving_percent_mean']:.3f}%, "
                     f"trajectory={r['trajectory_s_mean']:.6g} +/- {r['trajectory_s_sd']:.3g} s, "
                     f"delta={r['trajectory_difference_s_mean']:.6g} s, "
                     f"iterations={r['iterations_mean']:.3g}, projected N={r.get('conditional_break_even_steps')}, "
                     f"steady gain resolved={r.get('steady_gain_resolved')}, observed CI threshold={r.get('observed_ci_benefit_from_step')}")
    failures=[r for r in runs if r['status']!='PASS']
    lines.append(f'Failed trajectories retained: {len(failures)}; see runs.json.')
    (root/'report.txt').write_text('\n'.join(lines)+'\n')
    return summaries


def scaling(roots,output):
    rows=[]
    for root in roots:
        m=json.loads((root/'manifest.json').read_text())
        runs=json.loads((root/'runs.json').read_text())
        metadata=next((r for r in runs if r['status']=='PASS'),{})
        for row in summarize(root):
            if row['tier']!=0:continue
            rows.append(dict(root=str(root),scaling=m['scaling'],mesh=m['case'],nodes=m['nodes'],ranks=m['ranks'],steps=m['steps'],
                dofs=metadata.get('dofs'),subdomains=metadata.get('subdomains'),dofs_per_rank=metadata.get('dofs',0)/m['ranks'],subdomains_per_rank=metadata.get('subdomains',0)/m['ranks'],
                owned_dofs_rank_mean=metadata.get('owned_dofs_rank_mean'),owned_dofs_rank_min=metadata.get('owned_dofs_rank_min'),owned_dofs_rank_max=metadata.get('owned_dofs_rank_max'),**row))
    for r in rows:
        matching=[b for b in rows if b['configuration']==r['configuration'] and b['scaling']==r['scaling'] and b['steps']==r['steps']
                  and (b['mesh']==r['mesh'] if r['scaling']=='strong' else True)]
        base=min(matching,key=lambda b:b['ranks'])
        r['base_ranks']=base['ranks']
        r['trajectory_speedup']=base['trajectory_s_mean']/r['trajectory_s_mean']
        r['scaling_efficiency']=r['trajectory_speedup']*(base['ranks']/r['ranks'] if r['scaling']=='strong' else 1)
        r['work_per_rank_ratio_vs_base']=r['dofs_per_rank']/base['dofs_per_rank'] if base['dofs_per_rank'] else None
    write_csv(output,rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('roots',nargs='+',type=Path)
    p.add_argument('--scaling-output',type=Path);a=p.parse_args()
    if a.scaling_output:scaling(a.roots,a.scaling_output)
    else:
        for root in a.roots:summarize(root)
