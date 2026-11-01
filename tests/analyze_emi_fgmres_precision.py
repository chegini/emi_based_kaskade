#!/usr/bin/env python3
"""Read-only EMI campaign analysis. Writes a new report directory, never solver data.

Dependencies: numpy, matplotlib. Example:
python3 tests/analyze_emi_fgmres_precision.py --campaign output/emi-fgmres-...
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analyze_emi_potential_diagnostics import read, error_metrics

DEFAULTS = (
    "emi-fgmres-20261101T200034Z-14562",
    "emi-precision-followup-20261101T202243Z-14983/b4-cap500",
    "emi-fgmres-20261101T203933Z-15690",
    "emi-precision-followup-20261101T202243Z-14983/b6-ten-steps",
    "emi-fgmres-20261101T151610Z-14189",
    "emi-fgmres-20261101T151946Z-14375",
)
LABELS = {"uncompressed-fgmres":"Uncompressed", "full-fgmres16":"16/16",
          "full-fgmres8-split":"8/8 + split", "full-fgmres6-split":"6/6 + split",
          "full-fgmres4-split":"4/4 + split", "full-cg16":"CG 16/16"}
COLORS = {"uncompressed-fgmres":"#444444","full-fgmres16":"#0072B2",
          "full-fgmres8-split":"#009E73","full-fgmres6-split":"#E69F00",
          "full-fgmres4-split":"#D55E00","full-cg16":"#CC79A7"}
for pair,color in zip(((8,6),(6,8),(8,4),(4,8)),('#882255','#44AA99','#332288','#AA4499')):
    r,p=pair;mode=f'full-fgmres-r{r}-p{p}-split'
    LABELS[mode]=f'R{r}/P{p} + split';COLORS[mode]=color

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()

def write_csv(path, rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("x",newline="") as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def trend(rows, start=None):
    selected=[r for r in rows if start is None or int(r["iteration"])>=start]
    k=np.array([int(r["iteration"]) for r in selected])
    y=np.array([float(r["true_relative_residual"]) for r in selected])
    if len(y)<2:raise ValueError("Need at least two residual observations")
    return dict(first_iteration=int(k[0]),last_iteration=int(k[-1]),observations=len(y),
                first_residual=float(y[0]),last_residual=float(y[-1]),
                end_over_start=float(y[-1]/y[0]),
                increases=int(np.sum(y[1:]>y[:-1]*(1+1e-12))),
                best_residual=float(y.min()),
                log10_slope=float(np.polyfit(k,np.log10(np.maximum(y,1e-300)),1)[0]))

def table(rows, keys):
    return "<table><thead><tr>"+"".join("<th>"+html.escape(k)+"</th>" for k in keys)+"</tr></thead><tbody>"+"".join(
        "<tr>"+"".join("<td>"+html.escape(fmt(r.get(k)))+"</td>" for k in keys)+"</tr>" for r in rows)+"</tbody></table>"

def fmt(v):
    if v is None:return "not measured"
    if isinstance(v,float):return f"{v:.6g}"
    return str(v)

def git_id(repo):
    r=subprocess.run(["git","rev-parse","HEAD"],cwd=repo,capture_output=True,text=True)
    return r.stdout.strip() if r.returncode==0 else None

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign",type=Path,action="append")
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    repo=Path(__file__).resolve().parents[1]
    roots=[x.resolve() for x in args.campaign] if args.campaign else [repo/"output"/x for x in DEFAULTS]
    out=args.output or repo/"output"/("emi-precision-study-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+f"-{os.getpid()}")
    out.mkdir(parents=True,exist_ok=False)
    (out/"figures").mkdir()
    files={};allrows=[];steps=[];traces=[];operators=[];campaigns=[];cycles=[];four=[];titles={}
    def record(path):
        key=str(path.relative_to(repo)) if path.is_relative_to(repo) else str(path)
        if key not in files:files[key]=sha(path)
        return path
    for root in roots:
        print('Reading campaign:',root.name,flush=True)
        env=json.loads(record(root/"environment.json").read_text())
        summary=json.loads(record(root/"summary.json").read_text())
        case=env.get("case","10cells");duration=int(env["steps"])
        mixed=any('fgmres-r' in mode for mode in summary)
        titles[root.name]=f'{case} • {duration} step(s)'+(' • mixed' if mixed else (f" • cap {env.get('four_bit_limit',300)}" if duration==1 else ''))
        campaigns.append(dict(path=str(root),case=case,steps=duration,environment=env))
        base=root/case;op=read(record(base/"direct/potential-operator.csv"))
        a=np.array(op["a_one"]);at=np.array(op["a_transpose_one"]);diag=np.array(op["diagonal"])
        n=len(a)
        operators.append(dict(campaign=root.name,case=case,dofs=n,one_A_one=math.fsum(a),
            norm_A_one=float(np.linalg.norm(a)),constant_rayleigh=math.fsum(a)/n,
            mean_diagonal=float(diag.mean()),constant_rayleigh_over_mean_diagonal=math.fsum(a)/sum(diag),
            row_column_sum_discrepancy=float(np.linalg.norm(a-at))))
        refs={}
        for ref in ("direct","uncompressed-fgmres"):
            if summary[ref]["status"]!="PASS":raise ValueError("Missing completed reference")
            refs[ref]={s:read(record(base/ref/f"potential-state-step{s}.csv"))["u"] for s in range(duration+1)}
        for mode,row in summary.items():
            if mode=="direct":continue
            print('  ',mode,flush=True)
            folder=base/mode
            record(folder/"command.sh");record(folder/"run.log")
            m=row.get("metrics",{});success=row["status"]=="PASS"
            state0=read(record(folder/"potential-state-step0.csv"))["u"]
            if state0!=refs["direct"][0]:raise ValueError("Initial state differs")
            k=sum(m.get("iterations_per_step",[]))
            result=dict(campaign=root.name,case=case,steps=duration,mode=mode,status=row["status"],
                iterations=k,restriction_bits=row.get("restriction_bits"),
                prolongation_bits=row.get("prolongation_bits"),
                interface_bytes=row.get("interface_logical_bytes"),
                outer_bytes=row.get("outer_logical_bytes"),
                accounted_bytes=row.get("accounted_logical_bytes"),
                solver_seconds=m.get("solver_wrapper_max_seconds"),
                application_seconds=m.get("wall_seconds"),
                communication_seconds=None,measured_network_bytes=None)
            for s in range(1,duration+1):
                tracefile=folder/f"bddc-fgmres-step{s}.csv"
                if tracefile.exists():
                    with record(tracefile).open() as f:trace=list(csv.DictReader(f))
                    if [int(r["iteration"]) for r in trace]!=list(range(len(trace))):raise ValueError("Invalid trace order")
                    traces.append((root.name,case,duration,mode,s,trace))
                    if mode=="full-fgmres4-split":
                        item=dict(campaign=root.name,case=case,step=s,status=row['status'],
                            **trend(trace,max(0,int(trace[-1]["iteration"])-100)))
                        item["last10"]=trend(trace,max(0,int(trace[-1]["iteration"])-10))
                        abs_true=np.array([float(r["true_absolute_residual"]) for r in trace])
                        abs_est=np.array([float(r["estimated_absolute_residual"]) for r in trace])
                        item["max_estimate_true_relative_gap"]=float(np.max(np.abs(abs_true-abs_est)/np.maximum(abs_true,1e-300)))
                        four.append(item)
                        for begin in range(0,len(trace)-1,30):
                            chunk=trace[begin:min(begin+31,len(trace))]
                            cycles.append(dict(campaign=root.name,case=case,step=s,**trend(chunk)))
            if success:
                for s in range(1,duration+1):
                    u=(refs['uncompressed-fgmres'][s] if mode=='uncompressed-fgmres'
                       else read(record(folder/f"potential-state-step{s}.csv"))["u"])
                    item=dict(campaign=root.name,case=case,steps=duration,mode=mode,step=s,
                        iterations=m["iterations_per_step"][s-1],true_relative_residual=m["relative_residuals"][s-1])
                    for reference,values in refs.items():
                        err=error_metrics(u,values[s])
                        item.update({reference+"_"+key:v for key,v in err.items()})
                        # Algebraic coefficient-norm decomposition, not an FE norm.
                        e=np.asarray(u)-np.asarray(values[s])
                        if not np.isclose(np.dot(e,e),np.dot(e-e.mean(),e-e.mean())+len(e)*e.mean()**2,rtol=1e-11,atol=1e-28):
                            raise ValueError("Mean/error decomposition failed")
                    steps.append(item)
                result["final_true_relative_residual"]=m["relative_residuals"][-1]
                result.update({k:v for k,v in item.items() if k.startswith(("direct_","uncompressed-fgmres_"))})
            else:
                progress=row.get("residual_progress") or {}
                result["last_true_relative_residual"]=progress.get("last_true_relative")
                result["iterations"]=progress.get("last_iteration",k)
                if 'fgmres' in mode:
                    result['completed_application_outer_bytes_lower_bound']=7*result['iterations']*(2*n*8+12)
                    ledger=re.findall(r'^BDDC outer solver logical payload bytes: (\d+) collectives=(\d+)\s*$',
                                      (folder/'run.log').read_text(errors='replace'),re.M)
                    if len(ledger)==1:
                        payload,collectives=map(int,ledger[0]);k=result['iterations']
                        # Recover the single failed step's emitted ledger, not absent interface counters.
                        if duration==1 and payload==7*(k*(2*n*8+12)+4+56+n*8+16) and collectives==4*k+4:
                            result['outer_bytes']=payload;result['outer_collectives']=collectives
                            result['accounted_bytes_lower_bound']=payload
            for base_mode,short in (("uncompressed-fgmres","uncompressed"),("full-fgmres16","fgmres16"),("full-fgmres8-split","split8")):
                b=summary.get(base_mode,{})
                for name,field in (("interface","interface_logical_bytes"),("accounted","accounted_logical_bytes")):
                    result[f"{name}_saving_vs_{short}_percent"]=(100*(1-row[field]/b[field]) if success and b.get("status")=="PASS" and row.get(field) and b.get(field) else None)
            if success:
                result["mean_interface_bytes_per_application"]=result["interface_bytes"]/k
                baseline=summary["full-fgmres16"]
                kb=sum(baseline["metrics"]["iterations_per_step"])
                result["iteration_overhead_vs16_percent"]=100*(k/kb-1)
                result["break_even_iterations_vs16_at_observed_packet_mean"]=baseline["interface_logical_bytes"]/result["mean_interface_bytes_per_application"]
            allrows.append(result)
    # Verify the cap extension independently of the driver's summary flag.
    four_traces=[(c,t) for c,case,d,mode,s,t in traces if case=='10cells' and s==1 and mode=='full-fgmres4-split']
    prefix_checks=[]
    for old,old_trace in four_traces:
        for new,new_trace in four_traces:
            if len(old_trace)>=len(new_trace):continue
            identical=new_trace[:len(old_trace)]==old_trace
            if not identical:raise ValueError('4-bit prefix differs across cap experiments')
            prefix_checks.append(dict(original=old,extended=new,identical_rows=len(old_trace)))
    (out/'four-bit-prefix-checks.json').write_text(json.dumps(prefix_checks,indent=2)+'\n')
    # Hash implementation used for interpretation, not evidence it built every historical executable.
    for name in ("EMI_model.hh","EMI_model.cpp","EMI_fgmres.hh","EMI_fgmres_mpi.hh","tests/run_emi_fgmres_check.py"):
        record(repo/name)
    frozen=[]
    frozen_root=repo/'output/emi-frozen-probe-20261101T141220Z-13123'
    if (frozen_root/'summary.json').is_file():
        frozen=[dict(mode=m,**r) for m,r in json.loads(record(frozen_root/'summary.json').read_text()).items()]
        record(frozen_root/'environment.json')
        write_csv(out/'frozen-preconditioner.csv',frozen)
    for ancestor in repo.parents:
        library=ancestor/'kaskade7_new/mg/bddcSpaceTransfer.hpp'
        if library.is_file():
            record(library);record(library.with_suffix('.hh'));break
    record(Path(__file__).resolve())
    write_csv(out/"runs.csv",allrows);write_csv(out/"per-step.csv",steps)
    write_csv(out/"constant-mode.csv",operators);write_csv(out/"four-bit-cycles.csv",cycles)
    (out/"four-bit-analysis.json").write_text(json.dumps(four,indent=2,allow_nan=False)+"\n")
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),command=sys.argv,
        python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__,platform=platform.platform(),
        analysis_git_head=git_id(repo),campaigns=campaigns,input_sha256=files,
        note="Saved executable hashes remain in each campaign environment; current source hashes are not historical build provenance.")
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
        "axes.spines.right":False,"axes.grid":True,"grid.alpha":.2,"pdf.fonttype":42,
        "svg.fonttype":"none","savefig.dpi":300})
    figures=[]
    def save(fig,name,caption):
        for ext in ("png","pdf","svg"):fig.savefig(out/"figures"/f"{name}.{ext}",bbox_inches="tight")
        plt.close(fig);figures.append((name,caption))
    groups=list(dict.fromkeys((r["campaign"],r["case"],r["steps"]) for r in allrows))
    def panels():
        cols=min(3,len(groups));rows=math.ceil(len(groups)/cols)
        fig,axs=plt.subplots(rows,cols,figsize=(4.4*cols,4.5*rows),squeeze=False,layout="constrained")
        for ax in axs.flat[len(groups):]:ax.set_visible(False)
        return fig,list(axs.flat)[:len(groups)]
    metrics=[("iterations","FGMRES applications","iterations"),
        ("interface_bytes","BDDC interface [MB]","interface"),
        ("accounted_bytes","Accounted communication [MB]","accounted"),
        ("solver_seconds","Solver wrapper [s]","runtime")]
    for field,ylabel,name in metrics:
        fig,axs=panels()
        for ax,(campaign,case,duration) in zip(axs,groups):
            rr=[r for r in allrows if r["campaign"]==campaign and "fgmres" in r["mode"]]
            for i,r in enumerate(rr):
                v=r.get(field)
                lower=field=='accounted_bytes' and v is None and r.get('accounted_bytes_lower_bound') is not None
                if lower:v=r['accounted_bytes_lower_bound']
                if v is not None and (r["status"]=="PASS" or field=="iterations" or lower):
                    v=v/1e6 if "bytes" in field else v
                    ax.bar(i,v,color=COLORS.get(r["mode"],"#777777"),alpha=.85,
                           hatch="//" if r["status"]!="PASS" else None)
                    if field=='iterations' and r['status']=='PASS':
                        ax.annotate(str(r['iterations']),(i,v),ha='center',va='bottom',fontsize=8)
                    if r["status"]!="PASS":ax.annotate("lower bound; failed" if lower else "failed at cap",(i,v),ha="center",va="bottom",fontsize=8)
                else:ax.text(i,.05,"not measured",rotation=90,ha="center",transform=ax.get_xaxis_transform(),fontsize=8)
            ax.set_xticks(range(len(rr)),[LABELS.get(r["mode"],r["mode"]) for r in rr],rotation=50,ha="right")
            ax.set_title(titles[campaign]);ax.set_ylabel(ylabel)
        save(fig,name,f"{ylabel}. Separate campaigns and durations. Failed 4/4 is censored at its cap; missing measurements are not zero. Runtime includes validation and is a single run, not a speedup benchmark.")
    fig,axs=panels()
    for ax,(campaign,case,duration) in zip(axs,groups):
        rr=[r for r in allrows if r["campaign"]==campaign and "fgmres" in r["mode"]]
        for i,r in enumerate(rr):
            for off,key,marker in ((-.12,"relative_l2","o"),(.12,"centered_relative_l2","s")):
                v=r.get("uncompressed-fgmres_"+key)
                if v and v>0:ax.plot(i+off,v,marker,color=COLORS[r["mode"]],markersize=6)
            if r["mode"]=="uncompressed-fgmres":ax.text(i,.08,"self = 0",ha="center",fontsize=8,transform=ax.get_xaxis_transform())
            if r['status']!='PASS':ax.text(i,.08,'failed',ha='center',fontsize=8,transform=ax.get_xaxis_transform())
        ax.set_yscale("log");ax.set_ylim(1e-11,1e-3);ax.set_xlim(-.5,len(rr)-.5)
        ax.set_ylabel("Final relative coefficient error")
        ax.plot([],[],'ko',label='raw');ax.plot([],[],'ks',label='mean-centered');ax.legend(fontsize=8,loc='upper left')
        ax.set_xticks(range(len(rr)),[LABELS[r["mode"]] for r in rr],rotation=50,ha="right")
        ax.set_title(titles[campaign])
    save(fig,"accuracy","Final errors against uncompressed FGMRES: circles raw, squares mean-centered. Both use the uncentered reference norm. Direct comparisons and all per-step values are in CSV. Centering is diagnostic only.")
    fig,axs=panels()
    for ax,(campaign,case,duration) in zip(axs,groups):
        for r in allrows:
            if r["campaign"]!=campaign or r["status"]!="PASS" or "fgmres" not in r["mode"]:continue
            x=r["iteration_overhead_vs16_percent"];y=r["interface_saving_vs_fgmres16_percent"]
            dy=14 if r['mode']=='full-fgmres-r6-p8-split' else (-14 if r['mode']=='full-fgmres8-split' and 'mixed' in titles[campaign] else 4)
            ax.scatter(x,y,c=COLORS[r["mode"]]);ax.annotate(LABELS[r["mode"]],(x,y),xytext=(-4 if x>0 else 4,dy),
                ha='right' if x>0 else 'left',textcoords="offset points",fontsize=8)
        ax.axhline(0,color="gray",lw=.8);ax.axvline(0,color="gray",lw=.8)
        ax.margins(.2);ax.set_title(titles[campaign])
        ax.set_xlabel("Application overhead vs FGMRES16 [%]");ax.set_ylabel("Interface savings vs FGMRES16 [%]")
    save(fig,"tradeoff","Measured total interface saving versus iteration overhead. Lower-bit packets only help when the increased application count stays below break-even.")
    fig,axs=panels()
    for ax,(campaign,case,duration) in zip(axs,groups):
        for trace_campaign,c,d,mode,s,trace in traces:
            # match campaign as well as case and duration
            if trace_campaign!=campaign or s!=1:continue
            ax.semilogy([int(t["iteration"]) for t in trace],[float(t["true_relative_residual"]) for t in trace],
                        label=LABELS[mode],color=COLORS[mode])
        ax.axhline(1e-7,color="black",ls=":",lw=1)
        ax.set_title(titles[campaign]+'\nfirst-step history');ax.set_xlabel("FGMRES application")
        ax.set_ylabel("True relative residual");ax.legend(fontsize=8)
    save(fig,"residuals","True residual histories at the first step. Same plotted first step does not imply validation of all subsequent steps.")
    for campaign,case,duration,mode,s,trace in traces:
        if mode!="full-fgmres4-split":continue
        fig,axs=plt.subplots(1,2,figsize=(10,3.6),layout="constrained")
        for ax,start in zip(axs,(0,max(0,len(trace)-101))):
            tt=trace[start:]
            ax.semilogy([int(t["iteration"]) for t in tt],[float(t["true_relative_residual"]) for t in tt],color=COLORS[mode])
            ax.axhline(1e-7,color="black",ls=":",label="required tolerance")
            for k in range(start//30*30,len(trace),30):ax.axvline(k,color="gray",lw=.5,alpha=.4)
            ax.set_xlim(int(tt[0]['iteration']),int(tt[-1]['iteration']))
            ax.set_xlabel("FGMRES application");ax.set_ylabel("True relative residual")
            ax.set_title("Complete history" if start==0 else "Final 100 iteration intervals");ax.legend()
        save(fig,"four-bit-"+campaign,f"4/4, {case}, step {s}. Vertical lines mark restart boundaries. Short plateaus coexist with substantial reduction across cycles.")
    for campaign,case,duration in groups:
        if duration==1:continue
        fig,axs=plt.subplots(1,2,figsize=(10,3.6),layout="constrained")
        for mode in LABELS:
            rr=[r for r in steps if r["campaign"]==campaign and r["mode"]==mode]
            if not rr:continue
            axs[0].plot([r["step"] for r in rr],[r["iterations"] for r in rr],"o-",color=COLORS[mode],label=LABELS[mode])
            axs[1].semilogy([r["step"] for r in rr],[r["true_relative_residual"] for r in rr],"o-",color=COLORS[mode])
        axs[0].set_ylabel("Applications");axs[0].legend(fontsize=8)
        axs[1].set_ylabel("True relative residual");axs[1].axhline(1e-7,color="black",ls=":")
        for ax in axs:ax.set_xlabel("Time step")
        fig.suptitle(case+" • ten-step validation")
        save(fig,"steps-"+campaign,"Iterations and accepted residuals at every time step; numerical tolerances remain unchanged.")
        modes=list(dict.fromkeys(mode for c,_,_,mode,_,_ in traces if c==campaign))
        fig,axs=plt.subplots(1,len(modes),figsize=(3.5*len(modes),3.7),squeeze=False,layout='constrained')
        norm=matplotlib.colors.Normalize(1,duration);cmap=plt.get_cmap('viridis')
        for ax,mode in zip(axs[0],modes):
            for c,_,_,m,s,trace in traces:
                if c==campaign and m==mode:
                    ax.semilogy([int(t['iteration']) for t in trace],[float(t['true_relative_residual']) for t in trace],color=cmap(norm(s)),alpha=.8)
            ax.axhline(1e-7,color='black',ls=':',lw=1);ax.set_title(LABELS[mode])
            ax.set_xlabel('FGMRES application');ax.set_ylabel('True relative residual')
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axs[0],label='Time step',ticks=range(1,duration+1))
        fig.suptitle(case+' • all ten residual histories')
        save(fig,'all-histories-'+campaign,'Every recorded FGMRES residual history, colored by time step; no averaging of convergence trajectories.')
    sections=[]
    four_pass=[r for r in allrows if r['mode']=='full-fgmres4-split' and r['status']=='PASS']
    six_ten=[r for r in allrows if r['mode']=='full-fgmres6-split' and r['steps']==10 and r['status']=='PASS']
    sections.append("<h1>EMI precision, convergence and communication</h1><p>Reproducible analysis of supplied experiments. No numerical implementation was changed by this analysis. One-step evidence, ten-step validation, logical interface payload and modeled outer traffic are kept separate.</p>")
    if four_pass:
        r=four_pass[-1]
        sections.append(f"<p><strong>4/4 converges in {r['iterations']} applications</strong> on the one-step 10cells test, with true relative residual {r['final_true_relative_residual']:.7e}. The original 300-iteration prefix is checked independently against the extended run. The driver skipped 1,000 because convergence was already reached. This is one-step evidence, not ten-step validation.</p>")
    if six_ten:
        r=six_ten[-1]
        sections.append(f"<p><strong>6/6 passes all ten 10cells steps</strong>: {r['iterations']} applications, {r['interface_bytes']:,} interface bytes and {r['accounted_bytes']:,} accounted bytes. It uses {-r['interface_saving_vs_split8_percent']:.2f}% more interface bytes than the matched 8/8 control. Eight bits remains the best validated uniform interface choice over ten steps; the mixed one-step screen below is a separate result.</p>")
    sections.append("<h2>4/4-bit residual evidence</h2>"+table([{k:v for k,v in r.items() if k!="last10"} for r in four],
        ["campaign","status","first_iteration","last_iteration","first_residual","last_residual","end_over_start","increases","max_estimate_true_relative_gap"]))
    sections.append("<p>The original 300-iteration trace is not globally stagnating: iterations 200–300 decrease by about 75.7-fold without an increase. The last ten intervals are nearly flat, but preceding restart cycles show renewed reduction. Estimated and explicitly recomputed residuals agree closely, arguing against a misleading residual estimate. Restart-related plateaus are observed; quantization-degraded correction directions are a plausible mechanism, not separately established by these residuals.</p>")
    if four_pass:
        sections.append('<p>The cap extension crosses the tolerance on iteration 301, the first iteration of cycle 11: 1.417261494e-7 → 5.778563937e-8. This is renewed progress at a restart, not a demonstrated accuracy floor at 1.417e-7. The longer-cap run reproduces all 301 original observations (iterations 0–300). No tolerance, precision, pass-count or numerical algorithm change was needed.</p>')
    sections.append("<h2>Measured campaign totals</h2>"+table(allrows,["case","steps","mode","status","iterations","interface_bytes","outer_bytes","accounted_bytes","solver_seconds","final_true_relative_residual"]))
    sections.append('<p><strong>Experimental design:</strong> physics, mesh, tolerance and outer solver are matched within each FGMRES campaign. The established 16/16 control uses original restriction/prolongation quantizers without constant splitting. The 8/8, 6/6 and 4/4 candidates use exact-zero restriction and constant-split affine prolongation. Thus comparisons among the three low-bit modes isolate bit depth, while comparison with 16/16 evaluates a composite method change. A 16/16 control with the same zero/split policies is needed to isolate precision alone across that boundary. CG16 is an additional solver baseline, not an FGMRES precision data point.</p>')
    sections.append('<p>The failed 300-application 4/4 run emitted an outer ledger before aborting: 1,240,585,892 logical bytes in 1,204 collectives. This matches the adapter formula for 8 ranks and 36,860 DOFs. Its interface counter is missing, so that outer ledger is a lower bound on accounted traffic, not a full total or measured network traffic. It already exceeds successful FGMRES16 accounted traffic by over 13 times; extending the cap cannot make this run optimal under that accounting scope.</p>')
    sections.append("<p>MB in plots means 10^6 bytes. Interface bytes include metadata, padding, exact trailers and intra-rank logical deliveries. Accounted bytes add the existing outer-collective payload model: Bcast/Reduce=(P−1)×payload and Allreduce=2(P−1)×payload. This excludes coarse-grid traffic, setup, training, internal BDDC reductions and MPI framing. No hardware network measurement or isolated communication-time measurement is available. Those fields are null, not zero. Root owns Arnoldi and matrix products; these centralized experiments are not a distributed scalability benchmark.</p>")
    sections.append("<h2>Accuracy and the mean offset</h2>"+table(allrows,["case","steps","mode","direct_relative_l2","direct_centered_relative_l2","uncompressed-fgmres_relative_l2","uncompressed-fgmres_centered_relative_l2","direct_mean_fraction_of_squared_error"]))
    sections.append("<p>Errors are Euclidean norms of unique algebraic potential coefficients, not volume-weighted finite-element L2 errors. Centered errors remove only the arithmetic mean of the error, keeping the original reference norm as denominator. Raw errors are retained. Each time step is compared with the corresponding reference trajectory, so later errors include prior-step effects. No physical accuracy threshold was imposed.</p>")
    sections.append("<h3>The constant mode is weakly anchored, not a demonstrated gauge freedom</h3>"+table(operators,["case","dofs","one_A_one","norm_A_one","constant_rayleigh","mean_diagonal","constant_rayleigh_over_mean_diagonal"]))
    sections.append("<p>EMI_model.hh DomainCache uses bulk gradients; InnerBoundaryCache uses differences of potentials, which cancel a common constant. BoundaryCache adds an extracellular exterior term with gamma=10^-5 referenced to zero. The saved operator has nonzero A·1 and positive 1ᵀA1. Thus the common constant is not an exact null vector in these runs. With the usual connected, positively coupled regions this boundary anchor removes the common constant ambiguity; these diagnostics alone are not a proof of full matrix rank or a condition-number measurement.</p><p>For error e=c·1+e₀ with arithmetic mean(e₀)=0, ||e||²=Nc²+||e₀||², while A(c·1)=cA1. A small A1 allows appreciable mean error with a small residual. The coefficient-error decomposition was checked numerically. This explains the observed sensitivity without assuming the offset is physically irrelevant. A centered potential error is not by itself a membrane-voltage, current or gradient error certificate.</p>")
    sections.append("<h2>Why the split and flexible solver help</h2><p>The split transmits a float32 producer constant c separately and quantizes w(u−c), reconstructing Q[w(u−c)]+wc. It protects the weighted constant component from the low-bit grid; it does not enforce the global reference potential or remove a final solution offset. Scale and minimum remain float32, with an additional 8-byte header (float32 c and uint32 tag) already counted. Exact corner values remain exceptions to the quantized payload: 8/8 refers to the quantized values, not every transmitted field.</p><p>Quantization is input dependent and nonlinear; varying scales, rounding and any feedback can change the effective preconditioner. Lossless integer lifting/Huffman changes representation, not the decoded quantization error. FGMRES stores each actual preconditioned direction z_j and minimizes the residual over their span, accommodating varying/nonlinear preconditioning without requiring the CG positivity recurrence. It does not guarantee fast convergence when Az_j directions are poor. See <a href='https://epubs.siam.org/doi/abs/10.1137/0914028'>Saad (1993)</a> and <a href='https://petsc.org/release/manualpages/KSP/KSPFGMRES/'>PETSc FGMRES documentation</a>.</p>")
    sections.append("<h2>Which precision currently minimizes communication?</h2><p>Among the measured successful uniform one-step FGMRES configurations, 8/8 split minimizes interface traffic. The new mixed R6/P8 screen reduces this further, but has only one-step validation. 6/6 uses about 23.5% fewer interface bytes per application than 8/8 but needs 31 rather than 22 applications, so its solve communicates 7.83% more interface bytes. Its break-even count at the observed mean packet cost is about 28.7 applications. This is an accounting estimate, not a prediction for a changed trajectory. 16/16 minimizes the modeled interface-plus-outer total among these successful FGMRES runs. The optimum depends on the communication scope; neither a global optimum nor a hardware optimum has been established.</p><p>In the three available ten-step studies, 8/8 split saves approximately 84.92%, 77.36%, and 72.00% interface bytes versus uncompressed FGMRES for 10cells, weak_8_h16, and strong_h16, respectively. It uses more accounted bytes than FGMRES16 in all three. New 6/6 ten-step evidence covers 10cells only; no corresponding weak/strong 6/6 claim is made. The original failed 4/4 campaign lacks interface totals, while its successful extended-cap campaign reports complete totals.</p>")
    if four_pass:
        sections.append('<p>At 301 applications, 4/4 sends 9,570,528 interface bytes: 10.17% more than uncompressed FGMRES and 625.88% more than split8. Its smaller per-application representation is overwhelmed by the extra iterations. Demonstrating four-bit convergence does not make four bits an efficient operating point.</p>')
    if six_ten:
        sections.append('<p>Ten-step 6/6 requires [31,32,32,32,32,33,32,33,32,33] applications, all below the unchanged per-step cap. Its maximum accepted relative residual is 9.81994e-8. Versus split8, per-application interface bytes fall by about 23.7%, but applications rise by 41.9%, yielding 8.25% more total interface bytes and 40.49% more accounted bytes. The final direct-reference raw/centered errors are 1.898481e-5 and 1.728000e-9; the raw offset still dominates. The observed wrapper time is 63.878322 s versus 47.185566 s for split8, a single validation-run comparison rather than a timing benchmark.</p>')
    sections.append("<h2>Measurement gaps and controlled next experiments</h2><p>The current FGMRES traces do not contain per-transfer quantization/reconstruction errors, so no precision-wise reconstruction error curve is fabricated. Earlier same-input frozen probes measured whole preconditioner z and Az errors, not scalar codec errors: 8/8 without splitting had relative Az error ≈421, versus ≈0.017 with splitting, on that one input. This supports the split mechanism but does not diagnose late 4-bit iterations. Isolated MPI/codec time and actual network traffic are also unavailable. Existing wrapper timings include validation/output and have no replicate uncertainty estimates.</p><p>The completed one-step R8/P6 and R6/P8 screen changes one directional precision parameter at a time from a fresh 8/8 split control. R8/P6 needs 32 applications and 1,708,232 interface bytes, 29.56% more than 8/8. R6/P8 keeps 22 applications and sends 1,150,088 interface bytes: 12.77% less than 8/8 and 86.76% less than uncompressed FGMRES. It is the best measured one-step interface configuration in this study. This supports greater sensitivity to prolongation precision in this particular experiment, not a universal directional ordering. R6/P8 raw/centered direct-reference errors are 3.996068e-4 and 3.748169e-9, versus 5.440595e-5 and 1.115847e-9 for 8/8; the larger errors require continued reporting and scientific assessment. Its accounted total is 94,039,668 bytes, only 0.18% below 8/8 and still 3.13% above FGMRES16. At 22 applications, the outer payload alone (92,889,580 bytes) exceeds FGMRES16 total accounted bytes (91,186,880). The next controlled experiment validates R6/P8 over ten 10cells steps with fresh controls and exact first-step parity against its successful screen. Restart=30, cap=300, float32 metadata, one restriction pass, no fallback, physics and tolerance remain fixed. R8/P4 and R4/P8 are deferred until this validation. Neither ten-step R6/P8 reliability nor a global precision optimum has been established.</p>")
    if frozen:
        sections.append('<h3>Separate frozen-input evidence (one application only)</h3>'+table(frozen,
            ['mode','status','z_relative_error','Az_relative_error','cosine_v_Az','interface_bytes']))
    sections.append("<h2>Reproducibility</h2><p>manifest.json contains input SHA-256 hashes, campaign settings, executable hashes recorded by the original campaigns, software versions and the analysis command. Source hashes recorded now document interpretation, not proof of the source used to build historical binaries. runs.csv, per-step.csv, constant-mode.csv and four-bit-cycles.csv contain full numerical tables. Figures are supplied in PDF, SVG and 300-dpi PNG. Previous output directories are never overwritten. This report is generated from the listed inputs and must be regenerated after new experiments.</p>")
    for name,caption in figures:
        sections.append(f"<figure><img src='figures/{name}.png'><figcaption>{html.escape(caption)} <a href='figures/{name}.pdf'>PDF</a> · <a href='figures/{name}.svg'>SVG</a></figcaption></figure>")
    (out/"report.html").write_text("<!doctype html><html lang='en'><meta charset='utf-8'><title>EMI precision study</title><style>body{max-width:1200px;margin:3em auto;padding:0 1em;font:16px/1.6 system-ui;color:#18212b}table{display:block;overflow:auto;border-collapse:collapse;font:12px/1.5 monospace}td,th{padding:6px;border-bottom:1px solid #ddd;text-align:right;white-space:nowrap}h1,h2,h3{line-height:1.2}figure{margin:3em 0}img{width:100%}figcaption{font-size:14px;color:#43505c}a{color:#00659d}</style>"+"".join(sections)+"</html>")
    print(out)
    print(json.dumps(four,indent=2))
    return 0

if __name__=="__main__":raise SystemExit(main())
