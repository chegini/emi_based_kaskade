#!/usr/bin/env python3
"""Controlled 4-bit cap extension followed by 6-bit ten-step validation.

Only the existing diagnostic driver's arguments change. No numerical code changes.
Run inside the Kaskade build container. --dry-run prints commands without solving.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys

CONTROLS=("direct","full-cg16","uncompressed-fgmres","full-fgmres16","full-fgmres8-split")
MODE="full-fgmres4-split"

def can_extend(root, cap):
    """Conservative evidence gate: a cap exit, parity, controls, continuing decline."""
    records=json.loads((root/"summary.json").read_text())
    if any(records.get(m,{}).get("status")!="PASS" for m in CONTROLS):
        return False,"A fresh control failed."
    row=records.get(MODE,{})
    if row.get("status")=="PASS":return False,"4/4 already converged."
    if not row.get("saved_trace_prefix_rows"):return False,"Saved residual prefix parity is missing."
    folder=root/"10cells"/MODE
    text=(folder/"run.log").read_text(errors="replace").lower()
    if "iteration limit" not in text and "iteration-limit" not in text:
        return False,"Failure was not identified as an iteration-limit exit."
    with (folder/"bddc-fgmres-step1.csv").open() as f:rows=list(csv.DictReader(f))
    if [int(r["iteration"]) for r in rows]!=list(range(cap+1)):
        return False,"Trace does not reach the requested cap."
    values=[float(r["true_relative_residual"]) for r in rows[-101:]]
    if any(not math.isfinite(v) or v<=0 for v in values):return False,"Invalid residual."
    factor=values[-1]/values[0]
    # This is a declared scheduling criterion, not a solver stopping criterion.
    decreasing=factor<1-1e-8 and values[-1]<=min(values)*(1+1e-8)
    return decreasing,f"Final 100 intervals: end/start={factor:.12g}; extend only if decreasing beyond relative 1e-8 noise and last is best."

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference-campaign",type=Path,required=True,
                   help="successful ten-step 10cells split8 campaign")
    p.add_argument("--screening-campaign",type=Path,required=True,
                   help="one-step campaign containing successful 6/6 and failed 4/4 at 300")
    p.add_argument("--output",type=Path)
    p.add_argument("--dry-run",action="store_true")
    a=p.parse_args();repo=Path(__file__).resolve().parents[1]
    ref=a.reference_campaign.resolve();screen=a.screening_campaign.resolve()
    prior=json.loads((screen/"summary.json").read_text())
    if prior.get("full-fgmres6-split",{}).get("status")!="PASS":p.error("6/6 screening must pass")
    progress=prior.get(MODE,{}).get("residual_progress") or {}
    if progress.get("last_iteration")!=300 or progress.get("last_true_relative",0)<=1e-7:
        p.error("Expected failed 300-iteration 4/4 screening")
    root=(a.output or repo/"output"/("emi-precision-followup-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+f"-{os.getpid()}")).resolve()
    root.mkdir(parents=True,exist_ok=False)
    driver=repo/"tests/run_emi_fgmres_check.py"
    def command(bits,steps,folder,limit=300,previous=screen):
        cmd=[sys.executable,str(driver),"--reference-campaign",str(ref),
             "--screening-campaign",str(previous),"--constant-split","--case","10cells",
             "--bits",str(bits),"--steps",str(steps),"--output",str(folder)]
        if limit!=300:cmd+=["--four-bit-limit",str(limit)]
        return cmd
    commands=[command(4,1,root/"b4-cap500",500),
              command(4,1,root/"b4-cap1000",1000,root/"b4-cap500"),
              command(6,10,root/"b6-ten-steps")]
    plan=dict(reference_campaign=str(ref),screening_campaign=str(screen),
        commands=commands,dry_run=a.dry_run,tolerance=1e-7,restriction_passes=1,
        restart=30,metadata="float32 scale/minimum and split center; 8-byte split header",
        cap1000_gate="Controls PASS; saved prefix identical; cap exit; last-100 end/start<1-1e-8; last residual is best within relative 1e-8",
        note="Fresh starts at each cap, not checkpoint continuation. Controls retain cap 300. Six-bit uses cap 300 per step.")
    (root/"plan.json").write_text(json.dumps(plan,indent=2)+"\n")
    (root/"commands.txt").write_text("\n".join(shlex.join(c) for c in commands)+"\n")
    print("Follow-up campaign:",root,flush=True)
    if a.dry_run:
        print("NOT EXECUTED. Second command is conditional:\n"+(root/"commands.txt").read_text())
        return 0
    rc500=subprocess.run(commands[0],cwd=repo).returncode
    try:
        allow,reason=can_extend(root/"b4-cap500",500)
    except (OSError,ValueError,KeyError) as e:
        allow,reason=False,'Incomplete 500-cap evidence: '+str(e)
    decisions=dict(cap500_exit=rc500,cap1000_allowed=allow,reason=reason)
    print("1000-cap decision:",reason,flush=True)
    if allow:decisions["cap1000_exit"]=subprocess.run(commands[1],cwd=repo).returncode
    (root/"decisions.json").write_text(json.dumps(decisions,indent=2)+"\n")
    # Independent six-bit validation still runs when four bits fails to converge.
    decisions["six_bit_ten_step_exit"]=subprocess.run(commands[2],cwd=repo).returncode
    (root/"decisions.json").write_text(json.dumps(decisions,indent=2)+"\n")
    print("Results:",root)
    four_root=root/("b4-cap1000" if allow else "b4-cap500")
    four_summary=four_root/"summary.json"
    four_status=(json.loads(four_summary.read_text()).get(MODE,{}).get("status") if four_summary.is_file() else None)
    return int(four_status!="PASS" or decisions["six_bit_ten_step_exit"]!=0)

if __name__=="__main__":raise SystemExit(main())
