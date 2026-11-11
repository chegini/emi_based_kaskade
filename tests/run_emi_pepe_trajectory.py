#!/usr/bin/env python3
"""Matched raw/16-bit split FCG trajectories; no assumed error acceptance threshold."""
import csv
import json
import os
from pathlib import Path
import shlex
from run_emi_solver_study import arguments, make_command, execute
from run_emi_mpi_routing_check import digest
from collect_emi_scaling_smoke import collect
from analyze_emi_potential_diagnostics import read, error_metrics
from run_emi_histogram_start_check import require


def main():
    a = arguments()
    require(a.case == 'pepe' and a.solvers == ['fcg'] and a.bits == [16]
            and a.split == 'on', 'Require Pepe, FCG, 16 bits, split on')
    a.steps = 100
    root = a.output.resolve()
    if not a.plan_only:
        require(int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) == a.expected_nodes,
                'Allocation node count mismatch')
        root.mkdir(parents=True, exist_ok=False)
    cases = [dict(name='fcg-uncompressed', solver='fcg', bits=None, split='off', direct=False),
             dict(name='fcg-before-16-split-on', solver='fcg', bits=16, split='on', direct=False)]
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    for key in ('KASKADE_RESTRICTION_CAPTURE', 'KASKADE_FROZEN_TRANSFER_DIR',
                'KASKADE_TRANSFER_AUDIT_DIR', 'EMI_MATERIAL_TRACE'):
        env.pop(key, None)
    report = dict(complete=False, runs={}, note='Compression-induced global coefficient errors; no scientific threshold. dt=0.01; final time=1.0 in model units. VTK and snapshots ON.')
    if not a.plan_only:
        model_hash = digest(a.executable or a.repo/'emiModel')
        (root/'environment.json').write_text(json.dumps(dict(configuration=vars(a),
            executable_sha256=model_hash, nodes=os.environ.get('SLURM_JOB_NODELIST')),default=str,indent=2))
    for case in cases:
        folder = root/case['name']
        cmd = make_command(a, case, folder, ['mpirun'])
        cmd[cmd.index('--vtk')+1] = '1'
        if a.plan_only:
            print(shlex.join(cmd)); continue
        folder.mkdir()
        require(digest(a.executable or a.repo/'emiModel') == model_hash, 'Executable changed')
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print('Running '+case['name']+'; 100 steps, VTK and full-precision snapshots ON',flush=True)
        rc, timed_out, elapsed = execute(cmd, folder, a, env)
        row = collect((folder/'run.log').read_text(errors='replace'),100,a.ranks,rc,1e-7)
        row.update(timed_out=timed_out, elapsed_seconds=elapsed)
        report['runs'][case['name']] = row
        (root/'summary.json').write_text(json.dumps(report,indent=2))
        require(not timed_out and not row['errors'], str(row['errors']))
        print('PASS: all 100 true-residual checks',flush=True)
    if a.plan_only: return
    raw, compressed = [root/c['name'] for c in cases]
    require(read(raw/'potential-state-step0.csv')['u'] == read(compressed/'potential-state-step0.csv')['u'], 'Initial states differ')
    require(read(raw/'potential-state-step1.csv')['rhs'] == read(compressed/'potential-state-step1.csv')['rhs'], 'Initial physical RHS differs')
    rows=[]
    for step in range(101):
        reference = read(raw/f'potential-state-step{step}.csv')['u']
        current = read(compressed/f'potential-state-step{step}.csv')['u']
        row = dict(step=step,time=step*.01,**error_metrics(current,reference))
        for case in cases:
            record = report['runs'][case['name']]
            row[case['name']+'_iterations'] = record['iterations_per_step'][step-1] if step else 0
            row[case['name']+'_residual'] = record['relative_residuals'][step-1] if step else None
        rows.append(row)
    with (root/'compression-error.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (root/'difference-step100.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['dof','uncompressed','compressed','difference'])
        w.writerows((i,b,c,c-b) for i,(b,c) in enumerate(zip(reference,current)))
    report.update(complete=True, final_error=rows[-1],
                  peak_absolute_l2=max(rows,key=lambda x:x['absolute_l2']))
    (root/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print('Final error: '+json.dumps(rows[-1]),flush=True)
    print('Results: '+str(root),flush=True)

if __name__ == '__main__':
    main()
