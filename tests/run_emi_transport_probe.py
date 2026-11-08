#!/usr/bin/env python3
"""Bounded Slurm/Apptainer transport checks before running Pepe; no EMI rebuild."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shlex
import subprocess

from run_emi_mpi_routing_check import arguments, digest, slurm_container_command


def stage_arguments(repo, image, kaskade, nodes, ranks, btl):
    return arguments(['--repo',str(repo),'--slurm-container',str(image),
                      '--kaskade',str(kaskade),'--expected-nodes',str(nodes),
                      '--ranks',str(ranks),'--slurm-btl',btl,
                      '--slurm-vader-single-copy','none'])


def bounded(command):
    return ['--time=00:03:00' if x=='--time=00:10:00' else x for x in command]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--slurm-container',type=Path,required=True)
    p.add_argument('--kaskade',type=Path,required=True)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    a.repo=a.repo.resolve(); a.slurm_container=a.slurm_container.resolve(); a.kaskade=a.kaskade.resolve()
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))<2:
        p.error('Run inside an allocation with at least two nodes and two tasks per node')
    if not a.slurm_container.is_file(): p.error('Container image is missing')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=(a.output or a.repo/'output'/f'emi-transport-probe-{stamp}-{os.getpid()}').resolve()
    if not root.is_relative_to(a.repo.parent): p.error('Output must be beneath the bound BDDC parent')
    root.mkdir(parents=True,exist_ok=False)
    source=a.repo/'tests/emi_mpi_transport_probe.cpp'; binary=root/'transport-probe'
    one=stage_arguments(a.repo,a.slurm_container,a.kaskade,1,2,'self,vader')
    two=stage_arguments(a.repo,a.slurm_container,a.kaskade,2,4,'self,vader,tcp')
    report=dict(source_sha256=digest(source),job=os.environ['SLURM_JOB_ID'],runs={},
                note='Transport integrity only; successful tiny checks do not validate Pepe or establish performance.')
    def run(name,cmd,marker=None):
        cmd=bounded(cmd)
        (root/(name+'.command.sh')).write_text(shlex.join(cmd)+'\n')
        print('Running '+name+' ...',flush=True)
        with (root/(name+'.log')).open('w') as stream:
            rc=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT).returncode
        log=(root/(name+'.log')).read_text(errors='replace')
        ok=rc==0 and (marker is None or marker in log)
        report['runs'][name]=dict(exit_code=rc,status='PASS' if ok else 'FAIL',command=cmd)
        (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
        print('  '+('PASS' if ok else 'FAIL; see '+str(root/(name+'.log'))),flush=True)
        return ok
    print('Output: '+str(root),flush=True)
    build=slurm_container_command(one,['/bin/bash','-c',
        'set -eu; ompi_info --version; ompi_info --param btl vader --level 9; '
        'mpicxx -std=c++17 -O2 -Wall -Wextra "$1" -o "$2"',
        'build',str(source),str(binary)],1,1,mpi='none')
    if not run('build-and-mpi-info',build): return 1
    report['binary_sha256']=digest(binary)
    for name,config,nodes,ranks in [('one-node-no-tcp',one,1,2),('two-nodes',two,2,4)]:
        cmd=slurm_container_command(config,[str(binary),str(ranks),str(nodes)],nodes,ranks,mpi='pmix_v5')
        # Verbose BTL selection is diagnostic-only and never used for timings.
        script=cmd.index('-c')+1
        cmd[script]=cmd[script].replace('exec "$@"',
            'export OMPI_MCA_btl_base_verbose=30\nexec "$@"')
        if not run(name,cmd,f'PASS: transport payload integrity; ranks={ranks} nodes={nodes}'): return 1
    print('Transport probes passed. Next run the bounded Pepe setup diagnostic; no performance claim.',flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
