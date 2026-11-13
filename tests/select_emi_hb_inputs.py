#!/usr/bin/env python3
"""Reuse a complete, balanced prefix of frozen FCG inputs without recapturing."""
import argparse,hashlib,json,re
from pathlib import Path

def select(source,output,calls):
    source=source.resolve();groups={}
    for p in source.glob('*.buffer'):
        m=re.fullmatch(r'(restrict|prolongate)-(\d+)-sub(\d+)\.buffer',p.name)
        if not m or int(m[2])>calls:continue
        phase,call,producer=m[1],int(m[2]),int(m[3])
        with p.open() as f:header=f.readline().split()
        if header[:4]!=['KASKADE_FROZEN_V1',phase,str(call),str(producer)]:raise ValueError('Header mismatch: '+str(p))
        groups.setdefault((phase,call),{})[producer]=p
    expected={(phase,call) for phase in ['restrict','prolongate'] for call in range(1,calls+1)}
    if set(groups)!=expected:raise ValueError('Missing sampled phase/call')
    ids=set(next(iter(groups.values())))
    if not ids or any(set(g)!=ids for g in groups.values()):raise ValueError('Incomplete producer group')
    output.mkdir(parents=True,exist_ok=False);dest=output/'inputs';dest.mkdir()
    hashes={}
    for group in groups.values():
        for p in group.values():
            (dest/p.name).symlink_to(p)
            hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    (output/'sample-provenance.json').write_text(json.dumps({'source':str(source),'calls_per_phase':calls,'producers':sorted(ids),'files_sha256':hashes,'scope':'Early-call sample of existing 32-bit quantization-only FCG trajectory, split OFF; input files unchanged.'},indent=2)+'\n')
    print(f'Selected {len(hashes)} buffers: {calls} calls per phase, {len(ids)} producers.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--calls',type=int,default=3)
    a=p.parse_args()
    if a.calls<1 or a.calls>32:p.error('calls must be 1..32')
    select(a.source,a.output,a.calls)
