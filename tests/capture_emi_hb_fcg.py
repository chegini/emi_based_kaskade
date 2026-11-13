#!/usr/bin/env python3
"""Derive a 32-bit quantization-only FCG capture from the raw command with read-only frozen-input capture.
Runs under an eight-node Slurm allocation. Never overwrites the prior results.
"""
import argparse,hashlib,json,shlex,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-run',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
cmd=shlex.split((a.source_run/'command.sh').read_text())
i=next(i for i,x in enumerate(cmd) if x.endswith('/build/emiModel'))
model=Path(cmd[i]);assert cmd[cmd.index('--bddcFlexibleCg')+1]=='1'
assert cmd[cmd.index('--bddcCompression')+1]=='0'
assert cmd[cmd.index('--bddcProlongateConstantSplit')+1]=='0'
# The raw transfer class has no capture observer. Select the compressed class
# with 32-bit nodal quantization, transforms/encoding/splitting still disabled.
cmd[cmd.index('--bddcCompression')+1]='1'
for flag in ['--bddcCompressionBits','--bddcRestrictCompressionBits','--bddcProlongateCompressionBits']:
    cmd[cmd.index(flag)+1]='32'
cmd[cmd.index('--dir')+1]=str(a.output/'solve');(a.output/'solve').mkdir()
cmd[i:i]=['env','KASKADE_FROZEN_TRANSFER_DIR='+str(a.output/'inputs')]
(a.output/'capture-command.sh').write_text(shlex.join(cmd)+'\n')
(a.output/'provenance.json').write_text(json.dumps({'source_run':str(a.source_run.resolve()),'executable':str(model),'executable_sha256':hashlib.sha256(model.read_bytes()).hexdigest(),'trajectory':'32-bit quantization-only FCG, split OFF; first 32 calls per phase; not the split-ON compressed trajectory'},indent=2)+'\n')
with (a.output/'capture.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
assert list((a.output/'inputs').glob('*.buffer')),'No inputs captured'
