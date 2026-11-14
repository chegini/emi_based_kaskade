#!/usr/bin/env python3
"""Check full production R/P reconstruction and original uint32 codes on fixed inputs.
Run after compiling test_bddc_split_recenter.cpp inside the production container.
Raw audit files are removed after verification, keeping only bounded summaries.
"""
import argparse,csv,hashlib,os,shutil,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--executable',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--bits',nargs='+',type=int,default=[8,16,32]);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
rows=[]
for bits in a.bits:
 for depth in [1,2,3,4]:
  for split in [0,1]:
   for kind in range(5):
    name=f'b{bits}-d{depth}-s{split}-k{kind}'
    tmp=Path(tempfile.mkdtemp(prefix=name+'-',dir=a.output))
    try:
     states=[];packets=[];qhashes=[];grids=[];ancillary=[]
     for recenter in [0,1]:
      out=tmp/str(recenter);audit=out/'audit'
      env=dict(os.environ,KASKADE_TRANSFER_AUDIT_DIR=str(audit),KASKADE_LIFTING_PROFILE_PREFIX='',KASKADE_TRANSFER_TIMING='0')
      cmd=[str(a.executable),str(bits),str(depth),str(split),str(recenter),str(kind),str(out)]
      with (tmp/f'run-{recenter}.log').open('w') as log:
       subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
      states.append((out/'state.bin').read_bytes());packets.append((out/'packets.bin').read_bytes())
      codes={};checked=0
      for send in audit.glob('*_send_*_coefficients.txt'):
       phase,_,producer,consumer,_=send.stem.split('_')
       payload=send.read_bytes();codes[send.name]=hashlib.sha256(payload).hexdigest()
       recv=audit/f'{phase}_recv_{consumer}_{producer}_coefficients.txt'
       if recv.is_file():
        assert recv.read_bytes()==payload,(name,recenter,send,'Quantized codes changed in production receiver')
        checked+=1
      if kind not in [2,3]:assert checked==12,(name,'Expected all six directed interfaces for both R/P')
      qhashes.append(codes)
      grids.append({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in audit.glob('*_send_*_metadata.txt')})
      exact=67 if kind==2 else 3 if kind==1 else 0
      unchanged={}
      for f in audit.glob('*_send_*_wire.txt'):
       lines=f.read_text().splitlines();wire=bytes(map(int,lines[1:]))
       assert len(wire)==int(lines[0]),'Incomplete wire audit'
       if split and f.name.startswith('prolongate_'):unchanged[f.name+'-constant']=wire[:8]
       if exact:unchanged[f.name+'-corners']=wire[-8*exact:]
      ancillary.append(unchanged)
     assert states[0]==states[1],(name,'Recentering changed reconstructed values')
     assert qhashes[0]==qhashes[1],(name,'Recentering changed original quantized coefficients')
     assert grids[0]==grids[1],(name,'Recentering changed the wire quantizer grid')
     assert ancillary[0]==ancillary[1],(name,'Recentering changed the constant/corner bytes')
     if kind in [2,3]:assert packets[0]==packets[1],(name,'Raw fallback packets changed')
     rows.append(dict(bits=bits,depth=depth,split=split,input_kind=kind,state_sha256=hashlib.sha256(states[0]).hexdigest(),quantized_code_streams=len(qhashes[0]),status='PASS'))
    except Exception:
     print(f'FAIL: retained traces in {tmp}',flush=True)
     raise
    else:
     shutil.rmtree(tmp)
    print('PASS '+name,flush=True)
with (a.output/'summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(f'PASS: {len(rows)} fixed-input split/depth cases, each with recentering OFF/ON; exact quantized codes, reconstruction, constants/corners and compact packets.')
