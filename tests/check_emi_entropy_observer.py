#!/usr/bin/env python3
"""Full production-template parity of original vs observer, run in Dune container."""
import argparse,hashlib,itertools,json,os,re,shlex,subprocess
from pathlib import Path
def fixture_source(source):
 """Insert the A/B selector after the setter, regardless of line formatting."""
 matches=list(re.finditer(r'\bt\s*\.\s*setProlongateZigzagEncoding\s*\(\s*true\s*\)\s*;',source))
 if len(matches)!=1:
  raise ValueError(f'Expected one prolongation zigzag setter in production fixture; found {len(matches)}')
 end=matches[0].end()
 selector='\n   if(std::getenv("EMI_ENTROPY_TEST_A")){t.setRestrictTransformAfterQuantization(true);t.setProlongateTransformAfterQuantization(true);}'
 return source[:end]+selector+source[end:]


def main():
 p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True);p.add_argument('--overlay',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=False)
 source=fixture_source(Path(__file__).with_name('test_bddc_split_recenter.cpp').read_text())
 cpp=a.output/'observer-check.cpp';cpp.write_text(source)
 compiler=shlex.split(os.environ.get('CXX','c++'))
 for name,root in [('original',a.original),('observer',a.overlay)]:
  subprocess.run(compiler+['-std=c++20','-O2','-I',str(root),str(cpp),'-o',str(a.output/name)],check=True)
 env=dict(os.environ,KASKADE_LIFTING_OPTIMIZATION='4');env.pop('EMI_MESSAGE_ENTROPY_PREFIX',None);env.pop('EMI_ENTROPY_TEST_A',None)
 checks=[]
 cases=[('B',*c) for c in itertools.product([8,12,16,24,32],[1,4],[0,1],[0,1],range(5))]+[('A',16,d,s,0,k) for d,s,k in itertools.product([1,2,3,4],[0,1],range(5))]
 for method,bits,depth,split,recenter,kind in cases:
  name=f'{method}-b{bits}-d{depth}-s{split}-r{recenter}-k{kind}'
  paths=[]
  for mode in ('original','observer'):
   dest=a.output/name/mode;dest.mkdir(parents=True)
   runenv=dict(env)
   if method=='A':runenv['EMI_ENTROPY_TEST_A']='1'
   if mode=='observer':runenv.update(EMI_MESSAGE_ENTROPY_PREFIX=str(dest/'counts/entropy'),EMI_ENTROPY_MAX_BINS='500000')
   subprocess.run([str(a.output/mode),str(bits),str(depth),str(split),str(recenter),str(kind),str(dest)],env=runenv,check=True,stdout=subprocess.DEVNULL)
   paths.append(dest)
  for filename in ('state.bin','packets.bin'):
   assert paths[0].joinpath(filename).read_bytes()==paths[1].joinpath(filename).read_bytes(),(name,filename)
  checks.append(name)
 (a.output/'parity.json').write_text(json.dumps(dict(status='PASS',cases=len(checks),checks=checks,verified=['Exact reconstructed state bytes','Exact original encoded packet bytes'],fixture_sha256=hashlib.sha256(source.encode()).hexdigest()),indent=2)+'\n')
 print(f'PASS: entropy observer preserves exact production packets and reconstruction in {len(checks)} cases',flush=True)

if __name__=='__main__':main()
