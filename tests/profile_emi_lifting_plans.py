#!/usr/bin/env python3
"""Native production-method benchmark; interface storage is a Dune-free fixture.
Frozen graph connectivity/exact masks are real; weights are unit, as in the
existing frozen replay. This is not a live EMI performance/convergence result.
"""
import argparse,subprocess,shlex,os,json,hashlib,platform,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--kaskade',type=Path,required=True);p.add_argument('--test',type=Path,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--repeats',type=int,default=3);p.add_argument('--copies',type=int,default=1000)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
from lifting_performance_fixture import fixture
fixture=fixture(a.kaskade,a.test)
original=subprocess.run
extra=Path(__file__).with_name('profile_emi_lifting_main.cpp').read_text()
src=a.output/'profile.cpp';exe=a.output/'profile';src.write_text(fixture+extra)
headers={}
for name in ['bddcLiftingPlan.hh','bddcLiftingPerformance.hh','bddcTransferTiming.hh','bddcRecenter.hh']:
    header=a.kaskade/'mg'/name
    if not header.is_file():continue
    target=a.output/'include/mg'/name;target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(header,target)
    headers[name]=hashlib.sha256(target.read_bytes()).hexdigest()
(a.output/'manifest.json').write_text(json.dumps(dict(platform=platform.platform(),
    optimization=os.getenv('KASKADE_LIFTING_OPTIMIZATION','0'),
    timing=os.getenv('KASKADE_TRANSFER_TIMING','0'),
    input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),
    fixture_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),headers=headers,
    compiler=shlex.split(os.getenv('CXX','c++')),flags=['-std=c++20','-O2','-Wall','-Wextra','-Werror'],
    repeats=a.repeats,copies=a.copies,unit_weights=True,live_FCG=False),indent=2)+'\n')
original(shlex.split(os.getenv('CXX','c++'))+['-std=c++20','-O2','-Wall','-Wextra','-Werror','-I',str(a.output/'include'),'-I',str(a.kaskade),str(src),'-o',str(exe)],check=True)
with (a.output/'measurements.csv').open('w') as out:
 original([str(exe.resolve()),str(a.input.resolve()),str(a.repeats),str(a.copies)],cwd=a.output,stdout=out,check=True)
print(a.output/'measurements.csv')
