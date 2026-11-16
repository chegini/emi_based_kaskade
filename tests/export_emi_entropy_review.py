#!/usr/bin/env python3
"""Copy compact finished-job reports from scratch to one inode-friendly tar file."""
import argparse,datetime,hashlib,io,json,subprocess,tarfile
from pathlib import Path

DEFAULT_JOBS=[6184431,6184432,6184571,6184572,6184573,6184574,6184575]
PRECISIONS={6184431:12,6184432:8,6184571:8,6184572:12,6184573:16,6184574:24,6184575:32}
REPORTS=['report.txt','depth-split-report.txt','manifest.json','runs.json','errors.json',
 'combined-comparison.csv','depth-split-comparison.csv','phase-comparison.csv',
 'depth-split-paired.csv','timestep-comparison.csv','observed-prefixes.csv',
 'interface-depth-summary.csv','interface-level-shrinkage.csv',
 'components-by-step-level-rank.csv','message-entropy-summary.csv',
 'entropy-comparison.csv','status.json']
LIMIT=32*1024*1024


def export(scratch,repo,destination,jobs,include_sacct=True):
 scratch=Path(scratch);repo=Path(repo);destination=Path(destination)
 destination.parent.mkdir(parents=True,exist_ok=True)
 if destination.exists():raise FileExistsError('Refusing to replace existing review archive: '+str(destination))
 index=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scratch_root=str(scratch),
  note='Compact reports/logs and cached-validation residual histories only. Original results retained on scratch.',jobs=[],files=[],omitted=[])
 if include_sacct:
  try:
   result=subprocess.run(['sacct','-X','-j',','.join(map(str,jobs)),'--parsable2','--format=JobID,State,ExitCode,Elapsed'],capture_output=True,text=True,timeout=30)
   index['sacct']=dict(exit_code=result.returncode,stdout=result.stdout,stderr=result.stderr)
  except (OSError,subprocess.TimeoutExpired) as error:index['sacct']=dict(error=str(error))
 partial=destination.with_suffix(destination.suffix+'.partial')
 with partial.open('xb') as stream,tarfile.open(fileobj=stream,mode='w') as archive:
  included=set()
  def add(path,name,tail=False):
   path=Path(path)
   if name in included or not path.is_file():return False
   size=path.stat().st_size
   if size>LIMIT and not tail:
    index['omitted'].append(dict(path=str(path),reason='File exceeds 32 MiB review limit',bytes=size));return False
   with path.open('rb') as f:
    if tail and size>131072:f.seek(-131072,2)
    data=f.read()
   info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644
   archive.addfile(info,io.BytesIO(data));included.add(name)
   index['files'].append(dict(archive_path=name,source=str(path),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),truncated_to_tail=bool(tail and size>len(data))))
   return True
  for job in jobs:
   before=len(index['files']);j=dict(job_id=job,expected_bits=PRECISIONS.get(job),roots={})
   for label,prefix in [('performance','emi-hb-depth-split'),('entropy','emi-hb-entropy')]:
    root=scratch/f'{prefix}-{job}';results=root/'results';base=f'jobs/{job}/{label}'
    j['roots'][label]=dict(path=str(root),exists=root.is_dir(),results_exist=results.is_dir())
    for name in REPORTS:add(results/name,f'{base}/results/{name}')
    for name in ['source-sha256.txt','allocation.txt','performance-results.txt','mpi-preflight.log']:
     add(root/name,f'{base}/{name}',tail=name.endswith('.log'))
    for name in ['build/verified-build.json','observer-parity/parity.json']:
     add(root/name,f'{base}/{name}')
    add(root/'build/gccerr.txt',f'{base}/build/gccerr-tail.txt',tail=True)
    runs_file=results/'runs.json'
    if runs_file.is_file():
     for r in json.loads(runs_file.read_text()):
      # Only use a basename under this known job result root.
      folder=results/Path(r.get('folder',r.get('configuration',''))).name
      if label=='performance' and r.get('validation') and r.get('reuse') and r.get('status')=='PASS':
       for trace in sorted(folder.glob('bddc-convergence-step*.csv')):
        add(trace,f'{base}/results/{folder.name}/{trace.name}')
      if r.get('status')=='FAIL':
       add(folder/'run.log',f'{base}/results/{folder.name}/run-tail.log',tail=True)
    if label=='entropy' and results.is_dir():
     for file in sorted(results.glob('*/summary/entropy-provenance.json')):
      add(file,f'{base}/results/{file.relative_to(results).as_posix()}')
   for prefix in ['emi-hb-entropy','emi-hb-depth-split','emi-hb-8','emi-hb-12']:
    for root in [scratch,repo/'output']:
     add(root/f'{prefix}-{job}.log',f'jobs/{job}/logs/{prefix}-{job}-{root.name}-tail.log',tail=True)
   j['exported_files']=len(index['files'])-before;index['jobs'].append(j)
  if not index['files']:raise RuntimeError('No review files found under '+str(scratch))
  data=(json.dumps(index,indent=2)+'\n').encode();info=tarfile.TarInfo('review-index.json');info.size=len(data);info.mode=0o644
  archive.addfile(info,io.BytesIO(data))
 # Verify every archived report before making the completed archive visible.
 with tarfile.open(partial,'r:') as archive:
  for row in index['files']:
   with archive.extractfile(row['archive_path']) as f:
    if hashlib.sha256(f.read()).hexdigest()!=row['sha256']:raise RuntimeError('Review archive hash mismatch')
 partial.rename(destination)
 print('Review archive:',destination)
 print('Files:',len(index['files']),'Archive bytes:',destination.stat().st_size)
 for job in index['jobs']:print('Job',job['job_id'],':',job['exported_files'],'review files')
 return index


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--scratch-root',type=Path,default=Path('/scratch/htc/fchegini/emi-hb'))
 p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
 p.add_argument('--output',type=Path);p.add_argument('--jobs',type=int,nargs='+',default=DEFAULT_JOBS)
 p.add_argument('--no-sacct',action='store_true');a=p.parse_args()
 timestamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
 output=a.output or a.repo/'output'/f'emi-hb-entropy-review-{timestamp}.tar'
 export(a.scratch_root,a.repo,output,a.jobs,not a.no_sacct)
if __name__=='__main__':main()
