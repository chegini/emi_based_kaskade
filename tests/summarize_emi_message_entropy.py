#!/usr/bin/env python3
"""Exact empirical H0; bounded binary merge, no optional Python database modules."""
import argparse,csv,heapq,itertools,json,math,struct,tempfile
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
STAGES={0:'quantized_input',1:'encoder_coefficients',2:'huffman_bitlength',3:'routed_bytes'}
ROUTES={0:'same_rank',1:'same_node',2:'inter_node'}
MAGIC=b'EMIENT01'
WORD=struct.Struct('<QQ')


def records(stream,path):
 previous=-1
 while data:=stream.read(16*8192):
  if len(data)%16:raise ValueError('Truncated entropy record: '+str(path))
  for key,n in WORD.iter_unpack(data):
   if n<=0 or key<=previous:raise ValueError('Invalid/unsorted entropy chunk: '+str(path))
   previous=key
   yield key,n


def merge(paths):
 """Open at most the configured fan-in; aggregate equal full 64-bit keys."""
 with ExitStack() as stack:
  streams=[]
  for path in paths:
   stream=stack.enter_context(path.open('rb'))
   if stream.read(8)!=MAGIC:raise ValueError('Invalid entropy chunk header: '+str(path))
   streams.append(records(stream,path))
  for key,group in itertools.groupby(heapq.merge(*streams),key=lambda row:row[0]):
   n=sum(row[1] for row in group)
   if n>0xffffffffffffffff:raise OverflowError('Entropy count exceeds uint64')
   yield key,n


def bounded_merge(paths,scratch,fanin):
 """Multipass external merge; original chunks are always retained."""
 if fanin<2:raise ValueError('Entropy merge fan-in must be >=2')
 paths=list(paths);generation=0
 while len(paths)>fanin:
  next_paths=[]
  for start in range(0,len(paths),fanin):
   outpath=scratch/f'merge-{generation}-{start//fanin}.bin'
   with outpath.open('wb') as out:
    out.write(MAGIC)
    for key,n in merge(paths[start:start+fanin]):out.write(WORD.pack(key,n))
   next_paths.append(outpath)
  for old in paths:
   if old.parent==scratch:old.unlink()
  paths=next_paths;generation+=1
 yield from merge(paths)


def summarize(folder,output,merge_fanin=32):
 folder=Path(folder);output=Path(output);output.mkdir(parents=True,exist_ok=True)
 packet=defaultdict(lambda:dict(n=0,messages=0,weighted=0.));ranks=set()
 for file in sorted(folder.glob('entropy-rank-*-totals.csv')):
  ranks.add(int(file.name.split('-rank-')[1].split('-')[0]))
  with file.open() as f:
   for r in csv.DictReader(f):
    key=(int(r['bits'])<<56)|(int(r['flags'])<<48)|(int(r['phase'])<<40)|(int(r['stage'])<<32)
    s=packet[key];s['n']+=int(r['symbols']);s['messages']+=int(r['messages']);s['weighted']+=float(r['weighted_entropy_bits'])
 if not ranks:raise RuntimeError('No entropy totals: '+str(folder))
 paths=sorted(folder.glob('entropy-rank-*-part-*.bin'))
 frequency_records=0
 for file in paths:
  size=file.stat().st_size
  if size<8 or (size-8)%16:raise ValueError('Truncated entropy chunk: '+str(file))
  frequency_records+=(size-8)//16
 stats={key:dict(total=0,unique=0,terms=[]) for key in packet}
 with tempfile.TemporaryDirectory(prefix='entropy-merge-',dir=output) as td:
  for key,n in bounded_merge(paths,Path(td),merge_fanin):
   bucket=key&0xffffffff00000000
   if bucket not in stats:raise ValueError('Frequency bucket has no message totals')
   s=stats[bucket];s['total']+=n;s['unique']+=1;s['terms'].append(n*math.log2(n))
   if len(s['terms'])==8192:s['terms'][:]=[math.fsum(s['terms'])]
 rows=[]
 for bucket,t in sorted(packet.items()):
  s=stats[bucket];total=s['total'];unique=s['unique']
  if total!=t['n'] or total<=0:raise ValueError(f'Frequency total mismatch: {bucket}, {total}, {t}')
  h=math.log2(total)-math.fsum(s['terms'])/total;conditional=t['weighted']/total
  if not (h>=-1e-10 and conditional>=-1e-10 and h+1e-8>=conditional):raise ValueError(f'Invalid entropy: {h}, {conditional}')
  flags=(bucket>>48)&255;stage=(bucket>>32)&255;phase=(bucket>>40)&255;bits=(bucket>>56)&255
  rows.append(dict(phase=['restrict','prolongate'][phase],stage=STAGES[stage],route=ROUTES[flags>>3],
   source_bits=bits,floating=bool(flags&1),integer=bool(flags&2),encoded=bool(flags&4),
   symbols=total,messages=t['messages'],distinct_symbols=unique,
   pooled_shannon_bits_per_symbol=max(0.,h),message_weighted_shannon_bits_per_symbol=max(0.,conditional),
   unit='bits/byte' if stage==3 else 'bits/length_symbol' if stage==2 else 'bits/coefficient',
   ideal_message_entropy_bits=t['weighted']))
 with (output/'entropy.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['phase','stage','symbols']);w.writeheader();w.writerows(rows)
 (output/'entropy-provenance.json').write_text(json.dumps(dict(ranks=sorted(ranks),frequency_records=frequency_records,
  aggregation='Bounded sorted binary merge; Python standard library only, no SQLite.',merge_fanin=merge_fanin,
  definition='H0=-sum p*log2(p). Message-weighted H0=sum N_message*H0_message / sum N_message.',
  exclusions='Coefficient entropy excludes exact corners, split constants, headers and synthetic codebook seeds.',
  note='Routed-byte entropy includes real serialized R/P headers/corners/constants; its unit is bits/byte, not bits/coefficient.'),indent=2)+'\n')
 return rows,sorted(ranks)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();rows,ranks=summarize(a.folder,a.output);print('Entropy rows:',len(rows),'ranks:',ranks)
