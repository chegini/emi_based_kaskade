#!/usr/bin/env python3
"""Measured packet sizes by depth; coefficient ranges and raw tail bits by level.
Huffman headers, shared-codebook costs and byte padding are packet/group costs;
per-level tail bits alone are not the complete encoded byte count.
"""
import csv,collections,json,sys
from pathlib import Path
root=Path(sys.argv[1]);packet_rows=[];level_rows=[]
def read(p):
 with p.open() as f:yield from csv.DictReader(f)
for folder in sorted(root.glob('replay-b*-d*')):
 if not (folder/'level-stats.csv').exists():raise ValueError('Missing compact replay statistics: '+str(folder))
 bits,depth=folder.name.removeprefix('replay-b').split('-d')
 totals=collections.defaultdict(lambda:collections.Counter())
 for r in read(folder/'packets.csv'):
  t=totals[(r['mode'],r['phase'])]
  for k in ['raw_bytes','wire_bytes','huffman_bytes','tail_bytes','padding_bytes','exact_bytes','plan_bytes']:t[k]+=int(r[k])
  t['squared_error']+=float(r['squared_error']);t['squared_norm']+=float(r['squared_norm'])
 for (mode,phase),t in totals.items():packet_rows.append(dict(bits=bits,depth=depth,mode=mode,phase=phase,**t))
 stats={}
 for r in read(folder/'level-stats.csv'):
  if r['mode']=='quantization':continue
  key=(r['mode'],r['phase'],int(r['level']))
  v=stats.setdefault(key,dict(count=0,minimum=float(r['minimum']),maximum=float(r['maximum']),sum_squares=0.,tail_bits=0))
  v['count']+=int(r['count']);v['minimum']=min(v['minimum'],float(r['minimum']));v['maximum']=max(v['maximum'],float(r['maximum']));v['sum_squares']+=float(r['sum_squares']);v['tail_bits']+=int(r['tail_bits'])
 for (mode,phase,level),v in stats.items():
  import math
  level_rows.append(dict(bits=bits,depth=depth,mode=mode,phase=phase,level=level,count=v['count'],minimum=v['minimum'],maximum=v['maximum'],rms=math.sqrt(v['sum_squares']/v['count']),tail_bits=v['tail_bits'],scope='floating coefficients before quantization' if mode=='before' else 'signed modular integer details'))
for name,rows in [('packet-depths.csv',packet_rows),('coefficient-levels.csv',level_rows)]:
 with (root/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Audit tables:',root/'packet-depths.csv',root/'coefficient-levels.csv')
