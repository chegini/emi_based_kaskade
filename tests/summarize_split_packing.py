#!/usr/bin/env python3
"""Report actual routed-byte savings separately from logical inter-node savings."""
import csv
import json
import re
from pathlib import Path
import sys
root=Path(sys.argv[1])
rows=[]
for split in ('on','off'):
 for order in ('before','after'):
  report=json.loads((root/f'study-split-{split}'/order/'summary.json').read_text())
  runs=report['runs']
  baseline=runs['fcg-uncompressed'].get('packet_transport',{})
  raw=baseline.get('actual_bytes')
  for bits in (8,12,16,24,32):
   r=runs[f'fcg-{order}-{bits}-split-{split}']
   packet=r.get('packet_transport',{})
   valid=r['status']=='PASS' and runs['fcg-uncompressed']['status']=='PASS'
   actual=packet.get('actual_bytes') if valid else None
   original=packet.get('original_bytes') if valid else None
   old=None
   if valid:
    log=(root/f'study-split-{split}'/order/f'fcg-{order}-{bits}-split-{split}'/'run.log').read_text()
    matches=re.findall(r'^BDDC packet small metadata bytes: previous=(\d+) actual=(\d+) saved=(\d+) messages=(\d+)$',log,re.M)
    if len(matches)!=1: raise ValueError('Missing small-metadata summary')
    old,new,saved,messages=map(int,matches[0])
    if new!=actual or old-new!=saved or saved<0: raise ValueError('Invalid small-metadata savings')
   rows.append(dict(split=split,order=order,bits=bits,status=r['status'],
    raw_routed_bytes=raw,compressed_routed_before_packing=original,packed_routed_bytes=actual,
    total_routed_saving_percent=100*(1-actual/raw) if valid and raw and actual is not None else None,
    additional_packing_saving_percent=100*(1-actual/original) if valid and original and actual is not None else None,
    logical_inter_node_saving_percent=r.get('inter_node_bytes_saved_percent') if valid else None,
    previous_packed_routed_bytes=old,
    new_metadata_bytes_saved=old-actual if valid else None,
    new_metadata_saving_percent=100*(1-actual/old) if valid and old else None,
    total_saving_gain_percentage_points=100*(old-actual)/raw if valid and raw else None,
    iterations=r.get('iterations_per_step') if valid else None))
with (root/'packed-byte-summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Actual routed bytes (inter-rank, including intra-node MPI):',root/'packed-byte-summary.csv')
