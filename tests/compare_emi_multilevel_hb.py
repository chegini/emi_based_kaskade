#!/usr/bin/env python3
"""Compare HB depths at fixed solver, precision, packing and node count."""
import csv
import json
from pathlib import Path
import re
import sys

def summarize(root):
    rows=[]
    for depth in (1,0):
        for order in ('before','after'):
            folder=root/f'levels-{depth}'/order
            runs=json.loads((folder/'summary.json').read_text())['runs']
            base=runs['fcg-uncompressed']
            raw=base.get('packet_transport',{}).get('actual_bytes')
            for bits in (8,12,16,24,32):
                name=f'fcg-{order}-{bits}-split-on'
                r=runs[name]
                valid=r['status']=='PASS' and base['status']=='PASS' and raw
                actual=r.get('packet_transport',{}).get('actual_bytes') if valid else None
                log=(folder/name/'run.log').read_text(errors='replace')
                setup=re.findall(r'BDDC integer lifting setup rank=(\d+) local_plan_bytes=(\d+) exported_remote_bytes=(\d+)',log)
                rows.append(dict(levels=depth,order=order,bits=bits,status=r['status'],
                    iterations=r.get('iterations_per_step'),
                    relative_l2_direct=r.get('error_vs_direct',{}).get('relative_l2'),
                    actual_routed_bytes=actual,raw_routed_bytes=raw,
                    total_saving_percent=100*(1-actual/raw) if valid and actual is not None else None,
                    wrapper_seconds=r.get('solver_wrapper_max_seconds'),
                    exported_setup_bytes=sum(int(x[2]) for x in setup) if setup else None))
    for row in rows:
        old=next(r for r in rows if r['levels']==1 and r['order']==row['order'] and r['bits']==row['bits'])
        valid=row['actual_routed_bytes'] is not None and old['actual_routed_bytes'] is not None
        row['bytes_saved_vs_single']=old['actual_routed_bytes']-row['actual_routed_bytes'] if valid else None
        row['saving_gain_percentage_points']=(row['total_saving_percent']-old['total_saving_percent']) if valid else None
    with (root/'hb-comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('HB results:',root/'hb-comparison.csv')
    print('Routed bytes cover iterative R/P traffic. Setup bytes exclude collective replication. Timings include snapshots.')

if __name__=='__main__':
    summarize(Path(sys.argv[1]))
