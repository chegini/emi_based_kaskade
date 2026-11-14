#!/usr/bin/env python3
"""Check hierarchy reporting against hand-counted successful/terminal splits."""
import csv
from pathlib import Path
import tempfile
import unittest
from summarize_emi_hb_depth_split import inventory

class Inventory(unittest.TestCase):
 def test_shrinkage_and_terminal_stop(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);folder=root/'validate';folder.mkdir()
   common=dict(plan_id=7,producer=3,consumer=8,role='float_peer',requested_depth=4,
    interface_dofs=19,exact_dofs=2,phase='both',edges_after=8,max_degree=2,seconds=.001,
    process_peak_rss_bytes=99999,plan_storage_bytes=256)
   events=[]
   for kind,level,n,c,f in [('predictor_split',1,17,9,8),('initial_graph',1,17,9,8),
     ('elimination',2,17,5,4),('elimination',3,9,3,2),
     ('elimination',4,5,2,1),('graph_plan',4,17,2,15)]:
    events.append(dict(common,kind=kind,level=level,nodes_before=n,coarse=c,detail=f))
   # Another interface stops after one stage; a failed attempt must not count as level 2.
   for kind,level,c,f in [('predictor_split',1,1,1),('elimination',2,1,0),('graph_plan',1,1,1)]:
    events.append(dict(common,plan_id=8,kind=kind,level=level,nodes_before=2,coarse=c,detail=f))
   with (folder/'lifting-setup-test.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(events[0]));w.writeheader();w.writerows(events)
   runs=[dict(validation=True,reuse=True,status='PASS',variant='B',configuration='B-test',folder=str(folder))]
   plans=inventory(root,runs,required=True)
   self.assertEqual([p['actual_depth'] for p in plans],[4,1])
   self.assertEqual([p['reached_requested_depth'] for p in plans],[True,False])
   with (root/'interface-level-shrinkage.csv').open() as f:rows=list(csv.DictReader(f))
   self.assertEqual([int(r['active_dofs']) for r in rows[:4]],[17,9,5,3])
   self.assertEqual([int(r['coarse_dofs']) for r in rows[:4]],[9,5,3,2])
   self.assertAlmostEqual(float(rows[3]['coarse_fraction_of_original']),2/17)
 def test_missing_inventory_is_error(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   with self.assertRaisesRegex(RuntimeError,'Missing interface'):
    inventory(root,[dict(validation=True,reuse=True,status='PASS',variant='A',configuration='A',folder=str(root))],True)
if __name__=='__main__':unittest.main()
