"""Offline pilot configuration/topology checks; execution is mocked or refused."""
import contextlib
import csv
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_emi_partition_pilot as pilot


class Pilot(unittest.TestCase):
    def test_matrix_and_plan_does_not_execute(self):
        with patch.object(pilot.hb.study,'execute',side_effect=AssertionError('Unexpected launch')):
            output=io.StringIO()
            with contextlib.redirect_stdout(output):self.assertEqual(pilot.main(['--plan-only']),0)
        plan=json.loads(output.getvalue());c=plan['configurations']
        self.assertEqual(len(c),12)
        self.assertEqual({x['policy'] for x in c},{'P1','P2','P3'})
        self.assertTrue(all(x['bits']==16 and x['depth']==1 and x['split']=='on' for x in c if x['variant']!='raw'))
        self.assertEqual(plan['primary_trajectories']+plan['diagnostic_trajectories']+plan['validation_trajectories'],108)

    def test_refuses_outside_allocation(self):
        with patch.dict(os.environ,{},clear=True):
            with self.assertRaisesRegex(RuntimeError,'approved'):pilot.main([])

    def test_actual_command_options_and_pools(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=pilot.arguments(['--output',tmp+'/results','--intra-nodes','0,1,2,3','--extra-nodes','4,5,6,7'])
            for case in pilot.cases():
                cfg,cmd=pilot.command(a,case,Path(tmp)/case['name'],True)
                def value(flag):return cmd[cmd.index('--'+flag)+1]
                self.assertEqual(cmd[0],'srun');self.assertIn('--distribution=block:block',cmd)
                self.assertEqual(value('bddcPartitionMode'),case['partition_mode'])
                self.assertEqual(value('bddcFlexibleCg'),'1')
                self.assertEqual(value('bddcRelativeTolerance'),'1e-07')
                self.assertEqual(value('bddcTransformFirstRecentering'),str(int(case['variant']=='B-recenter')))
                self.assertEqual(value('bddcCgDiagnostics'),'0')
                self.assertEqual('--bddcIntracellularNodes' in cmd,case['policy']=='P3')
                self.assertEqual(cfg.ranks,16)

    def test_runtime_endpoint_reconciliation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            pilot.writecsv(p/'partition-interfaces.csv',[dict(rank_a=0,rank_b=1,route=2,interface_dofs=5)])
            rows=[dict(rank=rank,scope='inter_node',phase_applications=3,raw_reference_bytes=168,mpi_messages=3) for rank in [0,1]]
            pilot.validate_topology(p,rows)
            rows[0]['raw_reference_bytes']=167
            with self.assertRaisesRegex(AssertionError,'endpoints'):pilot.validate_topology(p,rows)


if __name__=='__main__':unittest.main()
