"""Offline tests for the three-case Pepe check; no simulation is executed."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_emi_partition_smoke as smoke


class Smoke(unittest.TestCase):
    def test_byte_observer_works_with_raw_and_compressed(self):
        with tempfile.TemporaryDirectory() as tmp:
            for raw in [False,True]:
                a=smoke.arguments(['--byte-entropy','--output',tmp+'/out']+(['--uncompressed'] if raw else []))
                for case in smoke.cases(raw):
                    cfg,cmd=smoke.command(a,case,Path(tmp)/case['name'])
                    self.assertEqual(cmd[cmd.index('--bddcSerializedByteEntropy')+1],'1')
                    self.assertEqual(cmd[cmd.index('--bddcCompression')+1],str(int(not raw)))
                    self.assertEqual(cmd[cmd.index('--bddcCgDiagnostics')+1],'0')

    def test_raw_only_plan(self):
        with patch.object(smoke.pilot.hb.study,'execute',side_effect=AssertionError('Unexpected launch')):
            out=io.StringIO()
            with contextlib.redirect_stdout(out):self.assertEqual(smoke.main(['--uncompressed','--plan-only']),0)
        plan=json.loads(out.getvalue())
        self.assertEqual(plan['trajectories'],3)
        self.assertEqual([c['name'] for c in plan['configurations']],['P2-raw','P1-raw','P3-raw'])
        self.assertTrue(all(c['bits'] is None and c['split']=='off' for c in plan['configurations']))

    def test_raw_flags_disable_entire_codec(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=smoke.arguments(['--uncompressed','--output',tmp+'/raw'])
            for case in smoke.cases(True):
                cfg,cmd=smoke.command(a,case,Path(tmp)/case['name'])
                for flag in ['bddcCompression','bddcGraphLifting','bddcGraphLiftingAfterQuantization',
                             'bddcProlongateConstantSplit','bddcTransformFirstRecentering',
                             'bddcHuffman','bddcBitlength','bddcZigzag','bddcCgDiagnostics']:
                    self.assertEqual(cmd[cmd.index('--'+flag)+1],'0',flag)
                self.assertNotIn('--bddcLosslessPacking',cmd)
                self.assertEqual(cmd[cmd.index('--bddcFlexibleCg')+1],'1')
                self.assertEqual(cmd[cmd.index('--bddcRelativeTolerance')+1],'1e-07')

    def test_only_three_cases_without_launch(self):
        with patch.object(smoke.pilot.hb.study,'execute',side_effect=AssertionError('Unexpected launch')):
            out=io.StringIO()
            with contextlib.redirect_stdout(out):self.assertEqual(smoke.main(['--plan-only']),0)
        plan=json.loads(out.getvalue())
        self.assertEqual((plan['nodes'],plan['ranks'],plan['steps'],plan['trajectories']),(2,4,1,3))
        self.assertEqual([c['policy'] for c in plan['configurations']],['P2','P1','P3'])
        self.assertTrue(all(c['variant']=='B-recenter' and c['bits']==16 and c['depth']==1 and c['split']=='on' for c in plan['configurations']))

    def test_exact_solver_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=smoke.arguments(['--output',tmp+'/results'])
            for case in smoke.cases():
                cfg,cmd=smoke.command(a,case,Path(tmp)/case['name'])
                def value(flag):return cmd[cmd.index('--'+flag)+1]
                self.assertIn('--nodes=2',cmd);self.assertIn('--ntasks=4',cmd)
                self.assertIn('--ntasks-per-node=2',cmd)
                for flag,expected in [('bddcFlexibleCg','1'),('bddcUseCg','1'),
                    ('bddcGraphLifting','1'),('bddcGraphLiftingLevels','1'),
                    ('bddcGraphLiftingAfterQuantization','0'),('bddcTransformFirstRecentering','1'),
                    ('bddcProlongateConstantSplit','1'),('bddcCompressionBits','16'),
                    ('maximumNumberOfTimeSteps','1'),('finalTime','0.01'),
                    ('bddcRelativeTolerance','1e-07'),('bddcLocalCompressionBypass','0'),
                    ('bddcCgDiagnostics','0'),('bddcConvergenceReport','1')]:
                    self.assertEqual(value(flag),expected,flag)
                self.assertEqual(value('bddcPartitionMode'),case['partition_mode'])
                if case['policy']=='P3':
                    self.assertEqual(value('bddcIntracellularNodes'),'0')
                    self.assertEqual(value('bddcExtracellularNodes'),'1')
                self.assertEqual(cfg.steps,1)

    def test_refuses_outside_allocation(self):
        with patch.dict(os.environ,{},clear=True):
            with self.assertRaisesRegex(RuntimeError,'existing SLURM allocation'):smoke.main([])

    def test_placement_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);a=smoke.arguments([])
            ranks=[dict(rank=r,node_group=r//2,hostname='host'+str(r//2)) for r in range(4)]
            owners=[dict(subdomain_id=s,owner_rank=s,node_group=s//2,region='I' if s<2 else 'E') for s in range(4)]
            smoke.pilot.writecsv(folder/'partition-ranks.csv',ranks)
            smoke.pilot.writecsv(folder/'partition-ownership.csv',owners)
            smoke.check_placement(folder,smoke.cases()[2],a)
            smoke.check_placement(folder,smoke.cases()[0],a)
            owners[0]['region']='E'
            smoke.pilot.writecsv(folder/'partition-ownership.csv',owners)
            with self.assertRaisesRegex(ValueError,'P3 does not separate'):
                smoke.check_placement(folder,smoke.cases()[2],a)


if __name__=='__main__':unittest.main()
