import json
import contextlib
import io
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_emi_hb_time_study as run
import summarize_emi_hb_time_study as report


class Analysis(unittest.TestCase):
    def test_campaign_orchestration_with_synthetic_files(self):
        # This tests orchestration/validation/accounting, NOT PDE performance.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);exe=root/'fake-model';exe.write_text('synthetic')
            payload=root/'input.txt';payload.write_text('synthetic input fingerprint')
            steps=2;ranks=2

            def fake_execute(cmd,folder,cfg,env):
                value=lambda flag:cmd[cmd.index('--'+flag)+1]
                compressed=value('bddcCompression')=='1'
                tier=int(env['KASKADE_TRANSFER_TIMING'])
                row_bytes=50 if compressed else 100
                (folder/'run.log').write_text(f'BDDC packet transport bytes: original=400 actual={row_bytes*4} packing={int(compressed)}\n')
                for rank in range(ranks):
                    rows=[]
                    for step in range(steps+1):
                        rows.append(dict(step=step,calls=2 if step else 0,wall_s=1.,fcg_s=.5 if step else 0,
                            relative_residual=1e-8 if step else 0,original_bytes=100 if step else 0,
                            packed_bytes=row_bytes if step else 0,pack_s=0,unpack_s=0,exchange_s=0,
                            verify_s=0,packet_solve_s=0,inter_node_packed_bytes=row_bytes if step else 0))
                    with (folder/f'hb-time-rank-{rank}.csv').open('w',newline='') as f:
                        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
                    (folder/f'hb-time-rank-{rank}-total.csv').write_text('simulation_s,application_compute_s,application_setup_s\n3,4,2\n')
                    stage='step,call,phase,stage,level,calls,seconds\n'
                    if tier and compressed:stage+='0,0,other,plan_setup,0,1,0.1\n'
                    (folder/f'hb-time-rank-{rank}-stages.csv').write_text(stage)
                if value('potentialDiagnostics')=='1':
                    for step in range(steps+1):
                        (folder/f'potential-state-step{step}.csv').write_text(f'dof,u\n0,{step+1}\n')
                        if step:(folder/f'bddc-convergence-step{step}.csv').write_text('iteration,relative\n2,1e-8\n')
                return 0,False,1.

            def fake_collect(*args):
                return dict(errors=[],iterations_per_step=[2,2],wall_seconds=4.,dofs=10,subdomains=4)

            with patch.dict('os.environ',SLURM_JOB_NUM_NODES='2'),patch.object(run.study,'execute',fake_execute),\
                 patch.object(run,'collect',fake_collect),patch.object(run,'routing_counts',return_value=[]),\
                 patch.object(run,'mesh_inputs',return_value={'input':str(payload)}),contextlib.redirect_stdout(io.StringIO()):
                status=run.main(['--executable',str(exe),'--kaskade','/tmp/k',
                    '--slurm-container','/tmp/s','--output',str(root/'results'),
                    '--ranks','2','--nodes','2','--steps','2'])
            self.assertEqual(status,0)
            self.assertTrue((root/'results/PASS').is_file())
            with (root/'results/combined-comparison.csv').open() as f:combined=list(csv.DictReader(f))
            self.assertEqual(len(combined),7)
            self.assertTrue(all(r['complete']=='True' for r in combined))
            (root/'results/steps.json').unlink();(root/'results/stages.json').unlink()
            recovered=report.summarize(root/'results')
            self.assertEqual(len(recovered),21)

    def test_optimization_campaign_with_synthetic_files(self):
        # This tests orchestration/validation/accounting, NOT PDE performance.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);exe=root/'fake-model';exe.write_text('synthetic')
            payload=root/'input.txt';payload.write_text('synthetic input fingerprint')
            steps=2;ranks=2

            def fake_execute(cmd,folder,cfg,env):
                value=lambda flag:cmd[cmd.index('--'+flag)+1]
                compressed=value('bddcCompression')=='1'
                tier=int(env['KASKADE_TRANSFER_TIMING'])
                row_bytes=50 if compressed else 100
                (folder/'run.log').write_text(f'BDDC packet transport bytes: original=400 actual={row_bytes*4} packing={int(compressed)}\n')
                for rank in range(ranks):
                    rows=[]
                    for step in range(steps+1):
                        rows.append(dict(step=step,calls=2 if step else 0,wall_s=1.,fcg_s=.5 if step else 0,
                            relative_residual=1e-8 if step else 0,original_bytes=100 if step else 0,
                            packed_bytes=row_bytes if step else 0,pack_s=0,unpack_s=0,exchange_s=0,
                            verify_s=0,packet_solve_s=0,inter_node_packed_bytes=row_bytes if step else 0))
                    with (folder/f'hb-time-rank-{rank}.csv').open('w',newline='') as f:
                        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
                    (folder/f'hb-time-rank-{rank}-total.csv').write_text('simulation_s,application_compute_s,application_setup_s\n3,4,2\n')
                    stage='step,call,phase,stage,level,calls,seconds\n'
                    if tier and compressed:stage+='0,0,other,plan_setup,0,1,0.1\n'
                    (folder/f'hb-time-rank-{rank}-stages.csv').write_text(stage)
                if value('potentialDiagnostics')=='1':
                    for step in range(steps+1):
                        (folder/f'potential-state-step{step}.csv').write_text(f'dof,u\n0,{step+1}\n')
                        if step:(folder/f'bddc-convergence-step{step}.csv').write_text('iteration,relative\n2,1e-8\n')
                return 0,False,1.

            def fake_collect(*args):
                return dict(errors=[],iterations_per_step=[2,2],wall_seconds=4.,dofs=10,subdomains=4)

            with patch.dict('os.environ',SLURM_JOB_NUM_NODES='2'),patch.object(run.study,'execute',fake_execute),\
                 patch.object(run,'collect',fake_collect),patch.object(run,'routing_counts',return_value=[]),\
                 patch.object(run,'mesh_inputs',return_value={'input':str(payload)}),contextlib.redirect_stdout(io.StringIO()):
                status=run.main(['--executable',str(exe),'--kaskade','/tmp/k',
                    '--slurm-container','/tmp/s','--output',str(root/'results'),
                    '--ranks','2','--nodes','2','--steps','2',
                    '--lifting-optimizations','0','1','4','5'])
            self.assertEqual(status,0)
            self.assertTrue((root/'results/PASS').is_file())
            with (root/'results/combined-comparison.csv').open() as f:combined=list(csv.DictReader(f))
            self.assertEqual(len(combined),25)
            self.assertTrue(all(r['complete']=='True' for r in combined))
            (root/'results/steps.json').unlink();(root/'results/stages.json').unlink()
            recovered=report.summarize(root/'results')
            self.assertEqual(len(recovered),75)
            with (root/'results/optimization-comparison.csv').open() as f:comparisons=list(csv.DictReader(f))
            self.assertEqual(len(comparisons),54)
            self.assertTrue(all(r['packed_bytes_equal']=='True' and r['iterations_equal']=='True' for r in comparisons))
            self.assertTrue(all(float(r['simulation_difference_s_mean'])==0 for r in comparisons))

    def test_thresholds_and_uncertainty(self):
        self.assertEqual(report.projected_break_even(10.,2.),6)
        self.assertEqual(report.projected_break_even(-1.,2.),1)
        self.assertIsNone(report.projected_break_even(10.,0.))
        self.assertEqual(report.projected_break_even(-1.,0.),1)
        self.assertIsNone(report.projected_break_even(10.,-2.))
        interval=report.mean_ci([1.,2.,3.])
        self.assertLess(interval['low'],0.)
        self.assertGreater(interval['high'],4.)
        self.assertIsNone(report.mean_ci([1.])['low'])

    def test_commands_and_weak_family(self):
        a=run.arguments(['--executable','/tmp/model','--kaskade','/tmp/kaskade',
            '--slurm-container','/tmp/image.sif','--output','/tmp/out','--steps','10'])
        cases=run.configurations(a)
        self.assertEqual(len(cases),7)
        for case in cases:
            _,cmd=run.command(a,case,Path('/tmp/out')/case['name'])
            self.assertIn(str(a.output.resolve().parent)+':'+str(a.output.resolve().parent),cmd)
            values=lambda flag:cmd[cmd.index('--'+flag)+1]
            self.assertEqual(values('maximumNumberOfTimeSteps'),'10')
            self.assertEqual(values('potentialDiagnostics'),'0')
            self.assertEqual(values('bddcConvergenceReport'),'0')
            self.assertEqual(values('bddcFlexibleCg'),'1')
            self.assertEqual(values('bddcProlongateConstantSplit'),'0')
            self.assertEqual(values('bddcReuseLiftingPlans'),'1')
            self.assertEqual(values('bddcTransformFirstRecentering'),str(int(case['variant']=='B-recenter')))
        with self.assertRaises(SystemExit),contextlib.redirect_stderr(io.StringIO()):
            run.arguments(['--executable','/tmp/e','--kaskade','/tmp/k',
                '--slurm-container','/tmp/s','--output','/tmp/o','--scaling','weak','--case','weak_32_h16'])

    def test_rank_max_after_accumulation_and_prefix(self):
        # Unequal rank setup/step work: adding separate maxima would fabricate
        # a speedup. The analyzer must take max only after each rank's prefix sum.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            manifest=dict(case='synthetic',nodes=2,ranks=2,steps=3,repeats=3,detail_repeats=3)
            runs=[];steps=[]
            for rep in range(3):
                for name in ['raw','B-b16-d4']:
                    run_name=f'r{rep}-{name}'
                    raw=name=='raw'
                    rows=[]
                    for rank in range(2):
                        setup=2 if rank==0 else 0
                        for step in range(4):
                            wall=setup if step==0 else ((2 if rank==0 else 3) if raw else 1)
                            row=dict(run=run_name,configuration=name,tier=0,repeat=rep,validation=False,warmup=False,
                                rank=rank,step=step,calls=0 if step==0 else 2,wall_s=wall,fcg_s=0 if step==0 else wall,
                                original_bytes=0 if step==0 else 100,packed_bytes=0 if step==0 else (100 if raw else 50),
                                inter_node_packed_bytes=0 if step==0 else (100 if raw else 50),
                                pack_s=0,unpack_s=0,exchange_s=0,verify_s=0)
                            rows.append(row)
                    steps.extend(rows)
                    total=max(sum(r['wall_s'] for r in rows if r['rank']==rank) for rank in range(2))
                    runs.append(dict(name=name,configuration=name,variant='raw' if raw else 'B',bits=None if raw else 16,depth=4,
                        tier=0,repeat=rep,validation=False,warmup=False,status='PASS',folder=str(root/run_name),
                        simulation_s=total,application_compute_s=total+1,application_setup_s=3,wall_seconds=total+2,packed_bytes=600 if raw else 300,
                        iterations_per_step=[2,2,2]))
            for name,data in [('manifest',manifest),('runs',runs),('steps',steps),('stages',[]),('errors',[])]:
                (root/(name+'.json')).write_text(json.dumps(data))
            summary=report.summarize(root)
            b=next(r for r in summary if r['configuration']=='B-b16-d4')
            self.assertEqual(b['trajectory_difference_s_mean'],-4.)
            self.assertEqual(b['saving_percent_mean'],50.)
            self.assertEqual(b['observed_ci_benefit_from_step'],1)

if __name__=='__main__':unittest.main()
