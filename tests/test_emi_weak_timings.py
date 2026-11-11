import copy
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run_emi_weak_timings as timing


class WeakTimings(unittest.TestCase):
    def fixture(self, repo):
        config = dict(case='weak_32_h16', solvers=['fcg'], bits=[12, 24], split='on',
                      transform_order='before', slurm_btl='self,vader,tcp',
                      slurm_vader_single_copy='none', ranks=4, expected_nodes=2,
                      iterations=1000, executable_sha256='hash')
        config['input_sha256'] = {p: 'hash' for p in timing.mesh_inputs(repo, config['case']).values()}
        routes = [dict(route=route, bytes=10, packets=1, phase=phase)
                  for phase in ('restrict', 'prolongate')
                  for route in ('same-rank', 'same-node', 'inter-node')]
        row = dict(status='PASS', errors=[], iterations_per_step=[40], relative_residuals=[1e-8],
                   cells=100, dofs=50, subdomains=64, coarse_interfaces=[1],
                   raw_sent_bytes=100, transmitted_sent_bytes=60, routes=routes,
                   solver_wrapper_max_seconds=2, wall_seconds=3)
        names = ('direct', 'fcg-uncompressed', 'fcg-before-12-split-on', 'fcg-before-24-split-on')
        validation = dict(complete=True, failed_cases=[], fatal_error=None,
                          runs={name: copy.deepcopy(row) for name in names})
        return config, validation, row

    def test_reference_and_command_flags(self):
        config, validation, _ = self.fixture(Path('/repo'))
        a, cases = timing.reference_arguments(config, validation, Path('/repo'), Path('/image'), Path('/kaskade'))
        for case in cases:
            cmd = timing.timing_command(a, case, Path('/tmp/out'))
            for option, expected in (('potentialDiagnostics','0'), ('vtk','0'), ('bddcHistogram','0'),
                                      ('profile','1'), ('bddcLocalCompressionBypass','0'),
                                      ('bddcRelativeTolerance','1e-07')):
                self.assertEqual(cmd[cmd.index('--'+option)+1], expected)
            self.assertIn('--nodes=2', cmd)
            self.assertIn('--ntasks=4', cmd)
        validation['complete'] = False
        with self.assertRaises(ValueError):
            timing.reference_arguments(config, validation, Path('/repo'), Path('/image'), Path('/kaskade'))

    def test_separate_executable_is_reused_from_validation(self):
        config, validation, _ = self.fixture(Path('/repo'))
        config['executable'] = '/repo/output/batch/build/emiModel'
        a, cases = timing.reference_arguments(config, validation, Path('/repo'),
                                              Path('/image'), Path('/kaskade'))
        for case in cases:
            cmd = timing.timing_command(a, case, Path('/tmp/out'))
            self.assertIn(config['executable'], cmd)
            self.assertNotIn('/repo/emiModel', cmd)

    def test_12_16_reference_keeps_selected_bits(self):
        config, validation, _ = self.fixture(Path('/repo'))
        config['bits'] = [12, 16]
        validation['runs']['fcg-before-16-split-on'] = validation['runs'].pop('fcg-before-24-split-on')
        a, cases = timing.reference_arguments(config, validation, Path('/repo'),
                                              Path('/image'), Path('/kaskade'))
        self.assertEqual(a.bits, [12, 16])
        self.assertEqual([c['name'] for c in cases],
                         ['fcg-uncompressed', 'fcg-before-12-split-on', 'fcg-before-16-split-on'])
        cmd = timing.timing_command(a, cases[2], Path('/tmp/out'))
        self.assertEqual(cmd[cmd.index('--bddcCompressionBits')+1], '16')

    def test_schedule_and_cpu_rejection(self):
        cases = [dict(name=s) for s in ('raw','12','24')]
        jobs = timing.schedule(cases, 1, 3, 42)
        self.assertEqual(len(jobs), 12)
        self.assertEqual(sum(j['warmup'] for j in jobs), 3)
        self.assertEqual(jobs, timing.schedule(cases, 1, 3, 42))
        for name in ('raw','12','24'):
            self.assertEqual(sum(j['case']==name and not j['warmup'] for j in jobs), 3)
        text = '\n'.join(json.dumps(dict(host=h,cpu='Intel')) for h in ('a','b'))
        self.assertEqual(len(timing.check_cpu_inventory(text, 2, 'Intel')), 2)
        for invalid in (text.replace('Intel', 'AMD', 1), text+'\n'+text.splitlines()[0], text.splitlines()[0]):
            with self.assertRaises(ValueError):
                timing.check_cpu_inventory(invalid, 2, 'Intel')

    def test_full_campaign_and_trajectory_failure(self):
        for mismatch in (False, True):
            with self.subTest(mismatch=mismatch), tempfile.TemporaryDirectory() as tmp:
                repo = Path(tmp).resolve()
                reference = repo/'reference'
                reference.mkdir()
                image = repo/'image.sif'
                image.touch()
                config, validation, row = self.fixture(repo)
                (reference/'environment.json').write_text(json.dumps(config))
                (reference/'summary.json').write_text(json.dumps(validation))
                if mismatch:
                    row['iterations_per_step'] = [41]
                def execute(cmd, folder, a, env):
                    (folder/'run.log').write_text('test')
                    return 0, False, 4
                inventory = '\n'.join(json.dumps(dict(host=h,cpu='Intel')) for h in ('a','b'))
                with patch.object(timing, '__file__', str(repo/'tests'/'runner.py')), \
                     patch.object(timing, 'digest', return_value='hash'), \
                     patch.object(timing.subprocess, 'run', return_value=SimpleNamespace(returncode=0, stdout=inventory)), \
                     patch.object(timing.study, 'execute', side_effect=execute), \
                     patch.object(timing, 'collect', side_effect=lambda *args: copy.deepcopy(row)), \
                     patch.object(timing, 'routing_counts', return_value=row['routes']), \
                     patch.dict(os.environ, SLURM_JOB_ID='1', SLURM_JOB_NUM_NODES='2', SLURM_NTASKS='4'), \
                     contextlib.redirect_stdout(io.StringIO()):
                    rc = timing.main(['--reference', str(reference), '--slurm-container', str(image),
                                      '--kaskade', str(repo/'kaskade'), '--cpu-model', 'Intel'])
                summary = json.loads(next((repo/'output').glob('*/summary.json')).read_text())
                self.assertEqual(rc, int(mismatch))
                self.assertEqual(summary['complete'], not mismatch)
                self.assertEqual(len(summary['runs']), 1 if mismatch else 12)
                if not mismatch:
                    for data in summary['aggregate'].values():
                        self.assertEqual(data['samples'], 3)
                        self.assertEqual(data['solver_wrapper_max_seconds']['median'], 2)


if __name__ == '__main__':
    unittest.main()
