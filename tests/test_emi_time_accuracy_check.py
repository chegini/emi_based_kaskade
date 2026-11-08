#!/usr/bin/env python3
import tempfile
import unittest
import os
import subprocess
from pathlib import Path

from run_emi_time_accuracy_check import (CASES, accuracy_arguments, accuracy_command,
                                         direct_check, trajectory_errors)


class AccuracyTests(unittest.TestCase):
    def config(self, nodes):
        return dict(case='pepe', ranks=2*nodes, expected_nodes=nodes, steps=1,
                    uncompressed_steps=0, split='on', bits=16, transform_order='before',
                    slurm_btl='self,vader,tcp', slurm_vader_single_copy='none',
                    iterations=1000, solvers=['cg', 'fcg', 'fgmres'])

    def test_physics_and_transport_for_all_layouts(self):
        for nodes in (2, 4, 8):
            a = accuracy_arguments(self.config(nodes), Path('/repo'), '/image', '/kaskade', 10)
            for name in CASES:
                cmd = accuracy_command(a, Path('/out'), name)
                def value(flag):
                    self.assertEqual(cmd.count(flag), 1)
                    return cmd[cmd.index(flag)+1]
                direct = name == 'direct'
                compressed = name.startswith('fgmres-compressed')
                self.assertIn('--nodes='+str(1 if direct else nodes), cmd)
                self.assertIn('--ntasks='+str(1 if direct else 2*nodes), cmd)
                self.assertIn('--time=00:40:00', cmd)
                self.assertEqual(value('--maximumNumberOfTimeSteps'), '10')
                self.assertEqual(value('--dt'), '0.01')
                self.assertEqual(value('--finalTime'), '0.1')
                self.assertEqual(value('--potentialDiagnostics'), '1')
                self.assertEqual(value('--bddcHistogram'), '0')
                self.assertNotIn('--bddcHistogramSteps', cmd)
                self.assertEqual(value('--bddcFgmres'), str(int(not direct)))
                self.assertEqual(value('--bddcCgConstantSplitExperiment'), '0')
                self.assertEqual(value('--bddcProlongateConstantSplit'), str(int(compressed)))
                self.assertEqual(value('--bddcLocalCompressionBypass'), str(int(name.endswith('-local-bypass'))))
                self.assertEqual(value('--bddcUncompressedSteps'), '0')
                self.assertEqual(value('--bddcGraphLifting'), str(int(compressed)))
                self.assertEqual(value('--bddcGraphLiftingAfterQuantization'), '0')

    def test_experimental_limits_and_reference_settings(self):
        for steps in (0, 11):
            with self.assertRaises(ValueError):
                accuracy_arguments(self.config(2), Path('/repo'), '/image', '/kaskade', steps)
        for changes in (dict(expected_nodes=3, ranks=6), dict(bits=12),
                        dict(steps=5), dict(uncompressed_steps=1), dict(transform_order='after')):
            with self.assertRaises(ValueError):
                accuracy_arguments(dict(self.config(2), **changes), Path('/repo'), '/image', '/kaskade', 10)

    def test_direct_requires_each_step_and_finite_residual(self):
        log = ('Direct global residual step 1: abs=1e-12 rel=1e-12\n'
               'Direct global residual step 2: abs=2e-12 rel=2e-12\nEMI model completed\n')
        self.assertEqual(direct_check(log, 0, 2), [1e-12, 2e-12])
        for bad, rc in [(log, 1), (log.replace('step 2', 'step 3'), 0),
                        (log.replace('rel=2e-12', 'rel=2e-2'), 0),
                        (log.replace('rel=2e-12', 'rel=nan'), 0)]:
            with self.assertRaises(ValueError):
                direct_check(bad, rc, 2)

    def test_raw_and_centered_errors_keep_mean_offset(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            for name, values in [('candidate', [2, 3]), ('direct', [1, 2]), ('raw', [1, 2])]:
                folder = base/name
                folder.mkdir()
                for step in (1, 2):
                    # Different later RHS values are legitimate in evolving trajectories.
                    (folder/f'potential-state-step{step}.csv').write_text(
                        'dof,u,rhs\n'+''.join(f'{i},{u},{step*u}\n' for i,u in enumerate(values)))
            rows = trajectory_errors(base/'candidate', base/'direct', base/'raw', 2)
            self.assertEqual([r['time'] for r in rows], [0.01, 0.02])
            for r in rows:
                self.assertEqual(r['vs_direct']['mean_error'], 1)
                self.assertGreater(r['vs_direct']['relative_l2'], 0)
                self.assertEqual(r['vs_direct']['centered_relative_l2'], 0)
                self.assertEqual(r['vs_direct'], r['vs_uncompressed'])

    def test_batch_gates_and_eight_node_layout(self):
        script = Path(__file__).with_name('submit_emi_time_accuracy.sbatch')
        for fail_first in (False, True):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root/'bin').mkdir()
                (root/'output').mkdir()
                stub = root/'bin/python3'
                stub.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$TEST_CALLS"\n'
                                'if [[ "$TEST_FAIL_FIRST" == 1 && "$2" == tests/run_emi_mpi_routing_check.py ]]; then exit 42; fi\n')
                stub.chmod(0o755)
                calls = root/'calls.txt'
                env = dict(os.environ, PATH=str(root/'bin')+':'+os.environ['PATH'],
                           SLURM_SUBMIT_DIR=str(root), SLURM_JOB_NUM_NODES='8',
                           SLURM_NTASKS='16', SLURM_JOB_ID='test', SLURM_JOB_NODELIST='test-nodes',
                           TEST_CALLS=str(calls), TEST_FAIL_FIRST=str(int(fail_first)))
                result = subprocess.run(['bash', str(script), '--with-timings'],
                                        env=env, text=True, capture_output=True)
                rows = calls.read_text().splitlines()
                self.assertEqual(result.returncode, 42 if fail_first else 0, result.stderr)
                self.assertEqual(len(rows), 1 if fail_first else 3)
                self.assertIn('--ranks 16 --expected-nodes 8', rows[0])
                if not fail_first:
                    self.assertIn('run_emi_time_accuracy_check.py', rows[1])
                    self.assertIn('--steps 10', rows[1])
                    self.assertIn('run_emi_routing_timings.py', rows[2])


if __name__ == '__main__':
    unittest.main()
