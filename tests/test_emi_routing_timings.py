#!/usr/bin/env python3
import unittest
from pathlib import Path

from run_emi_mpi_routing_check import arguments, make_command
from run_emi_routing_timings import aggregate, check_reference, schedule, timing_command, reference_arguments


class TimingTests(unittest.TestCase):
    def test_layout_inherited_from_reference(self):
        config = dict(case='pepe', ranks=4, expected_nodes=2, bits=16, iterations=1000,
                      slurm_btl='self,vader,tcp', slurm_vader_single_copy='none',
                      transform_order='before', solvers=['cg', 'fcg', 'fgmres'])
        for nodes, ranks in [(2, 4), (4, 8)]:
            a = reference_arguments(dict(config, ranks=ranks, expected_nodes=nodes),
                                    Path('/repo'), Path('/image.sif'), Path('/kaskade'))
            cmd = timing_command(a, Path('/out'), 'fcg-compressed-split-local-bypass')
            self.assertIn('--nodes='+str(nodes), cmd)
            self.assertIn('--ntasks='+str(ranks), cmd)
            self.assertIn('--ntasks-per-node=2', cmd)
            self.assertEqual(cmd[cmd.index('--bddcLocalCompressionBypass')+1], '1')
        for nodes, ranks in [(4, 6), (4, 2), (1, 2)]:
            with self.assertRaises(ValueError):
                reference_arguments(dict(config, ranks=ranks, expected_nodes=nodes),
                                    Path('/repo'), Path('/image.sif'), Path('/kaskade'))

    def test_snapshot_is_only_application_change(self):
        a = arguments(['--slurm-container', '/image.sif', '--case', 'pepe',
                       '--ranks', '4', '--expected-nodes', '2',
                       '--slurm-vader-single-copy', 'none'])
        for job in schedule(['cg', 'fcg', 'fgmres'], 1, 2, 7)[:9]:
            name = job['case']
            original = make_command(a, Path('/out'), ['mpirun'], name.split('-')[0],
                                    '-uncompressed-' not in name, True,
                                    bypass=name.endswith('-local-bypass'))
            timed = timing_command(a, Path('/out'), name)
            changed = [i for i, pair in enumerate(zip(original, timed)) if pair[0] != pair[1]]
            self.assertEqual(len(original), len(timed))
            self.assertEqual(changed, [original.index('--potentialDiagnostics')+1])
            self.assertEqual(timed[changed[0]], '0')

    def test_balanced_reproducible_schedule(self):
        jobs = schedule(['cg', 'fcg', 'fgmres'], 1, 3, 11)
        self.assertEqual(jobs, schedule(['cg', 'fcg', 'fgmres'], 1, 3, 11))
        self.assertEqual(len(jobs), 36)
        for offset in range(0, 36, 9):
            self.assertEqual(len({j['case'] for j in jobs[offset:offset+9]}), 9)
        self.assertTrue(all(j['warmup'] for j in jobs[:9]))
        self.assertFalse(any(j['warmup'] for j in jobs[9:]))

    def test_warmups_and_failures_excluded(self):
        rows = []
        for seconds, warmup, status in [(900, True, 'PASS'), (2, False, 'PASS'),
                                        (4, False, 'PASS'), (999, False, 'FAIL')]:
            rows.append(dict(case='cg-uncompressed-destination', warmup=warmup, status=status,
                             iterations_per_step=[41], relative_residuals=[1e-8],
                             solver_wrapper_max_seconds=seconds, wall_seconds=2*seconds,
                             launcher_seconds=3*seconds, transmitted_sent_bytes=300,
                             routes=[dict(route='same-rank', bytes=50),
                                     dict(route='same-node', bytes=50),
                                     dict(route='inter-node', bytes=200)]))
        data = aggregate(rows)['cg-uncompressed-destination']
        self.assertEqual(data['samples'], 2)
        self.assertEqual(data['solver_wrapper_max_seconds']['median'], 3)
        self.assertEqual(data['inter_node_saved_percent_vs_raw'], 0)

    def test_reference_mismatch_stops_measurement(self):
        reference = dict(cells=10, dofs=20, subdomains=3, coarse_interfaces={},
                         iterations_per_step=[4], raw_sent_bytes=40,
                         transmitted_sent_bytes=20, routes=[])
        check_reference(reference, reference)
        for key, changed in [('iterations_per_step', [5]), ('transmitted_sent_bytes', 21)]:
            with self.assertRaisesRegex(ValueError, key):
                check_reference(dict(reference, **{key: changed}), reference)


if __name__ == '__main__':
    unittest.main()
