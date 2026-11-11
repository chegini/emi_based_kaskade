#!/usr/bin/env python3
"""Test the experiment matrix, split gate, Slurm launch and failure handling."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_emi_solver_study as study


def value(cmd, key):
    return cmd[cmd.index('--'+key)+1]


class SolverStudy(unittest.TestCase):
    def test_separate_executable_in_serial_and_slurm_commands(self):
        for container in ([], ['--slurm-container', '/image.sif']):
            a = study.arguments(['--repo', '/repo', '--executable',
                                 '/repo/output/batch/build/emiModel', '--solvers', 'fcg',
                                 '--bits', '12', '24', '--split', 'on']+container)
            for case in study.planned_cases(a):
                cmd = study.make_command(a, case, Path('/tmp/result'), ['mpirun'])
                self.assertIn('/repo/output/batch/build/emiModel', cmd)
                self.assertNotIn('/repo/emiModel', cmd)

    def test_lego_weak_mesh_paths_and_cluster_layout(self):
        for size, nodes in ((32, 2), (64, 4), (128, 8), (256, 16)):
            mesh = f'weak_{size}_h16'
            a = study.arguments(['--case', mesh, '--solvers', 'fcg', '--split', 'on',
                                 '--bits', '12', '24', '--expected-nodes', str(nodes),
                                 '--ranks', str(2*nodes), '--slurm-container', '/image.sif'])
            inputs = study.mesh_inputs(a.repo, mesh)
            base = a.repo.parent/'input_emi_mesh_blocks'/'weak_strong_scaling'
            for key, suffix in (('input', '_mesh.vtu'), ('extra_set', '_list_extracellular.txt'),
                                ('intra_set', '_list_intracellular.txt'), ('excited', '_early_excited.txt')):
                self.assertEqual(inputs[key], str(base/(mesh+'_kaskade'+suffix)))
            cases = study.planned_cases(a)
            self.assertEqual(len(cases), 4)
            for case in cases:
                cmd = study.make_command(a, case, Path('/tmp/result'), ['mpirun'])
                self.assertIn('--ntasks='+str(1 if case['direct'] else 2*nodes), cmd)
                self.assertIn('--nodes='+str(1 if case['direct'] else nodes), cmd)
                for path in inputs.values():
                    self.assertIn(path, cmd)
        with self.assertRaises(ValueError):
            study.mesh_inputs(Path('/repo'), 'weak_256_h16_4_8_8')

    def test_setup_diagnostic_only_runs_uncompressed_with_trace(self):
        for container in ([], ['--slurm-container', '/image.sif']):
            a = study.arguments(['--case', 'kermit', '--solvers', 'fcg',
                                 '--setup-diagnostic-only', '--case-minutes', '60']+container)
            cases = study.planned_cases(a)
            self.assertEqual([c['name'] for c in cases], ['fcg-uncompressed'])
            cmd = study.make_command(a, cases[0], Path('/tmp/result'), ['mpirun'])
            self.assertEqual(value(cmd, 'bddcCompression'), '0')
            self.assertEqual(value(cmd, 'bddcFlexibleCg'), '1')
            self.assertEqual(value(cmd, 'bddcMpiRouting'), '1')
            self.assertEqual(value(cmd, 'bddcProlongateConstantSplit'), '0')
            if container:
                self.assertIn('--time=60', cmd)
                self.assertIn('export EMI_MATERIAL_TRACE=1', cmd[cmd.index('-c')+1])
            else:
                self.assertIn('EMI_MATERIAL_TRACE=1', cmd)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            study.arguments(['--setup-diagnostic-only'])
        self.assertEqual(study.failure_kind(
            'legacy BDDC local space omits a nonzero material coupling endpoint.', False),
            'material local-space support mismatch')

    def test_diagnostic_success_makes_no_reference_error_claim(self):
        a = study.arguments(['--solvers', 'fcg', '--setup-diagnostic-only'])
        case = study.planned_cases(a)[0]
        metrics = dict(errors=[], iterations_per_step=[44], relative_residuals=[1e-8],
                       transmitted_sent_bytes=50, dofs=2)
        routes = [dict(route='inter-node', bytes=50)]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path/'run.log').write_text('EMI model completed')
            for step in (0, 1):
                (path/f'potential-state-step{step}.csv').write_text('dof,u\n0,1\n1,2\n')
            with patch.object(study, 'collect', return_value=metrics), \
                 patch.object(study, 'routing_counts', return_value=routes):
                row = study.analyze(a, case, path, 0, False, {}, None)
                self.assertEqual(row['status'], 'PASS')
                self.assertNotIn('error_vs_direct', row)
                self.assertNotIn('logical_bytes_saved_percent', row)
                metrics['errors'] = ['Missing normal completion marker']
                row = study.analyze(a, case, path, 143, False, {}, None)
                self.assertEqual(row['status'], 'FAIL')

    def test_transform_orders_have_distinct_cases_and_correct_codec_flags(self):
        names = {}
        for order in ('before', 'after', 'both'):
            a = study.arguments(['--transform-order', order])
            cases = study.planned_cases(a)
            self.assertEqual(len(cases), 34)
            names[order] = {c['name'] for c in cases if c['bits'] is not None}
            for case in cases:
                compressed = case['bits'] is not None
                cmd = study.make_command(a, case, Path('/tmp/result'), ['mpirun'])
                expected = {'bddcGraphLifting': compressed and order in ('before', 'both'),
                            'bddcGraphLiftingAfterQuantization': compressed and order in ('after', 'both'),
                            'bddcGraphLiftingBothOrders': compressed and order == 'both'}
                for key, enabled in expected.items():
                    self.assertEqual(value(cmd, key), str(int(enabled)))
                for key in ('bddcZigzag', 'bddcBitlength', 'bddcHuffman'):
                    self.assertEqual(value(cmd, key), str(int(compressed)))
                self.assertEqual(value(cmd, 'bddcCgConstantSplitExperiment'), str(int(
                    compressed and case['split'] == 'on' and case['solver'] in ('cg', 'fcg'))))
                self.assertEqual(value(cmd, 'bddcLocalCompressionBypass'), '0')
        self.assertFalse(names['before'] & names['after'])
        self.assertFalse(names['before'] & names['both'])

    def test_geometry_selection_and_split_only(self):
        for mesh in ('pepe', 'rizzo', 'robin', 'kermit', 'gonzo', 'animal', '10cells'):
            a = study.arguments(['--case', mesh, '--bits', '8', '32', '--solvers', 'fcg', '--split', 'on'])
            cases = study.planned_cases(a)
            self.assertEqual(len(cases), 4)
            self.assertEqual([c['name'] for c in cases],
                ['direct', 'fcg-uncompressed', 'fcg-before-8-split-on', 'fcg-before-32-split-on'])
            expected = study.mesh_inputs(a.repo, mesh)
            if mesh in ('pepe', 'rizzo', 'robin'):
                self.assertTrue(expected['input'].endswith(mesh+'_combi_domi.vtu'))
                self.assertTrue(expected['excited'].endswith(mesh+'_combi_domi_excited.txt'))
            elif mesh != '10cells':
                self.assertTrue(expected['input'].endswith(mesh+'_mesh.vtu'))
                self.assertTrue(expected['excited'].endswith(mesh+'_early_excited.txt'))
            for c in cases:
                cmd = study.make_command(a, c, Path('/tmp/result'), ['mpirun'])
                for key,path in expected.items():
                    self.assertEqual(value(cmd,key),path)
                self.assertEqual(value(cmd,'bddcProlongateConstantSplit'),str(int(c['bits'] is not None)))
                self.assertEqual(value(cmd,'maximumNumberOfTimeSteps'),'1')
        with self.assertRaises(ValueError):
            study.mesh_inputs(Path('/repo'), 'not-a-mesh')
        a = study.arguments(['--solvers','fcg','--bits','8','--split','off'])
        self.assertEqual(len(study.planned_cases(a)),3)
        self.assertEqual(study.planned_cases(a)[-1]['name'],'fcg-before-8-split-off')

    def test_matrix_and_solver_controls(self):
        a = study.arguments([])
        cases = study.planned_cases(a)
        self.assertEqual(len(cases), 34)
        self.assertEqual(len({c['name'] for c in cases}), 34)
        for c in cases:
            cmd = study.make_command(a, c, Path('/tmp/result'), ['mpirun'])
            compressed = c['bits'] is not None
            split = compressed and c['split'] == 'on'
            self.assertEqual(value(cmd, 'bddcCompression'), str(int(compressed)))
            self.assertEqual(value(cmd, 'bddcProlongateConstantSplit'), str(int(split)))
            self.assertEqual(value(cmd, 'bddcCgConstantSplitExperiment'),
                             str(int(split and c['solver'] != 'fgmres')))
            self.assertEqual(value(cmd, 'bddcFlexibleCg'), str(int(not c['direct'] and c['solver'] == 'fcg')))
            self.assertEqual(value(cmd, 'bddcFgmres'), str(int(not c['direct'] and c['solver'] == 'fgmres')))
            for key in ('bddcGraphLifting', 'bddcZigzag', 'bddcBitlength', 'bddcHuffman'):
                self.assertEqual(value(cmd, key), str(int(compressed)))
            for key in ('bddcLocalCompressionBypass', 'bddcGraphLiftingAfterQuantization',
                        'bddcGraphLiftingBothOrders', 'bddcHistogram', 'vtk'):
                self.assertEqual(value(cmd, key), '0')
            self.assertEqual(value(cmd, 'potentialDiagnostics'), '1')
            self.assertEqual(value(cmd, 'maximumNumberOfTimeSteps'), '1')
            if compressed:
                self.assertEqual(value(cmd, 'bddcRestrictCompressionBits'), str(c['bits']))
                self.assertEqual(value(cmd, 'bddcProlongateCompressionBits'), str(c['bits']))

    def test_slurm_and_direct_layout(self):
        for nodes in (2, 4, 8):
            a = study.arguments(['--slurm-container', '/image.sif', '--expected-nodes', str(nodes),
                                 '--ranks', str(2*nodes)])
            cases = study.planned_cases(a)
            for c, expected in ((cases[0], (1, 1)), (cases[-1], (nodes, 2*nodes))):
                cmd = study.make_command(a, c, Path('/tmp/result'), ['mpirun'])
                self.assertEqual(cmd[0], 'srun')
                self.assertNotIn('mpirun', cmd)
                self.assertNotIn('--cleanenv', cmd)
                for flag in ('--mpi=pmix_v5', '--cpu-bind=cores', '--time=20',
                             '--nodes='+str(expected[0]), '--ntasks='+str(expected[1])):
                    self.assertIn(flag, cmd)
                script = cmd[cmd.index('-c')+1]
                self.assertIn('OMPI_MCA_btl=self,vader,tcp', script)
                self.assertIn('OMPI_MCA_btl_vader_single_copy_mechanism=none', script)

    def test_failed_compressed_continues_but_bad_reference_stops(self):
        for fail_reference in (False, True):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root/'emiModel').write_text('fixture')
                (root/'emiModel').chmod(0o700)
                a = study.arguments(['--repo', tmp, '--kaskade', tmp, '--case', '10cells',
                    '--ranks', '1', '--expected-nodes', '1', '--bits', '8', '--solvers', 'fcg',
                    '--output', str(root/'results')])
                called = []

                def fake_analyze(args, case, *unused):
                    called.append(case['name'])
                    fail = case['name'] == ('fcg-uncompressed' if fail_reference else 'fcg-before-8-split-off')
                    return dict(case, status='FAIL' if fail else 'PASS',
                                reason='fixture failure' if fail else None, errors=['fixture'] if fail else [])

                with patch.object(study, 'digest', return_value='fixture-hash'), \
                     patch.object(study, 'mesh_inputs', return_value={}), \
                     patch.object(study, 'execute', return_value=(0, False, 0.1)), \
                     patch.object(study, 'analyze', side_effect=fake_analyze), \
                     contextlib.redirect_stdout(io.StringIO()):
                    code = study.campaign(a)
                report = json.loads((root/'results/summary.json').read_text())
                self.assertEqual(code, 2 if fail_reference else 1)
                self.assertEqual(report['complete'], not fail_reference)
                self.assertEqual(len(called), 2 if fail_reference else 4)
                self.assertEqual(len(report['failed_cases']), 1)

    def test_analysis_and_matched_byte_denominators(self):
        a = study.arguments(['--ranks', '2', '--expected-nodes', '1'])
        case = study.planned_cases(a)[4]  # CG, 8 bits, split off
        metrics = dict(errors=[], iterations_per_step=[10], relative_residuals=[1e-8],
                       transmitted_sent_bytes=50, dofs=2)
        routes = [dict(route='same-node', bytes=50), dict(route='inter-node', bytes=0)]
        refs = dict(initial=[0, 0], direct=[1, 2], cg=[1, 2])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path/'run.log').write_text('EMI model completed')
            for step, values in [(0, [0, 0]), (1, [1, 2.01])]:
                (path/f'potential-state-step{step}.csv').write_text('dof,u\n0,'+str(values[0])+'\n1,'+str(values[1])+'\n')
            with patch.object(study, 'collect', return_value=metrics), \
                 patch.object(study, 'routing_counts', return_value=routes):
                row = study.analyze(a, case, path, 0, False, refs,
                    dict(transmitted_sent_bytes=200, **{'inter-node_bytes': 0}))
                self.assertEqual(row['status'], 'PASS')
                self.assertEqual(row['logical_bytes_saved_percent'], 75)
                self.assertIsNone(row['inter_node_bytes_saved_percent'])
                self.assertGreater(row['error_vs_direct']['relative_l2'], 0)
                (path/'potential-state-step1.csv').unlink()
                row = study.analyze(a, case, path, 0, False, refs, {})
                self.assertEqual(row['status'], 'FAIL')
        self.assertEqual(study.failure_kind('CANCELLED DUE TO TIME LIMIT', False), 'runtime limit')
        self.assertEqual(study.failure_kind('solver breakdown', False), 'solver breakdown')

    def test_master_worker_commands_and_guards(self):
        for mode in (0, 1):
            a = study.arguments(['--master-worker', str(mode), '--solvers', 'fcg',
                                 '--bits', '16', '--split', 'off', '--ranks', '4',
                                 '--expected-nodes', '2', '--slurm-container', '/image.sif'])
            for case in study.planned_cases(a):
                cmd = study.make_command(a, case, Path('/tmp/result'), ['mpirun'])
                self.assertEqual(value(cmd, 'masterWorker'), str(0 if case['direct'] else mode))
        for args in (['--master-worker', '1', '--solvers', 'fcg', '--split', 'on'],
                     ['--master-worker', '1', '--solvers', 'fgmres'],
                     ['--master-worker', '1', '--solvers', 'fcg', '--ranks', '1', '--expected-nodes', '1']):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                study.arguments(args)

    def test_invalid_matrix(self):
        for opts in (['--bits', '8', '8'], ['--ranks', '15'], ['--case-minutes', '0']):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                study.arguments(opts)


if __name__ == '__main__':
    unittest.main()
