"""Offline guards for the four-case Pepe validation; never launch a simulation."""
import contextlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import run_check as check


class WeightAwareCheck(unittest.TestCase):
    def test_plan_is_four_matched_trajectories_without_execution(self):
        with patch.object(check.pilot.hb.study, 'execute', side_effect=AssertionError('Launch')):
            text = io.StringIO()
            with contextlib.redirect_stdout(text):
                self.assertEqual(check.main([]), 0)
        plan = json.loads(text.getvalue())
        self.assertEqual((plan['nodes'], plan['ranks'], plan['timesteps']), (2, 4, 3))
        self.assertEqual(len(plan['commands']), 4)
        self.assertEqual(plan['status'], 'PREPARED_NOT_EXECUTED')

    def test_exact_scientific_options_and_isolated_changes(self):
        a = check.arguments([])
        cfg = check.configuration(a)
        flags = []
        for case in check.cases():
            settings, cmd = check.command(cfg, case, a.output / case['name'])
            self.assertEqual(settings.steps, 3)
            self.assertIn('--nodes=2', cmd)
            self.assertIn('--ntasks=4', cmd)
            self.assertEqual(cmd[cmd.index('--finalTime') + 1], '0.03')
            values = {flag: cmd[index + 1] for index, flag in enumerate(cmd[:-1])
                      if flag.startswith('--') and '=' not in flag}
            for name, value in [('bddcFlexibleCg', '1'), ('bddcUseCg', '1'),
                                ('bddcFgmres', '0'), ('bddcGraphLiftingAfterQuantization', '0'),
                                ('bddcPartitionMode', 'baseline'), ('bddcRelativeTolerance', '1e-07')]:
                self.assertEqual(values['--' + name], value)
            flags.append(values)
        # Recoding alone must not silently change bit depth, transformation, solver or mesh.
        del flags[2]['--dir']
        del flags[3]['--dir']
        flags[2]['--bddcTransformFirstRecentering'] = '1'
        self.assertEqual(flags[2], flags[3])

    def test_refuses_to_execute_outside_allocation(self):
        with patch.dict(os.environ, {}, clear=True), \
             patch.object(check.pilot.hb.study, 'execute', side_effect=AssertionError('Launch')):
            with self.assertRaisesRegex(ValueError, 'user-submitted'):
                check.main(['--execute'])

    def test_12_bit_plan_changes_only_precision_and_output_locations(self):
        plans = []
        for precision in [16, 12]:
            text = io.StringIO()
            with patch.object(check.pilot.hb.study, 'execute', side_effect=AssertionError('Launch')), \
                 contextlib.redirect_stdout(text):
                self.assertEqual(check.main(['--bits', str(precision)]), 0)
            plan = json.loads(text.getvalue())
            self.assertEqual(plan['bit_depth'], precision)
            self.assertEqual([c['bits'] for c in plan['configurations']],
                             [None, precision, precision, precision])
            plans.append(plan)
        self.assertEqual(plans[0]['commands'][0], plans[1]['commands'][0])
        for original, lower in zip(plans[0]['commands'][1:], plans[1]['commands'][1:]):
            expected = original.copy()
            for flag in ['bddcCompressionBits', 'bddcRestrictCompressionBits',
                         'bddcProlongateCompressionBits']:
                expected[expected.index('--' + flag) + 1] = '12'
            self.assertEqual(expected, lower)

    def test_snapshot_true_residual(self):
        self.assertAlmostEqual(check.snapshot_residual(dict(rhs=[3, 4], residual=[0, 1])), .2)
        self.assertEqual(check.snapshot_residual(dict(rhs=[0], residual=[0])), 0)
        self.assertIsNone(check.snapshot_residual(dict(rhs=[0], residual=[1])))

    def test_production_application_guards_accept_generated_commands(self):
        from test_emi_compression_report import extract, kaskade_include
        repo = Path(__file__).resolve().parents[2]
        source = (repo / 'EMI_model.cpp').read_text()
        options = extract(source, 'struct EmiOptions').replace('Kaskade::BDDC::EDGE', '0') + ';'
        start = source.index('  if (options.quiet != 0')
        end = source.index('  OptionalMpiSession mpiSession', start)
        fields = dict(re.findall(r'\("(\w+)",\s*options\.(\w+)\s*,', source))
        strings = set(re.findall(r'std::string\s+(\w+)', options))
        prefix = ('#include <cmath>\n#include <iostream>\n#include <string>\n'
                  '#include "EMI_partitioning.hh"\n#include "EMI_bddc_precision.hh"\n'
                  '#include "EMI_bddc_packet_packing.hh"\n#include "EMI_bddc_step_schedule.hh"\n' +
                  options + '\nvoid validate(EmiOptions& options){\n' + source[start:end] +
                  'EmiBddc::validateStepSchedule(options,options.maximumNumberOfTimeSteps);\n}\n')
        cfg = check.configuration(check.arguments([]))
        checks = []
        for case in check.cases(16) + check.cases(12):
            cmd = check.command(cfg, case, Path('/tmp') / case['name'])[1]
            payload = cmd[cmd.index('--mpi'):]
            assignments = []
            for flag, value in zip(payload[::2], payload[1::2]):
                field = fields[flag.removeprefix('--')]
                literal = json.dumps(value)
                if field not in strings:
                    literal = 'std::stod(' + literal + ')'
                assignments.append('options.' + field + '=' + literal + ';')
            checks.append('{EmiOptions options;\n' + '\n'.join(assignments) +
                          '\ntry{validate(options);}catch(std::exception const& e){'
                          'std::cerr<<e.what();return 1;}\n}\n')
        with tempfile.TemporaryDirectory(prefix='pepe-weight-aware-guards-') as tmp:
            cpp = Path(tmp) / 'guards.cpp'
            exe = Path(tmp) / 'guards'
            cpp.write_text(prefix + 'int main(){\n' + ''.join(checks) + '}\n')
            subprocess.run(shlex.split(os.environ.get('CXX', 'c++')) + [
                '-std=c++20', '-I', str(repo), '-I', kaskade_include(repo),
                str(cpp), '-o', str(exe)], check=True, timeout=180)
            subprocess.run([str(exe)], check=True, timeout=30)

    def test_topology_and_initial_state_are_validation_invariants(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ['partition-problem', 'partition-ownership', 'partition-ranks', 'partition-interfaces']:
                (root / (name + '.csv')).write_text('mode,planning_s,value\nbaseline,1,2\n')
            for name in ['partition-coarse.csv', 'potential-operator.csv', 'potential-state-step0.csv']:
                (root / name).write_text('same\n')
            first = check.invariants(root)
            (root / 'partition-problem.csv').write_text('mode,planning_s,value\nbaseline,4,2\n')
            self.assertEqual(first, check.invariants(root))
            (root / 'partition-interfaces.csv').write_text('mode,planning_s,value\nbaseline,1,3\n')
            self.assertNotEqual(first, check.invariants(root))


if __name__ == '__main__':
    unittest.main()
