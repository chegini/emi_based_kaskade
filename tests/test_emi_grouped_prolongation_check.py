import json
from pathlib import Path
import shlex
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import run_emi_grouped_prolongation_check as experiment
from run_emi_mixed_precision_validation import command
from run_emi_sd_cg_comparison import make_command as sd_command


def seed(model):
    saved = command(model, ['mpirun'], Path('/mesh/10cells'), 2, 10, 8,
                    'full-r16-p16', Path('/old'), matched_original_rho=True,
                    potential_diagnostics=True)
    return sd_command(saved, 'full-cg16', Path('/old'))


def log(steps=1, count=600, calls=21):
    rows = ['cells: 100', 'dofs: 2', 'BDDC subdomains: 20, shared dofs: 1',
            'BDDC coarse interfaces: corners=0, edges=1, faces=0']
    rows += [f'BDDC global residual step {k}: abs=1e-10 rel=7e-8 target=1e-7 status=converged iterations={calls}'
             for k in range(1, steps+1)]
    rows += ['raw sent bytes: 10000', f'transmitted sent bytes: {count}',
             'Profile (wall seconds, EMI end-to-end phases; 8 rank(s)):',
             '  solver wrapper (including setup) 1.0 1.0 1.0',
             'total cpu-time: 2s wall, 1s CPU', 'EMI model completed']
    return '\n'.join(rows)+'\n'


class GroupedExperimentTest(unittest.TestCase):
    def test_matched_commands(self):
        saved = seed(Path('/repo/emiModel'))
        for mode in experiment.MODES:
            cmd = experiment.make_command(saved, mode, Path('/new'))
            for flag in ('--input', '--extra_set', '--intra_set', '--excited', '--refine',
                         '--sigma_i', '--sigma_e', '--C_m', '--R', '--dt'):
                self.assertEqual(experiment.value(cmd, flag), experiment.value(saved, flag))
            for flag, val in (('--finalTime', '0.01'), ('--maximumNumberOfTimeSteps', '1'),
                              ('--bddcRestrictionPasses', '1'), ('--bddcConstantCorrection', '0'),
                              ('--bddcCgDiagnostics', '0'), ('--bddcProlongationDiagnostics', '0'),
                              ('--bddcProlongateZeroPreserving', '0')):
                self.assertEqual(experiment.value(cmd, flag), val)
            candidate = mode == experiment.MODES[-1]
            self.assertEqual(experiment.value(cmd, '--bddcProlongateGroupSize'), '32' if candidate else '0')
            for flag in ('--bddcRestrictCompressionBits', '--bddcProlongateCompressionBits'):
                self.assertEqual(experiment.value(cmd, flag), '8' if candidate else '16')
            self.assertEqual(experiment.value(cmd, '--bddcFlexibleCg'), str(int(mode in experiment.MODES[2:])))
            if mode == 'direct':
                self.assertEqual(experiment.value(cmd, '-np'), '1')
                for flag in ('--bddc', '--bddcCompression', '--bddcCgOriginalResidual', '--bddcMaterialAssembly'):
                    self.assertEqual(experiment.value(cmd, flag), '0')
            else:
                self.assertEqual(experiment.value(cmd, '--bddcRelativeTolerance'), '1e-7')
                self.assertEqual(experiment.value(cmd, '-np'), '8')
        broken = saved+['--bddcProlongateGroupSize', '16']
        with self.assertRaises(ValueError):
            experiment.make_command(broken, experiment.MODES[-1], Path('/new'))

    def test_group16_changes_only_group_size(self):
        saved = seed(Path('/repo/emiModel'))
        for a, b in zip(experiment.modes_for(32), experiment.modes_for(16)):
            cmd32 = experiment.make_command(saved, a, Path('/new'), 32)
            cmd16 = experiment.make_command(saved, b, Path('/new'), 16)
            changed = {cmd32[i-1] for i, (x, y) in enumerate(zip(cmd32, cmd16)) if x != y}
            self.assertEqual(changed, {'--bddcProlongateGroupSize'} if a == experiment.MODES[-1] else set())
        with self.assertRaises(ValueError): experiment.modes_for(8)

    def campaign(self, failure=None, missing_marker=False, parity=False, group_size=32):
        modes = experiment.modes_for(group_size)
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp).resolve(); model = repo/'emiModel'; model.write_text('mock')
            ref = repo/'old'; source = ref/'10cells/full-cg16'; source.mkdir(parents=True)
            (source/'command.sh').write_text(shlex.join(seed(model)))
            (source/'run.log').write_text(log(10))
            original = {p.name: p.read_bytes() for p in source.iterdir()}
            calls = []

            def launch(cmd, **kwargs):
                out = Path(experiment.value(cmd, '--dir')); mode = out.name; calls.append(mode)
                text = ''
                if mode in modes[2:]: text += experiment.MARKER+'\n'
                if mode == modes[-1] and not missing_marker:
                    text += experiment.group_marker(group_size)+'\n'
                if failure == mode:
                    text += 'BDDC step 1 failed after 17 solver calls: solver breakdown. State was not advanced.\n'
                elif mode == 'direct':
                    text += 'Direct global residual step 1: abs=1e-16 rel=1e-15\nEMI model completed\n'
                else:
                    text += log(count=400 if mode == modes[-1] else 600,
                                calls=22 if parity and mode == 'full-cg16' else 21)
                kwargs['stdout'].write(text)
                (out/'potential-state-step0.csv').write_text('dof,u\n0,0\n1,0\n')
                (out/'potential-state-step1.csv').write_text(
                    'dof,u\n0,1\n1,2\n' if mode == 'direct' else 'dof,u\n0,1.00000001\n1,2\n')
                return SimpleNamespace(returncode=134 if failure == mode else 0)

            argv = ['runner', '--repo', str(repo), '--reference-campaign', str(ref),
                    '--output', str(repo/'new'), '--group-size', str(group_size)]
            with patch.object(sys, 'argv', argv), patch.object(experiment.subprocess, 'run', side_effect=launch):
                rc = experiment.main()
                with self.assertRaises(FileExistsError): experiment.main()
            for p in source.iterdir(): self.assertEqual(original[p.name], p.read_bytes())
            rows = json.loads((repo/'new/summary.json').read_text())
            self.assertEqual(json.loads((repo/'new/environment.json').read_text())['prolongation_group_size'], group_size)
            self.assertIn(f'P8 groups of {group_size}', (repo/'new/summary.txt').read_text())
            return rc, rows, calls

    def test_group16_campaign(self):
        rc, rows, calls = self.campaign(group_size=16)
        self.assertEqual(rc, 0)
        self.assertEqual(calls, list(experiment.modes_for(16)))
        self.assertTrue(rows['full-fcg8-pgroup16']['beats_both_byte_baselines'])

    def test_fresh_controls_and_total_bytes(self):
        rc, rows, calls = self.campaign()
        self.assertEqual(rc, 0); self.assertEqual(calls, list(experiment.MODES))
        candidate = rows[experiment.MODES[-1]]
        self.assertTrue(candidate['beats_both_byte_baselines'])
        self.assertAlmostEqual(candidate['saved_percent']['full-cg16'], 100/3)
        self.assertGreater(candidate['full_precision_vs_direct']['relative_l2'], 0)
        self.assertTrue(rows['full-cg16']['source_first_step_parity'])

    def test_candidate_failure_retained_without_savings(self):
        rc, rows, calls = self.campaign(failure=experiment.MODES[-1])
        self.assertEqual(rc, 1); self.assertEqual(len(calls), 4)
        candidate = rows[experiment.MODES[-1]]
        self.assertEqual(candidate['failure_solver_call'], 17)
        self.assertIsNone(candidate['saved_percent']['full-cg16'])

    def test_control_failure_stops(self):
        rc, rows, calls = self.campaign(failure='full-cg16')
        self.assertEqual(rc, 1); self.assertEqual(calls, ['direct', 'full-cg16'])

    def test_missing_group_marker_and_parity_rejected(self):
        rc, rows, _ = self.campaign(missing_marker=True)
        self.assertEqual(rc, 1)
        self.assertIn('Grouped prolongation marker', ' '.join(rows[experiment.MODES[-1]]['errors']))
        rc, _, calls = self.campaign(parity=True)
        self.assertEqual(rc, 1); self.assertEqual(calls, ['direct', 'full-cg16'])


if __name__ == '__main__':
    unittest.main()
