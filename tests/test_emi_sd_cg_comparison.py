import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

from run_emi_flexible_cg_check import MARKER
from run_emi_mixed_precision_validation import command
from run_emi_sd_cg_comparison import MODES, make_command, value


class ComparisonTest(unittest.TestCase):
    def seed(self, model, launcher, out):
        return command(model, launcher, Path('/mesh'), 2, 2, 8, 'full-r16-p16', out,
                       matched_original_rho=True, matched_cg_diagnostics=True,
                       potential_diagnostics=True, matched_constant_correction=True)

    def test_matched_commands(self):
        seed = self.seed(Path('/model'), ['mpi'], Path('/old'))
        commands = {mode: make_command(seed, mode, Path('/new')) for mode in MODES}
        for cmd in commands.values():
            for flag in ['--input', '--extra_set', '--intra_set', '--excited', '--refine',
                         '--sigma_i', '--sigma_e', '--C_m', '--R', '--dt', '--finalTime',
                         '--maximumNumberOfTimeSteps', '-np']:
                self.assertEqual(value(cmd, flag), value(seed, flag))
            for flag, wanted in [('--bddcRelativeTolerance', '1e-7'), ('--bddcRestrictionPasses', '1'),
                                  ('--bddcConstantCorrection', '0'), ('--bddcCgDiagnostics', '0'),
                                  ('--bddcProlongateZeroPreserving', '0')]:
                self.assertEqual(value(cmd, flag), wanted)
        cg, sd = commands['full-cg8'], commands['full-sd8']
        changed = {cg[i-1] for i, (a, b) in enumerate(zip(cg, sd)) if a != b}
        self.assertEqual(changed, {'--bddcUseCg', '--bddcCgOriginalResidual'})
        fcg = commands['full-fcg8']
        changed = {cg[i-1] for i, (a, b) in enumerate(zip(cg, fcg)) if a != b}
        self.assertEqual(changed, {'--bddcFlexibleCg'})
        for mode in ['full-cg8', 'full-fcg8', 'full-sd8']:
            for flag in ['--bddcCompressionBits', '--bddcRestrictCompressionBits', '--bddcProlongateCompressionBits']:
                self.assertEqual(value(commands[mode], flag), '8')
            self.assertEqual(value(commands[mode], '--bddcHuffman'), '1')
        self.assertEqual(value(make_command(seed, 'full-sd8', Path('/new'), 1000), '--bddcIterations'), '1000')

    def test_failure_continues_to_sd_and_counts_total_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp).resolve(); ref = repo/'reference'; ref.mkdir()
            model = repo/'emiModel'; model.write_text('mock executable')
            fake = repo/'mpi'
            fake.write_text('''#!/usr/bin/env python3
from pathlib import Path
import sys
a=sys.argv[1:]
def v(k): return a[a.index(k)+1]
out=Path(v('--dir')); mode=out.name
with open('calls','a') as f:f.write(mode+'\\n')
if mode=='full-fcg8':print('''+repr(MARKER)+''')
if mode in ['full-cg8','full-fcg8']:
 print('BDDC step 1 failed after 7 solver calls: solver breakdown. State was not advanced.')
 sys.exit(134)
print('cells: 100\\ndofs: 2\\nBDDC subdomains: 2, shared dofs: 1')
print('BDDC coarse interfaces: corners=0, edges=1, faces=0')
for k in [1,2]:
 print(f'BDDC global residual step {k}: abs=1e-10 rel=7e-8 target=1e-7 status=converged iterations=11')
print('raw sent bytes: 10000')
print('transmitted sent bytes:',300 if mode=='full-sd8' else 600)
print('Profile (wall seconds, EMI end-to-end phases; 8 rank(s)):')
print('  solver wrapper (including setup) 1.0 1.0 1.0')
print('total cpu-time: 2s wall, 1s CPU\\nEMI model completed')
(out/'potential-state-step2.csv').write_text('dof,u\\n0,1.00000001\\n1,2\\n')
''')
            fake.chmod(0o755)
            for mode in ['direct', 'full-r16-p16']:
                folder = ref/'10cells'/mode; folder.mkdir(parents=True)
            seed = self.seed(model, [str(fake)], ref/'10cells/full-r16-p16')
            (ref/'10cells/full-r16-p16/command.sh').write_text(shlex.join(seed))
            (ref/'10cells/direct/potential-state-step2.csv').write_text('dof,u\n0,1\n1,2\n')
            (ref/'environment.json').write_text(json.dumps(dict(steps=2, bddc_ranks=8)))
            (ref/'summary.json').write_text(json.dumps({'10cells': {'runs': {
                m: {'status': 'PASS'} for m in ['direct', 'full-r16-p16']}}}))
            saved = (ref/'summary.json').read_bytes()
            script = Path(__file__).with_name('run_emi_sd_cg_comparison.py')
            args = [sys.executable, str(script), '--repo', str(repo), '--reference-campaign', str(ref),
                    '--output', str(repo/'new')]
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 1, result.stdout+result.stderr)
            rows = json.loads((repo/'new/summary.json').read_text())['10cells']
            self.assertEqual(list(rows), list(MODES))
            self.assertEqual((repo/'calls').read_text().splitlines(), list(MODES))
            self.assertEqual(rows['full-cg8']['failure_solver_call'], 7)
            self.assertEqual(rows['full-fcg8']['status'], 'FAIL')
            self.assertEqual(rows['full-sd8']['status'], 'PASS')
            self.assertEqual(rows['full-sd8']['saved_vs_cg16_percent'], 50)
            self.assertEqual(rows['full-sd8']['metrics']['iterations_per_step'], [11, 11])
            self.assertEqual((ref/'summary.json').read_bytes(), saved)
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)  # existing output is never overwritten


if __name__ == '__main__':
    unittest.main()
