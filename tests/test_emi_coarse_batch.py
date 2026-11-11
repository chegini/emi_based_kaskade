import json
from pathlib import Path
import tempfile
import unittest

from compare_emi_coarse_batch import CASES, compare


class Comparison(unittest.TestCase):
    def fixture(self, root):
        root.mkdir()
        config = dict(case='weak_32_h16', solvers=['fcg'], bits=[12, 24],
                      split='on', transform_order='before', ranks=4, expected_nodes=2,
                      iterations=1000, steps=1, input_sha256={'mesh': 'same'},
                      slurm_btl='self,vader,tcp', slurm_vader_single_copy='none',
                      executable_sha256=root.name)
        runs = {}
        for case in ('direct', *CASES):
            folder = root/case
            folder.mkdir()
            (folder/'potential-state-step0.csv').write_text('dof,u\n0,0\n1,0\n')
            (folder/'potential-state-step1.csv').write_text('dof,u\n0,1\n1,2\n')
            runs[case] = dict(status='PASS', cells=4, dofs=2, subdomains=2,
                coarse_interfaces={'faces': 1}, relative_residuals=[1e-8],
                iterations_per_step=[4], transmitted_sent_bytes=100)
        (root/'environment.json').write_text(json.dumps(config))
        (root/'summary.json').write_text(json.dumps(dict(complete=True, runs=runs,
                                                       failed_cases=[], fatal_error=None)))

    def test_reports_numerical_changes_without_hiding_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            old, new = Path(tmp)/'old', Path(tmp)/'new'
            for root in (old, new):
                self.fixture(root)
            path = new/CASES[1]/'potential-state-step1.csv'
            path.write_text('dof,u\n0,1.000001\n1,2\n')
            result = compare(old, new)
            self.assertEqual(result['rows'][0]['change_from_scalar']['relative_l2'], 0)
            self.assertGreater(result['rows'][1]['change_from_scalar']['relative_l2'], 0)
            self.assertGreater(result['rows'][1]['new_error']['centered_relative_l2'], 0)
            self.assertEqual(result['candidate_executable_sha256'], 'new')

    def test_missing_16_bit_baseline_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            old, new = Path(tmp)/'old', Path(tmp)/'new'
            for root in (old, new):
                self.fixture(root)
            config = json.loads((new/'environment.json').read_text())
            config['bits'] = [12, 16]
            (new/'environment.json').write_text(json.dumps(config))
            summary = json.loads((new/'summary.json').read_text())
            summary['runs']['fcg-before-16-split-on'] = summary['runs'].pop('fcg-before-24-split-on')
            (new/'summary.json').write_text(json.dumps(summary))
            (new/'fcg-before-24-split-on').rename(new/'fcg-before-16-split-on')
            rows = compare(old, new)['rows']
            self.assertTrue(rows[1]['baseline_available'])
            self.assertEqual(rows[2]['case'], 'fcg-before-16-split-on')
            self.assertFalse(rows[2]['baseline_available'])
            self.assertIsNone(rows[2]['change_from_scalar'])
            self.assertIsNone(rows[2]['old_error'])
            self.assertEqual(rows[2]['new_error']['relative_l2'], 0)

    def test_rejects_input_changes_failure_and_nonfinite_snapshot(self):
        for invalid in ('inputs', 'failed', 'nan'):
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as tmp:
                old, new = Path(tmp)/'old', Path(tmp)/'new'
                for root in (old, new):
                    self.fixture(root)
                if invalid == 'inputs':
                    p = new/'environment.json'
                    d = json.loads(p.read_text())
                    d['input_sha256']['mesh'] = 'changed'
                    p.write_text(json.dumps(d))
                elif invalid == 'failed':
                    p = new/'summary.json'
                    d = json.loads(p.read_text())
                    d['runs'][CASES[0]]['status'] = 'FAIL'
                    p.write_text(json.dumps(d))
                else:
                    (new/CASES[0]/'potential-state-step1.csv').write_text('dof,u\n0,nan\n1,2\n')
                with self.assertRaises(ValueError):
                    compare(old, new)


if __name__ == '__main__':
    unittest.main()
