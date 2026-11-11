import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import run_emi_solver_study as study
import run_emi_lossless_fullsolve as check


class Packing(unittest.TestCase):
    def test_option_only_compressed_unsplit_cases(self):
        cfg=study.arguments(['--solvers','fcg','--split','off','--bits','16','--lossless-packing','1'])
        for case in study.planned_cases(cfg):
            cmd=study.make_command(cfg,case,Path('/tmp/check'),['mpirun'])
            self.assertEqual('--bddcLosslessPacking' in cmd, case['bits'] is not None)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            study.arguments(['--lossless-packing','1','--split','on'])

    def test_initial_trace_has_no_step_length(self):
        first=dict(time_step='1', iteration='0', step_length='', residual_norm='1',
                   rhs_norm='1',relative_residual='1',target='1e-7')
        last=dict(first,iteration='1',step_length='0.5',relative_residual='1e-8')
        check.validate_trace([first,last])
        for broken in (dict(last,step_length=''),dict(last,residual_norm='nan'),
                       dict(last,iteration='3')):
            with self.assertRaises(ValueError): check.validate_trace([first,broken])

    def test_parity_rejects_numerical_or_counter_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name in ('off','on'):
                (root/name).mkdir()
                for filename in check.FILES:
                    (root/name/filename).write_text('same finite samples\n')
            old=dict(folder=root/'off', iterations_per_step=[40], relative_residuals=[1e-8],
                     original_bytes=100, actual_bytes=100, residual_lines=['residual'],
                     routes=[1], raw_sent_bytes=200, transmitted_sent_bytes=80)
            new=dict(old, folder=root/'on', actual_bytes=96)
            self.assertAlmostEqual(check.compare(old,new),4)
            for key in ('iterations_per_step','relative_residuals','original_bytes','residual_lines',
                        'routes','raw_sent_bytes','transmitted_sent_bytes'):
                with self.subTest(key=key), self.assertRaises(ValueError):
                    check.compare(old,dict(new,**{key:None}))
            for filename in check.FILES:
                path=new['folder']/filename
                original=path.read_bytes()
                path.write_text('changed')
                with self.assertRaises(ValueError): check.compare(old,new)
                path.write_bytes(original)


if __name__=='__main__': unittest.main()
