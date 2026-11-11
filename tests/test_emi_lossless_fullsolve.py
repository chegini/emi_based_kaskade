import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

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

    def test_length_accounting_must_match_transport(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            dictionary_line='BDDC packet dictionary reference bytes: wide=88 actual=88 saved=0 messages=0 byte_refs=0 wide_refs=0\n'
            record_line='BDDC packet record length bytes: wide=92 actual=88 saved=4 records=2 shortened=2 le255=2 le4095=0 le32767=0 larger=0\n'
            count_line='BDDC packet exact count omission bytes: with_counts=96 actual=88 saved=4\n'
            endpoint='BDDC packet endpoint plan bytes: explicit=104 actual=88 plan_ids=8 omitted_ids=8 plan_messages=1 reused_messages=1\n'
            for step in (0,1):
                (folder/f'potential-state-step{step}.csv').write_text('dof,u\n0,1\n')
            (folder/check.FILES[2]).write_text(
                'time_step,iteration,step_length,residual_norm,rhs_norm,relative_residual,target\n'
                '1,0,,1,1,1,1e-7\n1,1,0.5,1e-8,1,1e-8,1e-7\n')
            cfg=SimpleNamespace(ranks=4,expected_nodes=2,lossless_packing=1)
            metrics=dict(errors=[],transmitted_sent_bytes=20,iterations_per_step=[1])
            with patch.object(study,'collect',return_value=metrics), \
                 patch.object(study,'routing_counts',return_value=[dict(bytes=20)]):
                for summary in ('wide=112 actual=88','wide=91 actual=88',
                                'wide=112 actual=90',''):
                    (folder/'run.log').write_text(
                        'BDDC packet transport bytes: original=100 actual=88 packing=1\n'+
                        ('BDDC packet length compaction bytes: '+summary+'\n' if summary else '')+
                        'BDDC packet tail omission bytes: two_lengths=104 actual=88\n'+endpoint+count_line+record_line+dictionary_line)
                    if summary=='wide=112 actual=88':
                        self.assertEqual(check.inspect(folder,cfg,0,False)['wide_bytes'],112)
                    else:
                        with self.assertRaises(ValueError): check.inspect(folder,cfg,0,False)
                for tail in ('two_lengths=104 actual=88','two_lengths=91 actual=88',
                             'two_lengths=113 actual=88','two_lengths=104 actual=90',''):
                    (folder/'run.log').write_text(
                        'BDDC packet transport bytes: original=100 actual=88 packing=1\n'
                        'BDDC packet length compaction bytes: wide=112 actual=88\n'+
                        ('BDDC packet tail omission bytes: '+tail+'\n' if tail else '')+endpoint+count_line+record_line+dictionary_line)
                    if tail=='two_lengths=104 actual=88':
                        self.assertEqual(check.inspect(folder,cfg,0,False)['two_length_bytes'],104)
                    else:
                        with self.assertRaises(ValueError): check.inspect(folder,cfg,0,False)

                # Endpoint accounting must be present, consistent, and actually reused.
                for bad in ('', endpoint+endpoint, endpoint.replace('explicit=104','explicit=101'),
                            endpoint.replace('actual=88','actual=91'),
                            endpoint.replace('plan_ids=8','plan_ids=7'),
                            endpoint.replace('reused_messages=1','reused_messages=0')):
                    (folder/'run.log').write_text(
                        'BDDC packet transport bytes: original=120 actual=88 packing=1\n'
                        'BDDC packet length compaction bytes: wide=112 actual=88\n'
                        'BDDC packet tail omission bytes: two_lengths=104 actual=88\n'+bad+count_line+record_line+dictionary_line)
                    with self.assertRaises(ValueError): check.inspect(folder,cfg,0,False)

                for bad_count in ('', count_line+count_line, count_line.replace('saved=4','saved=8'),
                                  count_line.replace('actual=88','actual=91'),
                                  count_line.replace('with_counts=96','with_counts=92')):
                    (folder/'run.log').write_text(
                        'BDDC packet transport bytes: original=120 actual=88 packing=1\n'
                        'BDDC packet length compaction bytes: wide=112 actual=88\n'
                        'BDDC packet tail omission bytes: two_lengths=104 actual=88\n'+endpoint+bad_count+record_line+dictionary_line)
                    with self.assertRaises(ValueError): check.inspect(folder,cfg,0,False)

                for bad_record in ('', record_line+record_line, record_line.replace('saved=4','saved=6'),
                                   record_line.replace('wide=92','wide=96'),
                                   record_line.replace('shortened=2','shortened=1'),
                                   record_line.replace('records=2','records=3'),
                                   record_line.replace('le255=2','le255=0')):
                    (folder/'run.log').write_text(
                        'BDDC packet transport bytes: original=120 actual=88 packing=1\n'
                        'BDDC packet length compaction bytes: wide=112 actual=88\n'
                        'BDDC packet tail omission bytes: two_lengths=104 actual=88\n'+endpoint+count_line+bad_record+dictionary_line)
                    with self.assertRaises(ValueError): check.inspect(folder,cfg,0,False)

                for bad_dict in ('', dictionary_line+dictionary_line,
                                 dictionary_line.replace('saved=0','saved=3'),
                                 dictionary_line.replace('byte_refs=0','byte_refs=1'),
                                 dictionary_line.replace('messages=0','messages=1')):
                    (folder/'run.log').write_text(
                        'BDDC packet transport bytes: original=120 actual=88 packing=1\n'
                        'BDDC packet length compaction bytes: wide=112 actual=88\n'
                        'BDDC packet tail omission bytes: two_lengths=104 actual=88\n'+endpoint+count_line+record_line+bad_dict)
                    with self.assertRaises(ValueError): check.inspect(folder,cfg,0,False)

                # Positive dictionary savings keep earlier stage baselines separate.
                log=('BDDC packet transport bytes: original=120 actual=88 packing=1\n'
                     'BDDC packet length compaction bytes: wide=112 actual=88\n'
                     'BDDC packet tail omission bytes: two_lengths=104 actual=88\n'+endpoint+count_line+record_line)
                log=log.replace('actual=88','actual=85')
                log+='BDDC packet dictionary reference bytes: wide=88 actual=85 saved=3 messages=1 byte_refs=1 wide_refs=0\n'
                (folder/'run.log').write_text(log)
                result=check.inspect(folder,cfg,0,False)
                self.assertEqual(result['endpoint']['dictionary_reference_bytes_saved'],3)

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
