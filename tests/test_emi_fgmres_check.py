import contextlib
import io
import json
import re
from pathlib import Path
import shlex
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run_emi_fgmres_check as experiment
from test_emi_grouped_prolongation_check import seed, log


class FgmresDriverTest(unittest.TestCase):
    def test_commands(self):
        saved=seed(Path('/repo/emiModel'))
        for mode in experiment.MODES+(experiment.SPLIT_MODE,):
            cmd=experiment.make_command(saved,mode,Path('/new'))
            for flag in ('--input','--extra_set','--intra_set','--excited','--refine',
                         '--sigma_i','--sigma_e','--C_m','--R','--dt'):
                self.assertEqual(experiment.value(cmd,flag),experiment.value(saved,flag))
            for flag,wanted in (('--maximumNumberOfTimeSteps','1'),('--finalTime','0.01'),
                                ('--bddcRestrictionPasses','1'),('--bddcProlongateGroupSize','0'),
                                ('--bddcConstantCorrection','0'),('--bddcFlexibleCg','0'),
                                ('--bddcPreconditionerProbe','0'),('--bddcIterations','300'),
                                ('--bddcProlongateConstantSplit',str(int(mode==experiment.SPLIT_MODE)))):
                self.assertEqual(experiment.value(cmd,flag),wanted)
            if 'fgmres' in mode:
                for flag,wanted in (('--bddcFgmres','1'),('--bddcFgmresRestart','30'),
                                    ('--bddcCgOriginalResidual','0'),('--bddcUseCg','0'),
                                    ('--bddcRelativeTolerance','1e-7')):
                    self.assertEqual(experiment.value(cmd,flag),wanted)
            if mode in ('full-fgmres8',experiment.SPLIT_MODE):
                for flag in ('--bddcRestrictCompressionBits','--bddcProlongateCompressionBits'):
                    self.assertEqual(experiment.value(cmd,flag),'8')
                self.assertEqual(experiment.value(cmd,'--bddcRestrictZeroPreserving'),'1')
                self.assertEqual(experiment.value(cmd,'--bddcProlongateZeroPreserving'),'0')

    def test_new_flags_are_registered(self):
        source=(Path(__file__).resolve().parents[1]/'EMI_model.cpp').read_text()
        registered=set(re.findall(r'\(\s*"([^"\n]+)"\s*,\s*options\.',source))
        saved=seed(Path('/repo/emiModel'))
        cmd=experiment.make_command(saved,experiment.SPLIT_MODE,Path('/new'))
        for flag in cmd:
            if flag.startswith('--') and flag not in saved:self.assertIn(flag[2:],registered)

    def campaign(self,failure=None,missing_ledger=False,split=False,changed_control=False,missing_split_marker=False,steps=1,missing_snapshot=None,late_failure=False,changed_candidate=False,case="10cells",bits=None,pairs=None,changed_mixed=False,late_mixed_failure=False,reference_mixed_failure=False,reference_mixed_bits=6):
        with tempfile.TemporaryDirectory() as temp:
            repo=Path(temp).resolve();model=repo/'emiModel';model.write_text('mock')
            source=repo/'source/10cells/full-cg16';source.mkdir(parents=True)
            (source/'command.sh').write_text(shlex.join(seed(model)))
            (source/'run.log').write_text(log(1 if split else 10))
            if split:
                source_steps=10 if case!='10cells' else 1
                cmd=seed(model);experiment.set_value(cmd,'--maximumNumberOfTimeSteps',source_steps)
                (source/'command.sh').write_text(shlex.join(cmd))
                old_modes=experiment.MODES[:-1]+(experiment.SPLIT_MODE,) if steps==10 or bits is not None or pairs else experiment.MODES[:-1]
                if steps==10 and bits and 6 in bits:old_modes+=('full-fgmres6-split',)
                if steps==10 and pairs:old_modes+=('full-fgmres-r6-p8-split',)
                for mode in old_modes:
                    folder=repo/'source/10cells'/mode;folder.mkdir(exist_ok=True)
                    if mode in experiment.MIXED_BITS:
                        old_cmd=experiment.make_command(seed(model),mode,folder,source_steps)
                        experiment.set_value(old_cmd,'--bddcRestrictCompressionBits',reference_mixed_bits)
                        (folder/'command.sh').write_text(shlex.join(old_cmd))
                    k=30 if mode in experiment.SPLIT_PAIRS else 21
                    text=(''.join(f'Direct global residual step {step}: abs=1e-16 rel=1e-15\n'
                                  for step in range(1,source_steps+1))+'EMI model completed\n'
                          if mode=='direct' else log(source_steps,calls=k))
                    (folder/'run.log').write_text(text)
                    (folder/'potential-state-step0.csv').write_text('dof,u\n0,0\n1,0\n')
                    (folder/'potential-state-step1.csv').write_text('dof,u\n0,1\n1,2\n')
                    if 'fgmres' in mode:
                        for step in range(1,source_steps+1):
                            text+=f'BDDC FGMRES statistics step {step}: applications={k} cycles=1 root_matvecs={2*k} root_peak_basis_bytes=960\n'
                            (folder/f'bddc-fgmres-step{step}.csv').write_text('iteration,cycle,estimated_absolute_residual,true_absolute_residual,true_relative_residual\n'+
                                ''.join(f'{i},1,1e-8,1e-8,1e-8\n' for i in range(k+1)))
                        (folder/'run.log').write_text(text)
                    if mode in experiment.MIXED_BITS and reference_mixed_failure:
                        (folder/'run.log').write_text(text.replace('step 6: abs=1e-10 rel=7e-8 target=1e-7 status=converged',
                                                                 'step 6: abs=1e-10 rel=0.2 target=1e-7 status=not-converged'))
            original=(source/'run.log').read_bytes();calls=[]
            def run(cmd,**kw):
                out=Path(experiment.value(cmd,'--dir'));mode=out.name;calls.append(mode)
                self.assertEqual(out.parent.name,case)
                if case!='10cells':
                    self.assertEqual(experiment.value(cmd,'--refine'),'0')
                    for flag,suffix in experiment.INPUT_SUFFIXES.items():
                        self.assertEqual(experiment.value(cmd,flag),str(repo/'mesh'/(case+'_kaskade'))+suffix)
                fgmres='fgmres' in mode;candidate=mode in ('full-fgmres8',*experiment.SPLIT_PAIRS);k=30 if candidate else 21
                text=experiment.MARKER+'\n' if fgmres else ''
                if mode in experiment.SPLIT_PAIRS and not missing_split_marker:text+=experiment.SPLIT_MARKER+'\n'
                if mode==failure:
                    text+='BDDC step 1 failed after 5 solver calls: solver breakdown. State was not advanced.\n'
                    if fgmres:
                        (out/'bddc-fgmres-step1.csv').write_text('iteration,true_relative_residual\n0,1\n1,0.5\n2,0.7\n')
                elif mode=='direct':
                    text+=''.join(f'Direct global residual step {step}: abs=1e-16 rel=1e-15\n'
                                  for step in range(1,steps+1))+'EMI model completed\n'
                else:
                    text+=log(steps,count=400 if candidate else 600,calls=k)
                    outer=steps*(7*(k*(32+12)+4+56+16+16) if fgmres else 7*(k+1)*(32+16))
                    nc=steps*(4*k+4 if fgmres else 2*(k+1))
                    if not missing_ledger:text+=f'BDDC outer solver logical payload bytes: {outer} collectives={nc}\n'
                    if fgmres:
                        for step in range(1,steps+1):
                            text+=f'BDDC FGMRES statistics step {step}: applications={k} cycles=1 root_matvecs={2*k} root_peak_basis_bytes=960\n'
                            (out/f'bddc-fgmres-step{step}.csv').write_text('iteration,cycle,estimated_absolute_residual,true_absolute_residual,true_relative_residual\n'+
                                ''.join(f'{i},1,1e-8,1e-8,1e-8\n' for i in range(k+1)))
                    if (candidate and late_failure) or (mode in experiment.MIXED_BITS and late_mixed_failure):
                        text=text.replace('step 6: abs=1e-10 rel=7e-8 target=1e-7 status=converged',
                                          'step 6: abs=1e-10 rel=0.2 target=1e-7 status=not-converged')
                kw['stdout'].write(text)
                (out/'potential-state-step0.csv').write_text('dof,u\n0,0\n1,0\n')
                for step in range(1,steps+1):
                    if candidate and step==missing_snapshot:continue
                    delta=.01*step if candidate and step>1 else 0
                    (out/f'potential-state-step{step}.csv').write_text(f'dof,u\n0,{1+delta}\n1,{2+delta}\n')
                if (changed_candidate and candidate) or (mode in experiment.MIXED_BITS and changed_mixed):
                    (out/'potential-state-step1.csv').write_text('dof,u\n0,1.1\n1,2\n')
                if changed_control and mode=='full-fgmres16':
                    (out/'potential-state-step1.csv').write_text('dof,u\n0,1.1\n1,2\n')
                return SimpleNamespace(returncode=134 if mode==failure else 0)
            argv=['runner','--repo',str(repo),'--reference-campaign',str(repo/'source'),'--output',str(repo/'new')]
            if split:argv.append('--constant-split')
            if bits is not None:argv+=['--bits',*map(str,bits)]
            if pairs:argv+=['--pairs',*pairs]
            if steps==10 and ((bits and 6 in bits) or (pairs and case=='10cells')):argv+=['--screening-campaign',str(repo/'source')]
            argv+=['--steps',str(steps),'--case',case]
            if case!='10cells':
                (repo/'mesh').mkdir()
                for suffix in experiment.INPUT_SUFFIXES.values():
                    (repo/'mesh'/(case+'_kaskade'+suffix)).touch()
                argv+=['--mesh-dir',str(repo/'mesh')]
            with patch.object(sys,'argv',argv),patch.object(experiment.subprocess,'run',side_effect=run),contextlib.redirect_stdout(io.StringIO()):
                rc=experiment.main()
                with self.assertRaises(FileExistsError):experiment.main()
            self.assertEqual((source/'run.log').read_bytes(),original)
            env=json.loads((repo/'new/environment.json').read_text())
            self.assertEqual(env['case'],case)
            self.assertEqual(env['saved_first_step_parity_required'],case=='10cells')
            if pairs and case!='10cells':self.assertEqual(env['mixed_scaling_reference']['steps'],10)
            if steps==10 and not (failure or late_failure or missing_snapshot or changed_candidate or changed_mixed or late_mixed_failure):
                self.assertEqual(len((repo/'new/per-step-errors.csv').read_text().splitlines()),51 if (bits and 6 in bits) or pairs else 41)
            return rc,json.loads((repo/'new/summary.json').read_text()),calls

    def test_full_cost_can_lose_despite_interface_savings(self):
        rc,rows,calls=self.campaign()
        self.assertEqual(rc,0);self.assertEqual(calls,list(experiment.MODES))
        candidate=rows['full-fgmres8']
        self.assertGreater(candidate['saved_interface_logical_bytes_percent']['full-cg16'],0)
        self.assertLess(candidate['saved_accounted_logical_bytes_percent']['full-cg16'],0)
        self.assertFalse(candidate['beats_both_accounted_baselines'])

    def test_failures_are_not_savings(self):
        rc,rows,_=self.campaign(failure='full-fgmres8')
        self.assertEqual(rc,1)
        self.assertIsNone(rows['full-fgmres8']['saved_accounted_logical_bytes_percent']['full-cg16'])
        rc,_,calls=self.campaign(failure='uncompressed-fgmres')
        self.assertEqual(rc,1);self.assertEqual(calls,list(experiment.MODES[:3]))

    def test_missing_accounting_rejects_control(self):
        rc,rows,calls=self.campaign(missing_ledger=True)
        self.assertEqual(rc,1);self.assertEqual(calls,['direct','full-cg16'])
        self.assertIn('ledger',' '.join(rows['full-cg16']['errors']))

    def test_split_uses_fresh_matched_controls_and_keeps_total_cost(self):
        rc,rows,calls=self.campaign(split=True)
        self.assertEqual(rc,0)
        self.assertEqual(calls,list(experiment.MODES[:-1])+[experiment.SPLIT_MODE])
        self.assertNotIn('full-fgmres8',rows)
        self.assertTrue(all(rows[m]['saved_control_parity'] for m in experiment.MODES[:-1]))
        candidate=rows[experiment.SPLIT_MODE]
        self.assertGreater(candidate['saved_interface_logical_bytes_percent']['full-cg16'],0)
        self.assertLess(candidate['saved_accounted_logical_bytes_percent']['full-cg16'],0)
        self.assertFalse(candidate['beats_both_accounted_baselines'])
        rc,rows,calls=self.campaign(split=True,changed_control=True)
        self.assertEqual(rc,1);self.assertNotIn(experiment.SPLIT_MODE,calls)
        self.assertIn('parity',' '.join(rows['full-fgmres16']['errors']))

    def test_split_failure_has_progress_but_no_claimed_savings(self):
        rc,rows,_=self.campaign(split=True,failure=experiment.SPLIT_MODE)
        self.assertEqual(rc,1)
        row=rows[experiment.SPLIT_MODE]
        self.assertIsNone(row['saved_interface_logical_bytes_percent']['full-cg16'])
        self.assertEqual(row['residual_progress']['best_true_relative'],.5)
        self.assertEqual(row['residual_progress']['last_true_relative'],.7)
        rc,rows,_=self.campaign(split=True,missing_split_marker=True)
        self.assertEqual(rc,1)
        self.assertIn('Constant-split marker mismatch',rows[experiment.SPLIT_MODE]['errors'])

    def test_ten_steps_change_only_duration(self):
        saved=seed(Path('/repo/emiModel'))
        for mode in experiment.MODES[:-1]+(experiment.SPLIT_MODE,):
            one=experiment.make_command(saved,mode,Path('/new'),1)
            ten=experiment.make_command(saved,mode,Path('/new'),10)
            changed={one[i-1] for i,(a,b) in enumerate(zip(one,ten)) if a!=b}
            self.assertEqual(changed,{'--finalTime','--maximumNumberOfTimeSteps'})
            self.assertEqual(experiment.value(ten,'--finalTime'),'0.10')
        rc,rows,_=self.campaign(split=True,steps=10)
        self.assertEqual(rc,0)
        row=rows[experiment.SPLIT_MODE]
        self.assertTrue(row['saved_control_parity'])
        self.assertEqual(len(row['fgmres']['per_step']),10)
        self.assertEqual(row['fgmres']['applications'],300)
        self.assertEqual(row['outer_collectives'],1240)
        self.assertEqual(row['outer_logical_bytes'],10*7*(30*44+92))
        self.assertAlmostEqual(row['per_step_errors'][-1]['mean_error'],.1)
        self.assertGreater(row['full_precision_vs_direct']['relative_l2'],0)
        self.assertLess(row['full_precision_vs_direct']['centered_relative_l2'],1e-15)

    def test_later_step_failures_cannot_pass(self):
        for kwargs in ({'missing_snapshot':6},{'late_failure':True},{'changed_candidate':True}):
            rc,rows,_=self.campaign(split=True,steps=10,**kwargs)
            self.assertEqual(rc,1)
            self.assertEqual(rows[experiment.SPLIT_MODE]['status'],'FAIL')
            self.assertIsNone(rows[experiment.SPLIT_MODE]['saved_accounted_logical_bytes_percent']['full-cg16'])

    def test_variable_iterations_and_per_step_overhead(self):
        metrics={'dofs':2,'mpi_ranks':8,'iterations_per_step':[1,3,2],
                 'transmitted_sent_bytes':123}
        for mode,payload,calls in [('full-cg16',7*9*48,18),
                                   ('full-fgmres16',7*(6*44+3*92),36)]:
            line=f'BDDC outer solver logical payload bytes: {payload} collectives={calls}\n'
            result=experiment.inspect_accounting(line,mode,metrics)
            self.assertEqual(result['accounted_logical_bytes'],payload+123)
            with self.assertRaises(ValueError):
                experiment.inspect_accounting(line.replace(str(payload),str(payload-1)),mode,metrics)

    def test_scaling_cases_use_fresh_controls_without_cross_mesh_parity(self):
        for case in ('weak_8_h16','strong_h16'):
            rc,rows,calls=self.campaign(split=True,steps=10,case=case,changed_control=True)
            self.assertEqual(rc,0)
            self.assertEqual(calls,list(experiment.MODES[:-1])+[experiment.SPLIT_MODE])
            self.assertTrue(all('saved_control_parity' not in row for row in rows.values()))
            self.assertEqual(len(rows[experiment.SPLIT_MODE]['per_step_errors']),10)
        rc,rows,calls=self.campaign(split=True,steps=10,case='weak_8_h16',failure='full-cg16')
        self.assertEqual(rc,1)
        self.assertNotIn(experiment.SPLIT_MODE,calls)

    def test_scaling_commands_keep_codec_physics_and_solver(self):
        saved=seed(Path('/repo/emiModel'))
        for case in ('weak_8_h16','strong_h16'):
            for mode in experiment.MODES[:-1]+(experiment.SPLIT_MODE,):
                control=experiment.make_command(saved,mode,Path('/new'),10)
                candidate=experiment.make_command(saved,mode,Path('/new'),10,case,Path('/mesh'))
                changed={control[i-1] for i,(a,b) in enumerate(zip(control,candidate)) if a!=b}
                self.assertEqual(changed,set(experiment.INPUT_SUFFIXES)|{'--refine'})
        bad=list(saved);experiment.set_value(bad,'--sigma_i',99)
        with self.assertRaisesRegex(ValueError,'physical parameters'):
            experiment.make_command(bad,experiment.SPLIT_MODE,Path('/new'),10,'weak_8_h16',Path('/mesh'))

    def test_lower_bits_change_only_precision(self):
        saved=seed(Path('/repo/emiModel'))
        control=experiment.make_command(saved,experiment.SPLIT_MODE,Path('/new'))
        for mode,bits in experiment.SPLIT_BITS.items():
            cmd=experiment.make_command(saved,mode,Path('/new'))
            changed={control[i-1] for i,(a,b) in enumerate(zip(control,cmd)) if a!=b}
            expected=set() if bits==8 else {'--bddcCompressionBits',
                '--bddcRestrictCompressionBits','--bddcProlongateCompressionBits'}
            self.assertEqual(changed,expected)
            for flag in expected:self.assertEqual(experiment.value(cmd,flag),str(bits))
            for flag in ('--bddcProlongateConstantSplit','--bddcRestrictZeroPreserving'):
                self.assertEqual(experiment.value(cmd,flag),'1')
            self.assertEqual(experiment.value(cmd,'--bddcProlongateZeroPreserving'),'0')

    def test_extended_cap_changes_only_four_bit_limit(self):
        saved=seed(Path('/repo/emiModel'))
        for mode in experiment.MODES[:-1]+tuple(experiment.SPLIT_BITS):
            base=experiment.make_command(saved,mode,Path('/new'))
            extended=experiment.make_command(saved,mode,Path('/new'),four_bit_limit=500)
            changed={base[i-1] for i,(a,b) in enumerate(zip(base,extended)) if a!=b}
            self.assertEqual(changed,{'--bddcIterations'} if mode=='full-fgmres4-split' else set())
        one=experiment.make_command(saved,'full-fgmres6-split',Path('/new'))
        ten=experiment.make_command(saved,'full-fgmres6-split',Path('/new'),steps=10)
        self.assertEqual({one[i-1] for i,(a,b) in enumerate(zip(one,ten)) if a!=b},
                         {'--maximumNumberOfTimeSteps','--finalTime'})

    def test_directional_pairs_change_only_one_precision_argument(self):
        saved=seed(Path('/repo/emiModel'))
        control=experiment.make_command(saved,experiment.SPLIT_MODE,Path('/new'))
        for mode,(r,p) in experiment.MIXED_BITS.items():
            cmd=experiment.make_command(saved,mode,Path('/new'))
            expected='--bddcRestrictCompressionBits' if r!=8 else '--bddcProlongateCompressionBits'
            self.assertEqual({control[i-1] for i,(a,b) in enumerate(zip(control,cmd)) if a!=b},{expected})
            self.assertEqual(experiment.value(cmd,'--bddcRestrictCompressionBits'),str(r))
            self.assertEqual(experiment.value(cmd,'--bddcProlongateCompressionBits'),str(p))

    def test_r6_p8_ten_step_validation_preserves_first_step_and_checks_later_steps(self):
        mode='full-fgmres-r6-p8-split'
        rc,rows,calls=self.campaign(split=True,pairs=['6/8'],steps=10)
        self.assertEqual(rc,0);self.assertEqual(calls[-1],mode)
        self.assertTrue(rows[mode]['saved_control_parity'])
        self.assertEqual(len(rows[mode]['per_step_errors']),10)
        self.assertEqual((rows[mode]['restriction_bits'],rows[mode]['prolongation_bits']),(6,8))
        for failure in ({'changed_mixed':True},{'late_mixed_failure':True}):
            rc,rows,_=self.campaign(split=True,pairs=['6/8'],steps=10,**failure)
            self.assertEqual(rc,1);self.assertEqual(rows[mode]['status'],'FAIL')
            self.assertIsNone(rows[mode]['saved_interface_logical_bytes_percent'][experiment.SPLIT_MODE])
        saved=seed(Path('/repo/emiModel'))
        one=experiment.make_command(saved,mode,Path('/new'),steps=1)
        ten=experiment.make_command(saved,mode,Path('/new'),steps=10)
        self.assertEqual({one[i-1] for i,(a,b) in enumerate(zip(one,ten)) if a!=b},{'--finalTime','--maximumNumberOfTimeSteps'})

    def test_mixed_scaling_runs_fresh_controls_and_checks_all_steps(self):
        mode='full-fgmres-r6-p8-split'
        for case in ('weak_8_h16','strong_h16'):
            rc,rows,calls=self.campaign(split=True,pairs=['6/8'],steps=10,case=case,changed_mixed=True)
            self.assertEqual(rc,0)
            self.assertEqual(calls,list(experiment.MODES[:-1])+[experiment.SPLIT_MODE,mode])
            self.assertTrue(all('saved_control_parity' not in r for r in rows.values()))
            self.assertEqual(len(rows[mode]['per_step_errors']),10)
            self.assertEqual((rows[mode]['restriction_bits'],rows[mode]['prolongation_bits']),(6,8))
            saved=seed(Path('/repo/emiModel'))
            control=experiment.make_command(saved,experiment.SPLIT_MODE,Path('/new'),10,case,Path('/mesh'))
            mixed=experiment.make_command(saved,mode,Path('/new'),10,case,Path('/mesh'))
            self.assertEqual({control[i-1] for i,(a,b) in enumerate(zip(control,mixed)) if a!=b},{'--bddcRestrictCompressionBits'})
        rc,rows,_=self.campaign(split=True,pairs=['6/8'],steps=10,case='strong_h16',late_mixed_failure=True)
        self.assertEqual(rc,1);self.assertEqual(rows[mode]['status'],'FAIL')
        self.assertIsNone(rows[mode]['saved_interface_logical_bytes_percent'][experiment.SPLIT_MODE])
        rc,_,calls=self.campaign(split=True,pairs=['6/8'],steps=10,case='weak_8_h16',failure=experiment.SPLIT_MODE)
        self.assertEqual(rc,1);self.assertNotIn(mode,calls)

    def test_mixed_scaling_rejects_failed_or_wrong_precision_reference(self):
        for options in ({'reference_mixed_failure':True},{'reference_mixed_bits':8}):
            with self.assertRaises(ValueError):
                self.campaign(split=True,pairs=['6/8'],steps=10,case='weak_8_h16',**options)

    def test_directional_campaign_preserves_controls_and_keeps_failed_candidates(self):
        expected=list(experiment.MODES[:-1])+[experiment.SPLIT_MODE,'full-fgmres-r8-p6-split','full-fgmres-r6-p8-split']
        rc,rows,calls=self.campaign(split=True,pairs=['8/6','6/8'])
        self.assertEqual(rc,0);self.assertEqual(calls,expected)
        self.assertTrue(rows[experiment.SPLIT_MODE]['saved_control_parity'])
        for m,pair in [('full-fgmres-r8-p6-split',(8,6)),('full-fgmres-r6-p8-split',(6,8))]:
            self.assertEqual((rows[m]['restriction_bits'],rows[m]['prolongation_bits']),pair)
            self.assertTrue(rows[m]['constant_split'])
        rc,rows,calls=self.campaign(split=True,pairs=['8/6','6/8'],failure=expected[-2])
        self.assertEqual(rc,1);self.assertEqual(calls,expected)
        self.assertEqual(rows[expected[-1]]['status'],'PASS')
        self.assertIsNone(rows[expected[-2]]['saved_interface_logical_bytes_percent']['full-cg16'])
        rc,rows,calls=self.campaign(split=True,pairs=['8/6','6/8'],failure=experiment.SPLIT_MODE)
        self.assertEqual(rc,1);self.assertNotIn(expected[-2],calls)

    def test_six_bit_ten_steps_checks_first_step_and_all_ten_residuals(self):
        rc,rows,calls=self.campaign(split=True,steps=10,bits=[6])
        self.assertEqual(rc,0)
        candidate=rows['full-fgmres6-split']
        self.assertTrue(candidate['saved_control_parity'])
        self.assertEqual(len(candidate['per_step_errors']),10)
        self.assertEqual(candidate['restriction_bits'],6)
        self.assertEqual(candidate['prolongation_bits'],6)
        self.assertEqual(calls[-2:],['full-fgmres8-split','full-fgmres6-split'])

    def test_extended_trace_rejects_numerical_drift_or_shorter_history(self):
        with tempfile.TemporaryDirectory() as tmp:
            old=Path(tmp)/'old.csv';new=Path(tmp)/'new.csv'
            old.write_text('iteration,true_relative_residual\n0,1\n1,0.1\n')
            new.write_text(old.read_text()+'2,0.01\n')
            self.assertEqual(experiment.check_trace_prefix(old,new),2)
            for text in ('iteration,true_relative_residual\n0,1\n',
                         'iteration,true_relative_residual\n0,1\n1,0.11\n2,0.01\n'):
                new.write_text(text)
                with self.assertRaisesRegex(ValueError,'prefix'):experiment.check_trace_prefix(old,new)

    def test_low_bit_screen_keeps_8_control_and_runs_both_candidates(self):
        rc,rows,calls=self.campaign(split=True,bits=[6,4])
        expected=list(experiment.MODES[:-1])+list(experiment.SPLIT_BITS)
        self.assertEqual(rc,0);self.assertEqual(calls,expected)
        self.assertTrue(rows[experiment.SPLIT_MODE]['saved_control_parity'])
        for mode in ('full-fgmres6-split','full-fgmres4-split'):
            self.assertEqual(rows[mode]['restriction_bits'],experiment.SPLIT_BITS[mode])
            self.assertEqual(rows[mode]['prolongation_bits'],experiment.SPLIT_BITS[mode])
            self.assertNotIn('saved_control_parity',rows[mode])
            self.assertIn('uncompressed-fgmres',rows[mode]['saved_interface_logical_bytes_percent'])
            self.assertIn(experiment.SPLIT_MODE,rows[mode]['saved_interface_logical_bytes_percent'])
        rc,rows,calls=self.campaign(split=True,bits=[6,4],failure='full-fgmres6-split')
        self.assertEqual(rc,1);self.assertEqual(calls,expected)
        self.assertEqual(rows['full-fgmres4-split']['status'],'PASS')
        self.assertIsNone(rows['full-fgmres6-split']['saved_interface_logical_bytes_percent']['uncompressed-fgmres'])
        rc,rows,calls=self.campaign(split=True,bits=[6,4],failure=experiment.SPLIT_MODE)
        self.assertEqual(rc,1);self.assertNotIn('full-fgmres6-split',calls)
        self.assertNotIn('full-fgmres4-split',calls)

if __name__=='__main__':unittest.main()
