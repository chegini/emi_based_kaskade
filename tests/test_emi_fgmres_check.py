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

    def campaign(self,failure=None,missing_ledger=False,split=False,changed_control=False,missing_split_marker=False,steps=1,missing_snapshot=None,late_failure=False,changed_candidate=False):
        with tempfile.TemporaryDirectory() as temp:
            repo=Path(temp).resolve();model=repo/'emiModel';model.write_text('mock')
            source=repo/'source/10cells/full-cg16';source.mkdir(parents=True)
            (source/'command.sh').write_text(shlex.join(seed(model)))
            (source/'run.log').write_text(log(1 if split else 10))
            if split:
                cmd=seed(model);experiment.set_value(cmd,'--maximumNumberOfTimeSteps',1)
                (source/'command.sh').write_text(shlex.join(cmd))
                for mode in (experiment.MODES[:-1]+(experiment.SPLIT_MODE,) if steps==10 else experiment.MODES[:-1]):
                    folder=repo/'source/10cells'/mode;folder.mkdir(exist_ok=True)
                    k=30 if mode==experiment.SPLIT_MODE else 21
                    text=('Direct global residual step 1: abs=1e-16 rel=1e-15\nEMI model completed\n'
                          if mode=='direct' else log(calls=k))
                    (folder/'run.log').write_text(text)
                    (folder/'potential-state-step0.csv').write_text('dof,u\n0,0\n1,0\n')
                    (folder/'potential-state-step1.csv').write_text('dof,u\n0,1\n1,2\n')
                    if 'fgmres' in mode:
                        (folder/'bddc-fgmres-step1.csv').write_text('iteration,cycle,estimated_absolute_residual,true_absolute_residual,true_relative_residual\n'+
                            ''.join(f'{i},1,1e-8,1e-8,1e-8\n' for i in range(k+1)))
            original=(source/'run.log').read_bytes();calls=[]
            def run(cmd,**kw):
                out=Path(experiment.value(cmd,'--dir'));mode=out.name;calls.append(mode)
                fgmres='fgmres' in mode;candidate=mode in ('full-fgmres8',experiment.SPLIT_MODE);k=30 if candidate else 21
                text=experiment.MARKER+'\n' if fgmres else ''
                if mode==experiment.SPLIT_MODE and not missing_split_marker:text+=experiment.SPLIT_MARKER+'\n'
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
                    if candidate and late_failure:
                        text=text.replace('step 6: abs=1e-10 rel=7e-8 target=1e-7 status=converged',
                                          'step 6: abs=1e-10 rel=0.2 target=1e-7 status=not-converged')
                kw['stdout'].write(text)
                (out/'potential-state-step0.csv').write_text('dof,u\n0,0\n1,0\n')
                for step in range(1,steps+1):
                    if candidate and step==missing_snapshot:continue
                    delta=.01*step if candidate and step>1 else 0
                    (out/f'potential-state-step{step}.csv').write_text(f'dof,u\n0,{1+delta}\n1,{2+delta}\n')
                if changed_candidate and candidate:
                    (out/'potential-state-step1.csv').write_text('dof,u\n0,1.1\n1,2\n')
                if changed_control and mode=='full-fgmres16':
                    (out/'potential-state-step1.csv').write_text('dof,u\n0,1.1\n1,2\n')
                return SimpleNamespace(returncode=134 if mode==failure else 0)
            argv=['runner','--repo',str(repo),'--reference-campaign',str(repo/'source'),'--output',str(repo/'new')]
            if split:argv.append('--constant-split')
            argv+=['--steps',str(steps)]
            with patch.object(sys,'argv',argv),patch.object(experiment.subprocess,'run',side_effect=run),contextlib.redirect_stdout(io.StringIO()):
                rc=experiment.main()
                with self.assertRaises(FileExistsError):experiment.main()
            self.assertEqual((source/'run.log').read_bytes(),original)
            if steps==10 and not (failure or late_failure or missing_snapshot or changed_candidate):
                self.assertEqual(len((repo/'new/per-step-errors.csv').read_text().splitlines()),41)
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

if __name__=='__main__':unittest.main()
