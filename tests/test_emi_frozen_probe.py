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

import run_emi_frozen_preconditioner_probe as probe
from test_emi_grouped_prolongation_check import seed, log


class FrozenProbeTest(unittest.TestCase):
    def test_commands_preserve_physics_and_isolate_directions(self):
        saved=seed(Path('/repo/emiModel'))
        for mode in probe.MODES+('r8-p8-split',):
            cmd=probe.make_command(saved,mode,Path('/new'))
            for flag in ('--input','--refine','--sigma_i','--sigma_e','--C_m','--R','--dt',
                         '--extra_set','--intra_set','--excited'):
                self.assertEqual(probe.value(cmd,flag),probe.value(saved,flag))
            for flag,wanted in (('--bddcPreconditionerProbe','1'),('--bddcFgmres','1'),
                                ('--bddcUseCg','0'),('--bddcRestrictionPasses','1'),
                                ('--bddcRelativeTolerance','1e-7'),('--bddcProlongateGroupSize','0'),
                                ('--bddcProlongateZeroPreserving','0'),('--bddcConstantCorrection','0'),
                                ('--maximumNumberOfTimeSteps','1'),('--vtk','0')):
                self.assertEqual(probe.value(cmd,flag),wanted)
            r,p=(8,8) if mode=='r8-p8-split' else ((16,16) if mode=='uncompressed' else tuple(int(x[1:]) for x in mode.split('-')))
            self.assertEqual(probe.value(cmd,'--bddcRestrictCompressionBits'),str(r))
            self.assertEqual(probe.value(cmd,'--bddcProlongateCompressionBits'),str(p))
            self.assertEqual(probe.value(cmd,'--bddcRestrictZeroPreserving'),str(int(r==8)))
            self.assertEqual(probe.value(cmd,'--bddcProlongateConstantSplit'),str(int(mode=='r8-p8-split')))

    def test_added_flags_are_registered_by_application(self):
        source=(Path(__file__).resolve().parents[1]/'EMI_model.cpp').read_text()
        registered=set(re.findall(r'\(\s*"([^"\n]+)"\s*,\s*options\.',source))
        saved=seed(Path('/repo/emiModel'))
        for mode in probe.MODES+('r8-p8-split',):
            cmd=probe.make_command(saved,mode,Path('/new'))
            self.assertNotIn('--writeVTK',cmd)
            for flag in cmd:
                if flag.startswith('--') and flag not in saved:
                    self.assertIn(flag[2:],registered)
        self.assertIn('vtk',registered)

    def test_application_error_is_visible_and_deduplicated(self):
        line="  what():  unrecognised option '--writeVTK'\n"
        message=probe.application_failure(line*8,134)
        self.assertIn('exit 134',message)
        self.assertEqual(message.count("unrecognised option '--writeVTK'"),1)

    def test_small_correction_error_can_have_large_operator_error(self):
        base=dict(rhs=[1,0],v=[1,0],z=[1,0],Az=[1,0])
        data=dict(base,z=[1,0.001],Az=[1,10])
        result=probe.compare(data,base)
        self.assertAlmostEqual(result['z_relative_error'],.001)
        self.assertAlmostEqual(result['Az_relative_error'],10)
        self.assertGreater(result['best_one_direction_relative_residual'],.99)
        with self.assertRaises(ValueError):probe.compare(dict(data,v=[0,1]),base)

    def test_csv_rejects_nonfinite_and_wrong_normalization(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'probe.csv'
            for payload in ('0,1,1,nan,1\n','0,1,2,1,1\n','1,1,1,1,1\n'):
                p.write_text('dof,rhs,v,z,Az\n'+payload)
                with self.assertRaises(ValueError):probe.read_probe(p)

    def campaign(self,mismatch=False,split=False,changed_control=False):
        with tempfile.TemporaryDirectory() as temp:
            repo=Path(temp).resolve();model=repo/'emiModel';model.write_text('mock')
            source=repo/'source/10cells/full-cg16';source.mkdir(parents=True)
            (source/'command.sh').write_text(shlex.join(seed(model)))
            (source/'run.log').write_text(log(10));original=(source/'run.log').read_bytes()
            csv='dof,rhs,v,z,Az\n0,3,0.6,1,2\n1,4,0.8,2,3\n'
            for mode in probe.MODES:
                old=repo/'old-probe'/mode;old.mkdir(parents=True)
                (old/'bddc-frozen-probe.csv').write_text(csv)
            calls=[]
            def run(cmd,**kw):
                folder=Path(probe.value(cmd,'--dir'));calls.append(folder.name)
                kw['stdout'].write(probe.MARKER+'\n')
                if folder.name!='uncompressed':
                    header=8 if folder.name=='r8-p8-split' else 0
                    kw['stdout'].write(f'BDDC frozen probe traffic: restriction=80 prolongation={120+header} constant_header={header}\n')
                rhs='6,0.6' if mismatch and folder.name=='r16-p16' else '3,0.6'
                (folder/'bddc-frozen-probe.csv').write_text('dof,rhs,v,z,Az\n0,'+rhs+',1,2\n1,4,0.8,2,3\n')
                if changed_control and folder.name=='r16-p16':
                    (folder/'bddc-frozen-probe.csv').write_text(csv.replace('0.6,1,2','0.6,1,2.1'))
                return SimpleNamespace(returncode=0)
            argv=['probe','--repo',str(repo),'--reference-campaign',str(repo/'source'),'--output',str(repo/'new')]
            if split:argv+=['--constant-split','--reference-probe',str(repo/'old-probe')]
            with patch.object(sys,'argv',argv),patch.object(probe.subprocess,'run',side_effect=run),contextlib.redirect_stdout(io.StringIO()):
                rc=probe.main()
                with self.assertRaises(FileExistsError):probe.main()
            self.assertEqual((source/'run.log').read_bytes(),original)
            return rc,json.loads((repo/'new/summary.json').read_text()),calls

    def test_independent_jobs_and_preserved_results(self):
        rc,rows,calls=self.campaign()
        self.assertEqual(rc,0);self.assertEqual(calls,list(probe.MODES))
        self.assertTrue(all(r['status']=='MEASURED' for r in rows.values()))
        rc,rows,calls=self.campaign(mismatch=True)
        self.assertEqual(rc,1);self.assertEqual(calls,['uncompressed','r16-p16'])
        self.assertEqual(rows['r16-p16']['status'],'FAIL')

    def test_split_keeps_controls_and_counts_header(self):
        rc,rows,calls=self.campaign(split=True)
        self.assertEqual(rc,0);self.assertEqual(calls,list(probe.MODES)+['r8-p8-split'])
        self.assertTrue(all(rows[m]['saved_control_parity'] for m in probe.MODES))
        self.assertEqual(rows['r8-p8-split']['interface_bytes'],208)
        self.assertEqual(rows['r8-p8-split']['constant_header_bytes'],8)
        rc,rows,calls=self.campaign(split=True,changed_control=True)
        self.assertEqual(rc,1);self.assertEqual(calls,['uncompressed','r16-p16'])
        self.assertIn('parity',' '.join(rows['r16-p16']['errors']))

    def test_header_accounting_is_required(self):
        for text in ('', 'BDDC frozen probe traffic: restriction=1 prolongation=7 constant_header=8\n',
                     'BDDC frozen probe traffic: restriction=1 prolongation=80 constant_header=0\n'):
            with self.assertRaises(ValueError):probe.read_traffic(text,True)


if __name__=='__main__':unittest.main()
