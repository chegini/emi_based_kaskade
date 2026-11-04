#!/usr/bin/env python3
"""Check that paired commands change only the requested solver/transfer settings."""
import importlib.util
from pathlib import Path
import unittest
import subprocess
import sys

runner_path=Path(__file__).with_name('run_emi_cg_split_check.py')
spec=importlib.util.spec_from_file_location('cg_split_runner',runner_path)
runner=importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

class CommandChecks(unittest.TestCase):
    def setUp(self):
        names=('dir bddcHistogram bddcUseCg bddcFlexibleCg bddcFgmres '
               'bddcCgOriginalResidual bddcConvergenceReport bddcProlongateConstantSplit '
               'bddcCompression bddcGraphLifting bddcGraphLiftingAfterQuantization '
               'bddcZigzag bddcBitlength bddcHuffman bddcCompressionBits '
               'bddcRestrictCompressionBits bddcProlongateCompressionBits').split()
        self.source=['emiModel']+[item for name in names for item in ('--'+name,'0')]
        self.source+=['--bddcRestrictionPasses','1','--bddcRelativeTolerance','1e-7']
    def test_scaling_case_selection(self):
        for case,stamp in [('weak_8_h16','20261101T210711Z-16464'),
                           ('strong_h16','20261101T211157Z-16694')]:
            folder=runner.reference_directory(case)
            self.assertEqual(folder.name,case)
            self.assertEqual(folder.parent.name,'emi-fgmres-'+stamp)
            self.assertEqual(runner.reference_directory(case,Path('/custom')).parent,Path('/custom'))
            cmd=[]
            for flag,suffix in {'--input':'_mesh.vtu','--extra_set':'_list_extracellular.txt',
                                '--intra_set':'_list_intracellular.txt','--excited':'_early_excited.txt'}.items():
                cmd += [flag,'/meshes/'+case+'_kaskade'+suffix]
            for flag,val in {'--refine':0,'-np':8,'--sigma_i':3,'--sigma_e':20,
                             '--C_m':1,'--R':0.1,'--dt':0.01}.items():
                cmd += [flag,str(val)]
            runner.validate_scaling_source(case,cmd)
            wrong=cmd.copy();wrong[wrong.index('--input')+1]='/meshes/wrong.vtu'
            with self.assertRaises(ValueError):runner.validate_scaling_source(case,wrong)
            wrong=cmd.copy();wrong[wrong.index('--refine')+1]='2'
            with self.assertRaises(ValueError):runner.validate_scaling_source(case,wrong)

    def test_pairs(self):
        for solver in ('cg','fcg','fgmres'):
            for bits in (4,6,8,12,16,20,24,32):
                mode=solver+str(bits)
                off=runner.make_command(self.source,'fresh',mode)
                on=runner.make_command(self.source,'fresh',mode+'-split')
                expected=off.copy()
                expected[expected.index('--bddcProlongateConstantSplit')+1]='1'
                if solver!='fgmres':expected+=['--bddcCgConstantSplitExperiment','1']
                self.assertEqual(on,expected)
                self.assertEqual(on[on.index('--bddcRestrictCompressionBits')+1],str(bits))
                self.assertEqual(on[on.index('--bddcProlongateCompressionBits')+1],str(bits))
                self.assertEqual(on[on.index('--bddcRestrictionPasses')+1],'1')
                self.assertEqual(on[on.index('--bddcRelativeTolerance')+1],'1e-7')
    def test_quantization_only(self):
        disabled=('--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                  '--bddcZigzag','--bddcBitlength','--bddcHuffman')
        for solver in ('cg','fcg','fgmres'):
            for bits in (4,6,8,12,16,20,24,32):
                for suffix in ('','-split'):
                    mode=solver+str(bits)+suffix
                    full=runner.make_command(self.source,'fresh',mode)
                    quant=runner.make_command(self.source,'fresh',mode,'quantization-only')
                    expected=full.copy()
                    for flag in disabled:expected[expected.index(flag)+1]='0'
                    self.assertEqual(quant,expected)
                    self.assertEqual(quant[quant.index('--bddcCompression')+1],'1')
            self.assertEqual(runner.make_command(self.source,'fresh','uncompressed-'+solver),
                             runner.make_command(self.source,'fresh','uncompressed-'+solver,'quantization-only'))
        with self.assertRaises(ValueError):runner.make_command(self.source,'fresh','cg16','unknown')

    def test_mode_selection(self):
        for split,count in [('off',27),('on',27),('both',51)]:
            modes=runner.experiment_modes(('cg','fcg','fgmres'),(4,6,8,12,16,20,24,32),split)
            self.assertEqual(len(modes),count)
            self.assertEqual(len(set(modes)),count)
            self.assertEqual(modes[:3],['uncompressed-cg','uncompressed-fcg','uncompressed-fgmres'])
            if split=='off':self.assertFalse(any(m.endswith('-split') for m in modes))
            if split=='on':self.assertTrue(all(m.endswith('-split') for m in modes[3:]))

    def test_cli_preflight(self):
        base=[sys.executable,str(runner_path),'--validate-only','--case','strong_h16',
              '--solvers','cg','fcg','fgmres']
        for bits,expected in [(['4'],0),(['6'],0),(['8'],0),(['4','6','8','12','16'],0),(['3'],2),(['4','4'],2)]:
            result=subprocess.run(base+['--bits']+bits,capture_output=True,text=True)
            self.assertEqual(result.returncode,expected,result.stderr)

    def test_controls(self):
        for solver in ('cg','fcg','fgmres'):
            cmd=runner.make_command(self.source,'fresh','uncompressed-'+solver)
            self.assertNotIn('--bddcCgConstantSplitExperiment',cmd)
            self.assertEqual(cmd[cmd.index('--bddcCompression')+1],'0')
            self.assertEqual(cmd[cmd.index('--bddcProlongateConstantSplit')+1],'0')
            self.assertEqual(cmd[cmd.index('--bddcFlexibleCg')+1],str(int(solver=='fcg')))
    def test_rejected_sources_and_modes(self):
        with self.assertRaises(ValueError):
            runner.make_command(self.source+['--bddcCgConstantSplitExperiment','1'],'fresh','cg16')
        for mode in ('uncompressed-cg-split','cg3-split','fcg33','unknown'):
            with self.assertRaises(ValueError):runner.make_command(self.source,'fresh',mode)

if __name__=='__main__':unittest.main()
