#!/usr/bin/env python3
"""Check that paired commands change only the requested solver/transfer settings."""
import importlib.util
from pathlib import Path
import unittest
import subprocess
import sys
import tempfile

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

    def test_rizzo_mesh_selection(self):
        self.assertIsNone(runner.reference_directory('rizzo'))
        with tempfile.TemporaryDirectory() as directory:
            source=[]
            for flag,suffix in runner.MESH_INPUT_SUFFIXES.items():
                source += [flag,'/old/pepe_combi_domi'+suffix]
                (Path(directory)/('rizzo_combi_domi'+suffix)).touch()
            for flag,value in {'--refine':0,'-np':8,'--sigma_i':3,'--sigma_e':20,
                               '--C_m':1,'--R':0.1,'--dt':0.01,'--bddcIterations':300}.items():
                source += [flag,str(value)]
            result=runner.configure_mesh_case('rizzo',source,directory)
            expected=source.copy()
            for flag,suffix in runner.MESH_INPUT_SUFFIXES.items():
                expected[expected.index(flag)+1]=str((Path(directory)/('rizzo_combi_domi'+suffix)).resolve())
            self.assertEqual(result,expected)
            self.assertEqual(runner.configure_mesh_case('pepe',source,directory),source)
            wrong=source.copy();wrong[wrong.index('--refine')+1]='2'
            with self.assertRaises(ValueError):runner.configure_mesh_case('rizzo',wrong,directory)
            with self.assertRaises(ValueError):runner.configure_mesh_case('rizzo',source+['--input','duplicate'],directory)
            (Path(directory)/'rizzo_combi_domi_excited.txt').unlink()
            with self.assertRaises(ValueError):runner.configure_mesh_case('rizzo',source,directory)

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

    def test_rizzo_cli(self):
        result=subprocess.run([sys.executable,str(runner_path),'--validate-only',
                               '--case','rizzo','--solvers','cg','fcg','fgmres',
                               '--bits','12','16'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

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

    def test_explicit_modes_add_only_needed_controls(self):
        self.assertEqual(
            runner.selected_experiment_modes([], [], 'both', ['cg12-split','cg16-split','cg24-split']),
            ['uncompressed-cg','cg12-split','cg16-split','cg24-split'])
        self.assertEqual(
            runner.selected_experiment_modes([], [], 'off', ['cg12-split','fcg16','fgmres24-split']),
            ['uncompressed-cg','uncompressed-fcg','uncompressed-fgmres',
             'cg12-split','fcg16','fgmres24-split'])
        self.assertEqual(
            runner.selected_experiment_modes([], [], 'both', ['fcg12','uncompressed-fcg','fcg16-split']),
            ['uncompressed-fcg','fcg12','fcg16-split'])
        self.assertEqual(runner.selected_experiment_modes([], [], 'off', ['uncompressed-cg']),
                         ['uncompressed-cg'])

    def test_explicit_mode_solver_flags(self):
        for mode, flexible, fgmres in [('cg12-split',0,0),('fcg16',1,0),('fgmres24-split',0,1)]:
            cmd=runner.make_command(self.source,'fresh',mode)
            self.assertEqual(cmd[cmd.index('--bddcFlexibleCg')+1],str(flexible))
            self.assertEqual(cmd[cmd.index('--bddcFgmres')+1],str(fgmres))
            self.assertEqual(cmd[cmd.index('--bddcProlongateConstantSplit')+1],str(int(mode.endswith('-split'))))

    def test_explicit_modes_reject_invalid_and_duplicates(self):
        for modes in ([], ['cg12','cg12'], ['cg2'], ['gmres16'], ['uncompressed-cg-split']):
            with self.subTest(modes=modes), self.assertRaises(ValueError):
                runner.selected_experiment_modes([], [], 'both', modes)

    def test_legacy_selection_unchanged(self):
        for split in ('off','on','both'):
            self.assertEqual(runner.selected_experiment_modes(['cg','fcg'],[12,16],split),
                             runner.experiment_modes(['cg','fcg'],[12,16],split))

    def test_explicit_modes_cli_rank_overrides(self):
        for case in ('10cells','pepe','rizzo','robin','weak_8_h16','strong_h16'):
            cmd=[sys.executable,str(runner_path),'--validate-only','--case',case,
                 '--mpi-ranks','2','--modes','cg12-split','cg16-split','cg24-split']
            result=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('uncompressed-cg, cg12-split, cg16-split, cg24-split',result.stdout)
        base=[sys.executable,str(runner_path),'--validate-only','--modes','cg12-split']
        for extra in (['--bits','12'],['--solvers','cg'],['--split','on'],['--mpi-ranks','0']):
            result=subprocess.run(base+extra,capture_output=True,text=True)
            self.assertEqual(result.returncode,2,result.stdout)

if __name__=='__main__':unittest.main()
