import shlex
import unittest
from analyze_emi_potential_diagnostics import error_metrics

class TrajectoryChecks(unittest.TestCase):
    def test_synthetic_offset_and_maximum(self):
        m=error_metrics([2.,3.],[1.,2.])
        self.assertAlmostEqual(m['relative_l2'],(2/5)**.5)
        self.assertEqual(m['max_absolute_error'],1.)
        self.assertEqual(m['mean_error'],1.)
        self.assertEqual(m['centered_relative_l2'],0.)
    def test_identical_and_mismatched_states(self):
        self.assertEqual(error_metrics([1.,2.],[1.,2.])['absolute_l2'],0.)
        with self.assertRaises(ValueError): error_metrics([1.],[1.,2.])
    def test_generated_commands(self):
        from pathlib import Path
        from run_emi_solver_study import arguments, make_command
        a=arguments(['--case','pepe','--solvers','fcg','--split','on','--bits','16',
                     '--slurm-container','/tmp/image.sif','--kaskade','/tmp/kaskade'])
        a.steps=100
        commands=[]
        for bits in (None,16):
            case=dict(solver='fcg',bits=bits,split='on' if bits else 'off',direct=False)
            cmd=make_command(a,case,Path('/tmp/output'),['mpirun'])
            cmd[cmd.index('--vtk')+1]='1'
            commands.append(cmd)
        self.assertEqual(len(commands),2)
        for i,c in enumerate(commands):
            for flag,value in [('--maximumNumberOfTimeSteps','100'),('--finalTime','1'),
                               ('--dt','0.01'),('--vtk','1'),('--potentialDiagnostics','1'),
                               ('--bddcFlexibleCg','1'),('--bddcCompression',str(i)),
                               ('--bddcProlongateConstantSplit',str(i))]:
                self.assertEqual(c[c.index(flag)+1],value)
            self.assertIn('--ntasks=16',c)
            self.assertIn('--nodes=8',c)
if __name__=='__main__': unittest.main()
