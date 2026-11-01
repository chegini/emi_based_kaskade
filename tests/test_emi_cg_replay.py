from pathlib import Path
import unittest
from replay_emi_cg_diagnostics import observed_command


class Replay(unittest.TestCase):
    def test_only_observer_and_output_change(self):
        old=['mpirun','-np','8','./emiModel','--bddc','1','--direct','0',
             '--bddcUseCg','1','--bddcConvergenceReport','1','--bddcCgDiagnostics','0',
             '--dir','old','--maximumNumberOfTimeSteps','10',
             '--bddcRelativeTolerance','1e-7','--bddcRestrictCompressionBits','12',
             '--bddcProlongateCompressionBits','16','--bddcRestrictionPasses','1']
        saved=list(old)
        new=observed_command(old,Path('new'))
        self.assertEqual(old,saved)
        changes=[i for i,(a,b) in enumerate(zip(old,new)) if a!=b]
        self.assertEqual(changes,[old.index('--bddcCgDiagnostics')+1,old.index('--dir')+1])
        self.assertEqual(new[old.index('--bddcCgDiagnostics')+1],'1')
        old[old.index('--bddcUseCg')+1]='0'
        with self.assertRaises(ValueError): observed_command(old,Path('new'))


if __name__=='__main__': unittest.main()
