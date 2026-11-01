import csv
import json
import tempfile
import unittest
from pathlib import Path
from analyze_emi_potential_diagnostics import analyze,error_metrics
from run_emi_mixed_precision_validation import command,matched_command


def write(path,header,rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)


class PotentialDiagnostics(unittest.TestCase):
    def test_error_decomposition_and_commands(self):
        e=error_metrics([1.25,2.25],[1.,2.])
        self.assertEqual(e['mean_error'],.25)
        self.assertEqual(e['centered_relative_l2'],0.)
        self.assertAlmostEqual(e['mean_fraction_of_squared_error'],1.)
        for mode in ['direct','uncompressed-cg','full-r16-p16','full-zero-r14-p16']:
            args=(Path('/model'),['mpirun'],Path('/mesh'),0,10,8,mode,Path('/out'))
            a=command(*args,candidate_bits=14)
            b=command(*args,candidate_bits=14,potential_diagnostics=True)
            self.assertTrue(matched_command(a,b))
            i=b.index('--potentialDiagnostics')
            self.assertEqual(b[i+1],'1')
            self.assertEqual(a,b[:i]+b[i+2:])

    def test_per_step_constant_error_and_balance(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);runs={'direct':{},'full-zero-r14-p16':{}}
            (root/'environment.json').write_text(json.dumps({'steps':2,'physics':{'dt':.01}}))
            (root/'summary.json').write_text(json.dumps({'mesh':{'runs':runs}}))
            for mode in runs:
                folder=root/'mesh'/mode;folder.mkdir(parents=True)
                write(folder/'potential-operator.csv',
                      ['dof','a_one','a_transpose_one','diagonal'],[[0,1,1,2],[1,.1,.1,1.1]])
                for step in range(3):
                    c=.25 if mode!='direct' and step else 0.
                    u=[step*(1+c),step*(2+c)]
                    du=[1+c,2+c] if step else [0,0]
                    rhs=[0,1.2] if step else [0,0]
                    r=[-c,-.1*c]
                    write(folder/f'potential-state-step{step}.csv',
                          ['dof','u','correction','rhs','residual'],
                          [[i,u[i],du[i],rhs[i],r[i]] for i in range(2)])
            report=analyze(root)
            self.assertTrue(report['complete'])
            with (root/'potential-analysis/per-step-errors.csv').open() as f:
                rows=list(csv.DictReader(f))
            last=rows[-1]
            self.assertAlmostEqual(float(last['mean_error']),.5)
            self.assertAlmostEqual(float(last['mean_error_increment']),.25)
            self.assertAlmostEqual(float(last['mean_correction_error']),.25)
            self.assertAlmostEqual(float(last['constant_balance_closure']),0,places=14)
            self.assertEqual(float(last['centered_relative_l2']),0)
            with self.assertRaises(FileExistsError): analyze(root)


if __name__=='__main__': unittest.main()
