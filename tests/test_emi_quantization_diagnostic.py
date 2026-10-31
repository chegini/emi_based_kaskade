import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class QuantizationDiagnostic(unittest.TestCase):
    def test_reference_order_uniform_bits_and_failed_sweep_continuation(self):
        runner = Path(__file__).with_name('run_emi_quantization_diagnostic.py')
        with tempfile.TemporaryDirectory(prefix='emi-quantization-control-') as folder:
            repo = Path(folder)
            (repo/'input').mkdir()
            for suffix in ['_mesh.vtu','_list_extracellular.txt','_list_intracellular.txt','_early_excited.txt']:
                (repo/'input'/('10Cells3d_10extra'+suffix)).touch()
            model = repo/'emiModel'
            model.write_text('#!/bin/sh\nexit 0\n')
            model.chmod(0o755)
            fake = repo/'fake-mpi'
            fake.write_text('''#!/usr/bin/env python3
import csv,json,os,sys
from pathlib import Path
args=sys.argv[1:]
def value(key):return args[args.index(key)+1]
with open('capture.jsonl','a') as f:f.write(json.dumps(args)+'\\n')
out=Path(value('--dir'))
if value('--direct')=='1':
    print('Direct global residual step 1: abs=1e-18 rel=4e-16')
    (out/'emiLast.vtu').touch()
    if os.environ.get('FAIL_DIRECT')=='1':sys.exit(17)
else:
    compressed=value('--bddcCompression')=='1'
    full=value('--bddcHuffman')=='1'
    bad=compressed and value('--bddcCompressionBits')=='8' and not full
    rel=0.2 if bad else 7e-8
    print('cells: 2560\\ndofs: 1156')
    print('BDDC subdomains: 20, shared dofs: 216')
    print('BDDC coarse interfaces: corners=0, edges=27, faces=28')
    print(f"BDDC global residual step 1: abs=1e-10 rel={rel} target={value('--bddcRelativeTolerance')} status={'not-converged' if bad else 'converged'} iterations={300 if bad else 10}")
    with (out/'bddc-convergence-step1.csv').open('w') as f:
        w=csv.writer(f)
        w.writerow(['iteration','relative_residual'])
        w.writerow([300 if bad else 10,rel])
    if bad:
        print('BDDC true relative residual target not reached: iteration limit')
        sys.exit(134)
    print('Compression report (BDDC):')
    print('  raw sent bytes: '+str(30000 if full else 10000))
    print('  transmitted sent bytes: '+str(12000 if full else 10000))
    print(f"Profile (wall seconds, EMI end-to-end phases; {value('-np')} rank(s)):")
    print('  solver wrapper (including setup) 1.0 1.1 1.2')
    (out/('emiBDDCLastCompression.vtu' if compressed else 'emiBDDCLast.vtu')).touch()
print('total cpu-time:  2.300000s wall, 3s user + 1s system = 4s CPU')
print('EMI model completed')
''')
            fake.chmod(0o755)
            (repo/'compare_emi_vtu.py').write_text('''import sys
from pathlib import Path
assert all(Path(p).is_file() for p in sys.argv[1:])
print('same coordinates: yes')
print('relative L2 error: 2e-9')
print('max absolute error: 1e-8')
print('mean signed offset (first - reference): 1e-10')
print('relative L2 after removing mean offset: 1e-9')
''')
            env = dict(os.environ, MPIEXEC=str(fake),RESTRICT_BITS='32',PROLONGATE_BITS='8')
            command = [sys.executable,str(runner),'--repo',str(repo),'--output','sweep']
            result = subprocess.run(command,cwd=repo,env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            calls = [json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()]
            self.assertEqual(len(calls),7)
            self.assertTrue(all(row[row.index('--refine')+1]=='0' for row in calls))
            self.assertEqual(calls[0][calls[0].index('--direct')+1],'1')
            self.assertEqual(calls[0][calls[0].index('-np')+1],'1')
            self.assertEqual(calls[1][calls[1].index('--bddcCompression')+1],'0')
            for row,bits in zip(calls[2:],[8,12,16,24,32]):
                def value(key):return row[row.index(key)+1]
                self.assertEqual(value('--bddcUseCg'),'1')
                self.assertEqual(value('--bddcRestrictionPasses'),'1')
                for flag in ['--bddcCompressionBits','--bddcRestrictCompressionBits','--bddcProlongateCompressionBits']:
                    self.assertEqual(value(flag),str(bits))
                for flag in ['--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                             '--bddcZigzag','--bddcBitlength','--bddcHuffman']:
                    self.assertEqual(value(flag),'0')
            summary = json.loads((repo/'sweep/summary.json').read_text())
            self.assertEqual([r['status'] for r in summary],['PASS','PASS','FAIL','PASS','PASS','PASS','PASS'])
            self.assertEqual(summary[2]['relative_residual'],0.2)
            self.assertEqual(summary[2]['reason'],'iteration limit')
            self.assertNotIn('comparison_to_direct',summary[2])
            self.assertIn('comparison_to_direct',summary[3])
            self.assertIn('comparison_to_uncompressed',summary[3])
            self.assertNotIn('Missing normal completion marker',result.stderr)
            refined = subprocess.run(command[:-1]+['refined','--bits','8','--restriction-passes','4'],
                                     cwd=repo,env=env,capture_output=True,text=True)
            self.assertEqual(refined.returncode,1)
            refined_summary=json.loads((repo/'refined/summary.json').read_text())
            self.assertEqual(len(refined_summary),3)
            self.assertEqual(refined_summary[1]['restriction_passes'],1)
            self.assertEqual(refined_summary[2]['restriction_passes'],4)
            self.assertEqual(refined_summary[2]['restriction_bits'],8)
            self.assertEqual(refined_summary[2]['prolongation_bits'],8)
            captured=[json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()]
            last=captured[-1]
            self.assertEqual(last[last.index('--bddcRestrictionPasses')+1],'4')
            self.assertEqual(last[last.index('--bddcRestrictCompressionBits')+1],'8')
            self.assertEqual(last[last.index('--bddcProlongateCompressionBits')+1],'8')
            full_result = subprocess.run(
                command[:-1]+['full','--bits','8','--restriction-passes','4',
                              '--full-compression','--ranks','8','--refine','3'],
                cwd=repo,env=env,capture_output=True,text=True)
            self.assertEqual(full_result.returncode,1,full_result.stdout+full_result.stderr)
            full_summary=json.loads((repo/'full/summary.json').read_text())
            self.assertEqual([r['case'] for r in full_summary],
                             ['direct','uncompressed-cg','quantization-8','full-8'])
            full_row=full_summary[-1]
            self.assertEqual(full_row['status'],'PASS')
            self.assertEqual(full_row['ranks'],8)
            self.assertEqual(full_summary[0]['ranks'],1)
            self.assertAlmostEqual(full_row['net_logical_byte_savings_vs_uncompressed_percent'],-20)
            self.assertAlmostEqual(full_row['solver_metrics']['reported_transfer_savings_percent'],60)
            captured=[json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()]
            self.assertTrue(all(row[row.index('--refine')+1]=='3' for row in captured[-4:]))
            full_environment=json.loads((repo/'full/environment.json').read_text())
            self.assertEqual(full_environment['mesh_refinements'],3)
            last=captured[-1]
            for flag in ['--bddcGraphLifting','--bddcGraphLiftingAfterQuantization',
                         '--bddcZigzag','--bddcBitlength','--bddcHuffman']:
                self.assertEqual(last[last.index(flag)+1],'1')
            self.assertEqual(last[last.index('--bddcRestrictionPasses')+1],'4')
            self.assertIn('comparison_to_direct',full_row)
            self.assertIn('comparison_to_uncompressed',full_row)
            invalid = subprocess.run(command+['--refine','-1'],cwd=repo,env=env,
                                     capture_output=True,text=True)
            self.assertEqual(invalid.returncode,2)
            self.assertIn('--refine must be nonnegative',invalid.stderr)
            self.assertEqual((repo/'capture.jsonl').read_text().splitlines(),
                             [json.dumps(row) for row in captured])
            stopped = subprocess.run(command[:-1]+['direct-failure'],cwd=repo,
                                     env=dict(env,FAIL_DIRECT='1'),capture_output=True,text=True)
            self.assertEqual(stopped.returncode,1)
            stopped_summary=json.loads((repo/'direct-failure/summary.json').read_text())
            self.assertEqual(len(stopped_summary),1)
            self.assertFalse((repo/'direct-failure/uncompressed-cg').exists())


if __name__=='__main__':
    unittest.main()
