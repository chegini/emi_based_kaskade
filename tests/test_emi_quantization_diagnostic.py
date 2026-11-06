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
with open('range-env.jsonl','a') as f:f.write(json.dumps(os.environ.get('KASKADE_RESTRICTION_CAPTURE'))+'\\n')
out=Path(value('--dir'))
if value('--direct')=='1':
    print('Direct global residual step 1: abs=1e-18 rel=4e-16')
    (out/'emiLast.vtu').touch()
    if os.environ.get('FAIL_DIRECT')=='1':sys.exit(17)
else:
    compressed=value('--bddcCompression')=='1'
    full=value('--bddcHuffman')=='1'
    bad=compressed and value('--bddcCompressionBits')=='8' and not full
    family='fgmres' if value('--bddcFgmres')=='1' else 'fcg' if value('--bddcFlexibleCg')=='1' else 'cg'
    if value('--bddcGraphLiftingBothOrders')=='1':
        assert value('--bddcGraphLifting')==value('--bddcGraphLiftingAfterQuantization')=='1'
        print('BDDC combined lifting: floating transform, quantization, reversible integer transform')
    if os.environ.get('FAMILY_TEST'):
        assert value('--bddcUseCg')==str(int(family!='fgmres'))
        assert value('--bddcConvergenceReport')==str(int(family!='fgmres'))
        assert value('--bddcCgOriginalResidual')==str(int(family!='fgmres'))
        if family=='fgmres':
            assert '--bddcCgConstantSplitExperiment' not in args
            assert value('--bddcFgmresRestart')=='30'
        elif value('--bddcProlongateConstantSplit')=='1':
            assert value('--bddcCgConstantSplitExperiment')=='1'
        bad=compressed and family=='cg' and value('--bddcCompressionBits')=='12' and value('--bddcGraphLifting')=='1' and value('--bddcProlongateConstantSplit')=='0'
        if os.environ.get('FAIL_FCG_BASELINE') and family=='fcg' and not compressed:
            sys.exit(19)

    rel=0.2 if bad else 7e-8
    print('cells: 2560\\ndofs: 1156')
    print('BDDC subdomains: 20, shared dofs: 216')
    print('BDDC coarse interfaces: corners=0, edges=27, faces=28')
    print(f"BDDC global residual step 1: abs=1e-10 rel={rel} target={value('--bddcRelativeTolerance')} status={'not-converged' if bad else 'converged'} iterations={300 if bad else 10}")
    with (out/('bddc-fgmres-step1.csv' if family=='fgmres' else 'bddc-convergence-step1.csv')).open('w') as f:
        w=csv.writer(f)
        w.writerow(['iteration','relative_residual'])
        w.writerow([300 if bad else 10,rel])
    if bad:
        print('BDDC true relative residual target not reached: iteration limit')
        sys.exit(134)
    print('Compression report (BDDC):')
    print('  raw sent bytes: '+str(30000 if full else 10000))
    transmitted=(12000 if full else 10000)
    if os.environ.get('FAMILY_TEST'):
        baseline={'cg':10000,'fcg':20000,'fgmres':40000}[family]
        transmitted=baseline//2 if compressed else baseline
    print('  transmitted sent bytes: '+str(transmitted))
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
reference = Path(sys.argv[2]).parent.name
error = {'uncompressed-fcg':3e-9, 'uncompressed-fgmres':4e-9}.get(reference,2e-9)
if reference.startswith('before-'): error=0
print('relative L2 error: '+str(error))
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
            observed = subprocess.run(
                command[:-1]+['observed','--bits','16','--cg-diagnostics'],
                cwd=repo,env=env,capture_output=True,text=True)
            self.assertEqual(observed.returncode,0,observed.stdout+observed.stderr)
            observed_calls=[json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()][-3:]
            self.assertNotIn('--bddcCgDiagnostics',observed_calls[0])
            for row in observed_calls[1:]:
                self.assertEqual(row[row.index('--bddcCgDiagnostics')+1],'1')
                self.assertEqual(row[row.index('--bddcRelativeTolerance')+1],'1e-07')
            self.assertTrue(json.loads((repo/'observed/environment.json').read_text())['cg_diagnostics'])
            zero_result = subprocess.run(
                command[:-1]+['zero-test','--bits','12','--full-compression',
                              '--zero-preserving-restriction','--refine','2','--ranks','8'],
                cwd=repo,env=env,capture_output=True,text=True)
            self.assertEqual(zero_result.returncode,0,zero_result.stdout+zero_result.stderr)
            zero_calls=[json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()][-4:]
            for row in zero_calls[:2]:
                self.assertNotIn('--bddcRestrictZeroPreserving',row)
            for row in zero_calls[2:]:
                self.assertEqual(row[row.index('--bddcRestrictZeroPreserving')+1],'1')
                self.assertEqual(row[row.index('--bddcRestrictionPasses')+1],'1')
                self.assertEqual(row[row.index('--bddcProlongateCompressionBits')+1],'12')
            zero_summary=json.loads((repo/'zero-test/summary.json').read_text())
            self.assertEqual([r['zero_preserving_restriction'] for r in zero_summary],
                             [False,False,True,True])
            self.assertTrue(json.loads((repo/'zero-test/environment.json').read_text())['zero_preserving_restriction'])
            range_result = subprocess.run(
                command[:-1]+['ranges','--bits','16','--restriction-range-capture'],
                cwd=repo,env=dict(env,KASKADE_RESTRICTION_CAPTURE='/should-not-leak'),
                capture_output=True,text=True)
            self.assertEqual(range_result.returncode,0,range_result.stdout+range_result.stderr)
            range_env=[json.loads(line) for line in (repo/'range-env.jsonl').read_text().splitlines()][-3:]
            self.assertEqual(range_env[:2],[None,None])
            self.assertEqual(Path(range_env[2]),(repo/'ranges/quantization-16/restriction-buffers').resolve())
            range_calls=[json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()][-3:]
            self.assertNotIn('KASKADE_RESTRICTION_CAPTURE',range_calls[0])
            self.assertNotIn('KASKADE_RESTRICTION_CAPTURE',range_calls[1])
            self.assertIn('KASKADE_RESTRICTION_CAPTURE',range_calls[2])
            grouped_result = subprocess.run(
                command[:-1]+['grouped','--bits','12','--full-compression',
                              '--zero-preserving-restriction','--restriction-group-size','128',
                              '--prolongation-bits','16'],
                cwd=repo,env=env,capture_output=True,text=True)
            self.assertEqual(grouped_result.returncode,0,grouped_result.stdout+grouped_result.stderr)
            grouped_calls=[json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()][-4:]
            for row in grouped_calls[:2]:
                self.assertNotIn('--bddcRestrictGroupSize',row)
            for row in grouped_calls[2:]:
                self.assertEqual(row[row.index('--bddcRestrictGroupSize')+1],'128')
                self.assertEqual(row[row.index('--bddcRestrictCompressionBits')+1],'12')
                self.assertEqual(row[row.index('--bddcProlongateCompressionBits')+1],'16')
                self.assertEqual(row[row.index('--bddcRestrictionPasses')+1],'1')
            grouped_summary=json.loads((repo/'grouped/summary.json').read_text())
            self.assertEqual(grouped_summary[-1]['prolongation_bits'],16)
            self.assertEqual(grouped_summary[-1]['restriction_group_size'],128)
            self.assertNotIn('Uniform bit depths apply', (repo/'grouped/summary.txt').read_text())
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
            # New mesh selection changes all four inputs, and each explicit
            # transform order has an independently checked flag combination.
            # Use a private parent, avoiding fixture files outside this test.
            pepe_repo=repo/'pepe-repo'
            pepe_repo.mkdir()
            shutil.copy2(model,pepe_repo/'emiModel')
            shutil.copy2(repo/'compare_emi_vtu.py',pepe_repo/'compare_emi_vtu.py')
            meshes=repo/'input_emi_mesh'
            meshes.mkdir()
            for suffix in ('.vtu','_extracellular.txt','_intracellular.txt','_excited.txt'):
                (meshes/('pepe_combi_domi'+suffix)).touch()
            paired=subprocess.run([sys.executable,str(runner),'--repo',str(pepe_repo),
                '--case','pepe','--ranks','2','--bits','12','16','24',
                '--cg-original-residual','1','--transform-orders','before','after',
                '--output','orders'],cwd=pepe_repo,env=env,capture_output=True,text=True)
            self.assertEqual(paired.returncode,0,paired.stdout+paired.stderr)
            order_rows=json.loads((pepe_repo/'orders/summary.json').read_text())
            self.assertEqual([r['case'] for r in order_rows],['direct','uncompressed-cg']+
                [f'{mode}-{b}' for b in (12,16,24) for mode in ('quantization','before','after')])
            cmds=[json.loads(line) for line in (pepe_repo/'capture.jsonl').read_text().splitlines()]
            for i,c in enumerate(cmds):
                def val(flag): return c[c.index(flag)+1]
                for flag,suffix in zip(('--input','--extra_set','--intra_set','--excited'),
                                      ('.vtu','_extracellular.txt','_intracellular.txt','_excited.txt')):
                    self.assertEqual(Path(val(flag)),(meshes/('pepe_combi_domi'+suffix)).resolve())
                if i:
                    self.assertEqual(val('--bddcCgOriginalResidual'),'1')
                    for flag in ('--bddcFlexibleCg','--bddcFgmres','--bddcProlongateConstantSplit'):
                        self.assertEqual(val(flag),'0')
                    self.assertEqual(val('-np'),'2')
                if i>=2:
                    mode=('quantization','before','after')[(i-2)%3]
                    self.assertEqual(val('--bddcGraphLifting'),str(int(mode=='before')))
                    self.assertEqual(val('--bddcGraphLiftingAfterQuantization'),str(int(mode=='after')))
                    for flag in ('--bddcHuffman','--bddcBitlength','--bddcZigzag'):
                        self.assertEqual(val(flag),str(int(mode!='quantization')))
            all_command=[sys.executable,str(runner),'--repo',str(pepe_repo),
                '--case','pepe','--ranks','2','--bits','12','16','24','--cg-original-residual','1',
                '--transform-orders','before','after','--solvers','cg','fcg','fgmres','--split','both',
                '--iterations','1000']
            matrix=subprocess.run(all_command+['--output','matrix'],cwd=pepe_repo,
                env=dict(env,FAMILY_TEST='1'),capture_output=True,text=True)
            self.assertEqual(matrix.returncode,1,matrix.stdout+matrix.stderr)
            matrix_rows=json.loads((pepe_repo/'matrix/summary.json').read_text())
            self.assertEqual(len(matrix_rows),58)
            self.assertEqual([r['case'] for r in matrix_rows if r['status']=='FAIL'],
                             ['before-cg12'])
            for row in matrix_rows[4:]:
                self.assertEqual(row['uncompressed_reference'],'uncompressed-'+row['solver'])
                if row['status']=='PASS':
                    self.assertAlmostEqual(row['net_logical_byte_savings_vs_uncompressed_percent'],50)
                    self.assertIn('comparison_to_direct',row)
                    self.assertIn('comparison_to_uncompressed',row)
                    self.assertEqual(row['comparison_to_direct']['relative L2 error'],2e-9)
                    self.assertEqual(row['comparison_to_uncompressed']['relative L2 error'],
                                     {'cg':2e-9,'fcg':3e-9,'fgmres':4e-9}[row['solver']])
            matrix_calls=[json.loads(line) for line in (pepe_repo/'capture.jsonl').read_text().splitlines()][-58:]
            for row,c in zip(matrix_rows,matrix_calls):
                if row['case']=='direct': continue
                self.assertEqual(c[c.index('--bddcIterations')+1],'1000')
                self.assertEqual(c[c.index('--bddcProlongateConstantSplit')+1],str(int(row['constant_split'])))
                self.assertEqual(c[c.index('--bddcFgmres')+1],str(int(row['solver']=='fgmres')))
                self.assertEqual(c[c.index('--bddcFlexibleCg')+1],str(int(row['solver']=='fcg')))
                if row['restriction_bits'] is not None:
                    self.assertEqual(c[c.index('--bddcGraphLifting')+1],str(int(row['transform_order']=='before')))
                    self.assertEqual(c[c.index('--bddcGraphLiftingAfterQuantization')+1],str(int(row['transform_order']=='after')))
            # Failure of one baseline must skip only its own family.
            combined=subprocess.run(all_command+['--output','combined','--bits','16',
                '--skip-quantization','--transform-orders','before','both'],cwd=pepe_repo,
                env=dict(env,FAMILY_TEST='1'),capture_output=True,text=True)
            self.assertEqual(combined.returncode,0,combined.stdout+combined.stderr)
            combined_rows=json.loads((pepe_repo/'combined/summary.json').read_text())
            self.assertEqual(len(combined_rows),16)
            pairs=[r['both_vs_before'] for r in combined_rows if r.get('transform_order')=='both']
            self.assertEqual(len(pairs),6)
            self.assertTrue(all(p['status']=='EXACT_EXPORTED_SOLUTION_AND_TRACE' for p in pairs))
            self.assertTrue(all(p['iteration_difference']==0 for p in pairs))
            combined_calls=[json.loads(line) for line in (pepe_repo/'capture.jsonl').read_text().splitlines()][-16:]
            for row,c in zip(combined_rows,combined_calls):
                self.assertEqual(c[c.index('--bddcGraphLiftingBothOrders')+1],
                                 str(int(row.get('transform_order')=='both')))
            # Failure of one baseline must skip only its own family.
            baseline_failed=subprocess.run(all_command+['--output','one-baseline-failed',
                '--bits','16','--skip-quantization'],cwd=pepe_repo,
                env=dict(env,FAMILY_TEST='1',FAIL_FCG_BASELINE='1'),capture_output=True,text=True)
            self.assertEqual(baseline_failed.returncode,1)
            failed_rows=json.loads((pepe_repo/'one-baseline-failed/summary.json').read_text())
            self.assertEqual(len(failed_rows),16)
            self.assertEqual(len([r for r in failed_rows if r['status']=='SKIP']),4)
            self.assertTrue(all(r['solver']=='fcg' for r in failed_rows if r['status']=='SKIP'))
            self.assertEqual(failed_rows[-1]['status'],'PASS')
            rejected=subprocess.run([sys.executable,str(runner),'--repo',str(pepe_repo),
                '--case','pepe','--transform-orders','before','--zero-preserving-restriction',
                '--output','invalid-before'],cwd=pepe_repo,env=env,capture_output=True,text=True)
            self.assertEqual(rejected.returncode,2)
            self.assertFalse((pepe_repo/'invalid-before').exists())
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
            # Exercise the audit driver with mock launches only. A failed
            # single-pass 12-bit case must not prevent the refinement check.
            testdir=repo/'tests'
            testdir.mkdir()
            for name in ['run_emi_quantization_diagnostic.py','collect_emi_scaling_smoke.py',
                         'run_emi_cg_diagnostic_check.sh']:
                shutil.copy2(runner.with_name(name),testdir/name)
            fake.write_text(fake.read_text().replace(
                "bad=compressed and value('--bddcCompressionBits')=='8' and not full",
                "bad=compressed and (value('--bddcCompressionBits')=='8' or "
                "(value('--bddcCompressionBits')=='12' and value('--bddcRestrictionPasses')=='1'))"))
            driver=['bash',str(testdir/'run_emi_cg_diagnostic_check.sh')]
            audit=subprocess.run(driver,cwd=repo,env=dict(env,OUTPUT='audit',RANKS='8',REFINE='2'),
                                 capture_output=True,text=True)
            self.assertEqual(audit.returncode,0,audit.stdout+audit.stderr)
            self.assertTrue((repo/'audit/observer-parity.txt').is_file())
            failed=json.loads((repo/'audit/b12-passes1/summary.json').read_text())
            passed=json.loads((repo/'audit/b12-passes2/summary.json').read_text())
            self.assertEqual(failed[-1]['status'],'FAIL')
            self.assertEqual(passed[-1]['status'],'PASS')
            # Deliberate instrumentation-dependent history change must stop
            # the driver before the low-bit experiments.
            fake.write_text(fake.read_text().replace(
                'rel=0.2 if bad else 7e-8',
                "rel=0.2 if bad else (8e-8 if '--bddcCgDiagnostics' in args else 7e-8)"))
            changed=subprocess.run(driver,cwd=repo,
                                   env=dict(env,OUTPUT='audit-changed',RANKS='8',REFINE='2'),
                                   capture_output=True,text=True)
            self.assertNotEqual(changed.returncode,0)
            self.assertIn('convergence history changed',changed.stdout+changed.stderr)
            self.assertFalse((repo/'audit-changed/b12-passes1').exists())


if __name__=='__main__':
    unittest.main()
