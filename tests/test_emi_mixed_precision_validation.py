#!/usr/bin/env python3
"""Mock launches validate experiment controls; no PDE solves are performed."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from run_emi_mixed_precision_validation import objective, direct_errors, command, matched_command, history_parity, correction_metrics, logical_bytes


class MixedValidation(unittest.TestCase):
    def test_constant_correction_flags_and_payload(self):
        for mode in ['direct','uncompressed-cg','full-r16-p16','full-zero-r14-p16']:
            args=(Path('/model'),['mpirun'],Path('/mesh'),0,10,8,mode,Path('/out'))
            old=command(*args,candidate_bits=14,matched_original_rho=True)
            new=command(*args,candidate_bits=14,matched_original_rho=True,matched_constant_correction=True)
            i=new.index('--bddcConstantCorrection')
            self.assertEqual(new[i+1],'1' if mode.startswith('full-') else '0')
            self.assertEqual(old,new[:i]+new[i+2:])
        text=('BDDC constant correction step 1: delta=1e-8 sum_r=1e-14 one_A_one=1e-6 pre_rel=8e-8 post_rel=7e-8\n'
              'BDDC constant correction additional broadcast payload bytes: 224\n')
        rows,extra,errors=correction_metrics(text,1,8)
        self.assertFalse(errors); self.assertEqual(extra,224); self.assertEqual(len(rows),1)
        self.assertTrue(correction_metrics(text.replace('post_rel=7e-8','post_rel=2e-7'),1,8)[2])
        self.assertTrue(correction_metrics(text.replace('224','0'),1,8)[2])
        self.assertTrue(correction_metrics(text,2,8)[2])
        self.assertEqual(logical_bytes({'metrics':{'transmitted_sent_bytes':100,'constant_correction_payload_bytes':224}}),324)

    def test_fourteen_bit_changes_only_candidate_restriction(self):
        args=(Path('/repo/emiModel'),['mpirun'],Path('/mesh'),0,10,8)
        old=command(*args,'full-zero-r12-p16',Path('/out'),
                    matched_original_rho=True,matched_cg_diagnostics=True)
        new=command(*args,'full-zero-r14-p16',Path('/out'),
                    matched_original_rho=True,matched_cg_diagnostics=True,candidate_bits=14)
        changed=[i for i,(a,b) in enumerate(zip(old,new)) if a!=b]
        self.assertEqual(changed,[old.index('--bddcRestrictCompressionBits')+1])
        self.assertEqual(new[changed[0]],'14')
        self.assertFalse(matched_command(old,new))
        self.assertTrue(matched_command(old,new,allow_restriction_bits=True))
        new[new.index('--bddcProlongateCompressionBits')+1]='12'
        self.assertFalse(matched_command(old,new,allow_restriction_bits=True))
        for mode in ['direct','uncompressed-cg','full-r16-p16']:
            self.assertEqual(command(*args,mode,Path('/out')),
                             command(*args,mode,Path('/out'),candidate_bits=14))

    def test_matched_scalar_and_observer_flags(self):
        for mode in ['direct','uncompressed-cg','full-r16-p16','full-zero-r12-p16']:
            args=(Path('/repo/emiModel'),['mpirun'],Path('/mesh'),0,10,8,mode,Path('/out'))
            old=command(*args)
            new=command(*args,matched_original_rho=True,matched_cg_diagnostics=True)
            compressed=mode.startswith('full-')
            self.assertEqual(new[new.index('--bddcCgOriginalResidual')+1],str(int(compressed)))
            self.assertEqual(new[new.index('--bddcCgDiagnostics')+1],str(int(compressed)))
            self.assertTrue(matched_command(old,new))
            for flag in ['--bddcRestrictCompressionBits','--bddcProlongateCompressionBits',
                         '--bddcRestrictZeroPreserving','--bddcRestrictionPasses',
                         '--bddcRelativeTolerance','--sigma_i','--sigma_e','--dt']:
                self.assertEqual(new[new.index(flag)+1],old[old.index(flag)+1])

    def test_opt_in_scalar_only_changes_candidate(self):
        for mode in ['direct','uncompressed-cg','full-r16-p16','full-zero-r12-p16']:
            args=(Path('/repo/emiModel'),['mpirun'],Path('/mesh'),0,10,8,mode,Path('/out'))
            old=command(*args)
            new=command(*args,candidate_original_rho=True,candidate_cg_diagnostics=True)
            self.assertNotIn('--bddcCgOriginalResidual',old)
            candidate=mode=='full-zero-r12-p16'
            self.assertEqual(new[new.index('--bddcCgOriginalResidual')+1],str(int(candidate)))
            self.assertEqual(new[new.index('--bddcCgDiagnostics')+1],str(int(candidate)))
            self.assertTrue(matched_command(old,new))
            changed=list(new);changed[changed.index('--bddcRestrictionPasses')+1]='2'
            self.assertFalse(matched_command(old,changed))
        with tempfile.TemporaryDirectory() as folder:
            a=Path(folder)/'a'; b=Path(folder)/'b';a.mkdir();b.mkdir()
            self.assertFalse(history_parity(a,b))
            for d in [a,b]: (d/'bddc-convergence-step1.csv').write_text('same\n')
            self.assertTrue(history_parity(a,b))
            (b/'bddc-convergence-step1.csv').write_text('changed\n')
            self.assertFalse(history_parity(a,b))

    def test_accuracy_and_byte_gates(self):
        baseline=dict(status='PASS',metrics=dict(transmitted_sent_bytes=100),
                      comparison_to_direct={'relative L2 error':2e-6,'max absolute error':3e-7})
        candidate=dict(status='PASS',metrics=dict(transmitted_sent_bytes=90),
                       comparison_to_direct={'relative L2 error':1e-6,'max absolute error':1e-7})
        self.assertEqual(objective(baseline,candidate),[])
        candidate['comparison_to_direct']['max absolute error']=4e-7
        self.assertTrue(objective(baseline,candidate))
        candidate['comparison_to_direct']['max absolute error']=1e-7
        candidate['metrics']['transmitted_sent_bytes']=101
        self.assertTrue(objective(baseline,candidate))
        self.assertTrue(direct_errors('Direct global residual step 2: abs=1e-15 rel=1e-15\nEMI model completed\n',2,0))

    def test_multistep_campaign_and_early_failure(self):
        runner=Path(__file__).with_name('run_emi_mixed_precision_validation.py')
        with tempfile.TemporaryDirectory(prefix='emi-mixed-test-') as temp:
            repo=Path(temp).resolve()
            for base in [repo/'input/10Cells3d_10extra',repo/'meshes/weak_8_h16_kaskade',repo/'meshes/strong_h16_kaskade']:
                base.parent.mkdir(exist_ok=True)
                for suffix in ['_mesh.vtu','_list_extracellular.txt','_list_intracellular.txt','_early_excited.txt']:
                    Path(str(base)+suffix).touch()
            model=repo/'emiModel';model.write_text('#!/bin/sh\nexit 0\n');model.chmod(0o755)
            fake=repo/'fake-mpi'
            fake.write_text("""#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
a=sys.argv[1:]
def v(k): return a[a.index(k)+1]
with open('commands.jsonl','a') as f: f.write(json.dumps(a)+'\\n')
assert 'KASKADE_RESTRICTION_CAPTURE' not in os.environ
steps=int(v('--maximumNumberOfTimeSteps')); out=Path(v('--dir'))
direct=v('--direct')=='1'; compressed=v('--bddcCompression')=='1'
candidate=v('--bddcRestrictZeroPreserving')=='1'
if '--potentialDiagnostics' in a:
    import csv
    with (out/'potential-operator.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['dof','a_one','a_transpose_one','diagonal'])
        w.writerows([i,1,1,1] for i in range(30))
    for step in range(steps+1):
        with (out/f'potential-state-step{step}.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['dof','u','correction','rhs','residual'])
            w.writerows([i,step*.1,.1 if step else 0,.1 if step else 0,0] for i in range(30))
if direct:
    assert v('--directBackend')=='mumps-spd'
    assert v('--bddc')=='0' and v('-np')=='1'
    assert v('--bddcRelativeTolerance')=='0'
    assert v('--bddcMaterialAssembly')=='0' and v('--bddcInterfaceAssembly')=='0'
    for k in range(1,steps+1): print(f'Direct global residual step {k}: abs=1e-15 rel=1e-15')
    (out/'emiLast.vtu').touch()
else:
    assert '--directBackend' not in a, 'Direct-only backend leaked into a BDDC command'
    assert v('--bddcRelativeTolerance')=='1e-7'
    print('cells: 100\\ndofs: 30\\nBDDC subdomains: 4, shared dofs: 10')
    print('BDDC coarse interfaces: corners=0, edges=2, faces=3')
    for k in range(1,steps+1):
        rel=2e-6 if candidate and k==1 and os.environ.get('FAIL_EARLY') else 7e-8
        print(f'BDDC global residual step {k}: abs=1e-10 rel={rel} target=1e-7 status=converged iterations={11 if candidate else 10}')
        (out/f'bddc-convergence-step{k}.csv').write_text(f'{rel}\\n')
        if v('--bddcCgDiagnostics')=='1':
            (out/f'bddc-cg-diagnostics-step{k}.csv').write_text('mock observer\\n')
        if '--bddcConstantCorrection' in a and v('--bddcConstantCorrection')=='1':
            print(f'BDDC constant correction step {k}: delta=1e-8 sum_r=1e-14 one_A_one=1e-6 pre_rel=8e-8 post_rel=7e-8')
    if '--bddcConstantCorrection' in a and v('--bddcConstantCorrection')=='1':
        print('BDDC constant correction additional broadcast payload bytes: '+str(32*(int(v('-np'))-1)*steps))
    print('raw sent bytes: 10000')
    print('transmitted sent bytes: '+str(800 if candidate else 1000 if compressed else 10000))
    print(f'Profile (wall seconds, EMI end-to-end phases; {v("-np")} rank(s)):')
    print('  solver wrapper (including setup) 1.0 1.0 1.0')
    (out/('emiBDDCLastCompression.vtu' if compressed else 'emiBDDCLast.vtu')).touch()
print('total cpu-time: 2s wall, 1s CPU')
print('EMI model completed')
""")
            fake.chmod(0o755)
            (repo/'compare_emi_vtu.py').write_text("""import sys
from pathlib import Path
assert all(Path(x).is_file() for x in sys.argv[1:])
print('same coordinates: yes')
print('relative L2 error: '+('1e-6' if 'zero' in sys.argv[1] else '2e-6'))
print('max absolute error: '+('1e-7' if 'zero' in sys.argv[1] else '3e-7'))
print('mean signed offset (first - reference): 0')
print('relative L2 after removing mean offset: 1e-6')
""")
            env=dict(os.environ,MPIEXEC=str(fake),KASKADE_RESTRICTION_CAPTURE='/should-not-leak')
            cmd=[sys.executable,str(runner),'--repo',str(repo),'--mesh-dir','meshes',
                 '--steps','2','--ranks','8']
            run=subprocess.run(cmd+['--output','success'],env=env,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            results=json.loads((repo/'success/summary.json').read_text())
            self.assertTrue(all(r['objective_passed'] for r in results.values()))
            calls=[json.loads(s) for s in (repo/'commands.jsonl').read_text().splitlines()]
            self.assertEqual(len(calls),12) # Direct/uncompressed references run only once per mesh.
            for a in calls:
                v=lambda k:a[a.index(k)+1]
                self.assertEqual(v('--finalTime'),'0.02')
                self.assertEqual(v('--maximumNumberOfTimeSteps'),'2')
                self.assertEqual(v('--bddcRestrictionPasses'),'1')
                self.assertEqual(v('--bddcRestrictGroupSize'),'0')
                self.assertEqual(v('--bddcProlongateCompressionBits'),'16')
                self.assertEqual(v('--refine'),'2' if '10Cells' in v('--input') else '0')
                self.assertEqual(v('-np'),'1' if v('--direct')=='1' else '8')
                if v('--bddcRestrictZeroPreserving')=='1':
                    self.assertEqual(v('--bddcRestrictCompressionBits'),'12')
                    self.assertEqual(v('--bddcHuffman'),'1')
            experiment=subprocess.run(cmd+['--output','rho','--cases','weak_8_h16',
                '--candidate-original-rho','--candidate-cg-diagnostics',
                '--reference-campaign',str(repo/'success')],
                env=env,capture_output=True,text=True)
            self.assertEqual(experiment.returncode,0,experiment.stdout+experiment.stderr)
            rho=json.loads((repo/'rho/summary.json').read_text())['weak_8_h16']['runs']
            self.assertTrue(rho['full-r16-p16']['legacy_history_parity'])
            self.assertTrue(rho['full-zero-r12-p16']['original_residual_rho'])
            self.assertEqual(len(rho['full-zero-r12-p16']['cg_diagnostic_files']),2)
            matched=subprocess.run(cmd+['--output','matched','--cases','weak_8_h16',
                '--matched-original-rho','--matched-cg-diagnostics',
                '--reference-campaign',str(repo/'rho')],
                env=env,capture_output=True,text=True)
            self.assertEqual(matched.returncode,0,matched.stdout+matched.stderr)
            runs=json.loads((repo/'matched/summary.json').read_text())['weak_8_h16']['runs']
            for mode in ['full-r16-p16','full-zero-r12-p16']:
                self.assertTrue(runs[mode]['original_residual_rho'])
                self.assertEqual(len(runs[mode]['cg_diagnostic_files']),2)
            self.assertIn('history_parity_note',runs['full-r16-p16'])
            self.assertNotIn('saved_history_parity',runs['full-r16-p16'])
            self.assertTrue(runs['full-zero-r12-p16']['saved_history_parity'])
            self.assertTrue(runs['uncompressed-cg']['saved_history_parity'])
            fourteen=subprocess.run(cmd+['--output','fourteen','--cases','weak_8_h16',
                '--candidate-restriction-bits','14',
                '--matched-original-rho','--matched-cg-diagnostics',
                '--reference-campaign',str(repo/'matched')],
                env=env,capture_output=True,text=True)
            self.assertEqual(fourteen.returncode,0,fourteen.stdout+fourteen.stderr)
            runs14=json.loads((repo/'fourteen/summary.json').read_text())['weak_8_h16']['runs']
            self.assertNotIn('full-zero-r12-p16',runs14)
            self.assertTrue(runs14['full-r16-p16']['saved_history_parity'])
            self.assertIn('history_parity_note',runs14['full-zero-r14-p16'])
            self.assertNotIn('saved_history_parity',runs14['full-zero-r14-p16'])
            self.assertIn('R14/P16', (repo/'fourteen/summary.txt').read_text())
            self.assertEqual(json.loads((repo/'fourteen/environment.json').read_text())['candidate_restriction_bits'],14)
            import shlex
            c14=shlex.split((repo/'fourteen/weak_8_h16/full-zero-r14-p16/command.sh').read_text())
            self.assertEqual(c14[c14.index('--bddcRestrictCompressionBits')+1],'14')
            observed=subprocess.run(cmd+['--output','potential','--cases','weak_8_h16',
                '--candidate-restriction-bits','14','--potential-diagnostics',
                '--matched-original-rho','--matched-cg-diagnostics',
                '--reference-campaign',str(repo/'fourteen')],
                env=env,capture_output=True,text=True)
            self.assertEqual(observed.returncode,0,observed.stdout+observed.stderr)
            self.assertTrue(json.loads((repo/'potential/potential-analysis/interpretation.json').read_text())['complete'])
            observed_runs=json.loads((repo/'potential/summary.json').read_text())['weak_8_h16']['runs']
            self.assertTrue(observed_runs['full-zero-r14-p16']['saved_history_parity'])
            corrected=subprocess.run(cmd+['--output','corrected','--cases','weak_8_h16',
                '--candidate-restriction-bits','14','--potential-diagnostics',
                '--matched-original-rho','--matched-cg-diagnostics','--matched-constant-correction',
                '--reference-campaign',str(repo/'potential')],
                env=env,capture_output=True,text=True)
            self.assertEqual(corrected.returncode,0,corrected.stdout+corrected.stderr)
            corrected_runs=json.loads((repo/'corrected/summary.json').read_text())['weak_8_h16']['runs']
            for mode in ['full-r16-p16','full-zero-r14-p16']:
                r=corrected_runs[mode]
                self.assertTrue(r['first_step_cg_parity'])
                self.assertNotIn('saved_history_parity',r)
                self.assertTrue(r['constant_correction'])
                self.assertEqual(r['metrics']['constant_correction_payload_bytes'],448)
            self.assertTrue(corrected_runs['uncompressed-cg']['saved_history_parity'])
            # Reject changes in a saved same-scalar history; do not silently
            # waive reproducibility just because the experiment is matched.
            (repo/'rho/weak_8_h16/uncompressed-cg/bddc-convergence-step1.csv').write_text('different\n')
            mismatch=subprocess.run(cmd+['--output','mismatch','--cases','weak_8_h16',
                '--matched-original-rho','--matched-cg-diagnostics',
                '--reference-campaign',str(repo/'rho')],
                env=env,capture_output=True,text=True)
            self.assertEqual(mismatch.returncode,1,mismatch.stdout+mismatch.stderr)
            stopped=json.loads((repo/'mismatch/summary.json').read_text())['weak_8_h16']['runs']
            self.assertNotIn('full-r16-p16',stopped)
            failed=subprocess.run(cmd+['--output','early-failure','--cases','10cells'],
                                  env=dict(env,FAIL_EARLY='1'),capture_output=True,text=True)
            self.assertEqual(failed.returncode,1,failed.stdout+failed.stderr)
            result=json.loads((repo/'early-failure/summary.json').read_text())['10cells']
            self.assertFalse(result['objective_passed'])
            errors=result['runs']['full-zero-r12-p16']['errors']
            self.assertTrue(any('Step 1 did not meet' in s for s in errors))
            before=(repo/'commands.jsonl').read_text()
            reused=subprocess.run(cmd+['--output','success'],env=env,capture_output=True,text=True)
            self.assertNotEqual(reused.returncode,0)
            self.assertEqual(before,(repo/'commands.jsonl').read_text())


if __name__=='__main__':
    unittest.main()
