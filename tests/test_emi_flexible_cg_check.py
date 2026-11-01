import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

from run_emi_flexible_cg_check import make_command,MARKER,ZERO_P_MARKER
from run_emi_mixed_precision_validation import command


class FlexibleCheck(unittest.TestCase):
    def test_only_authorized_changes(self):
        seed=command(Path('/model'),['mpi'],Path('/mesh'),2,10,8,'full-zero-r14-p16',Path('/old'),
                     matched_original_rho=True,matched_cg_diagnostics=True,candidate_bits=14,
                     potential_diagnostics=True,matched_constant_correction=True)
        new=make_command(seed,'flex-zero-r12-p12',Path('/new'))
        self.assertEqual(new[-2:],['--bddcFlexibleCg','1'])
        differences={seed[i-1] for i,(a,b) in enumerate(zip(seed,new)) if a!=b}
        self.assertEqual(differences,{'--bddcRestrictCompressionBits','--bddcProlongateCompressionBits','--dir'})
        self.assertEqual(new[new.index('--bddcRestrictionPasses')+1],'1')
        self.assertEqual(new[new.index('--bddcRelativeTolerance')+1],'1e-7')
        eight=make_command(seed,'flex-zero-r8-p8',Path('/eight'))
        self.assertEqual(eight[eight.index('--bddcRestrictCompressionBits')+1],'8')
        self.assertEqual(eight[eight.index('--bddcProlongateCompressionBits')+1],'8')
        self.assertEqual(eight[eight.index('--bddcRestrictionPasses')+1],'1')
        zero_p=make_command(seed,'flex-zero-r8-p8-zeroP',Path('/eight'))
        i=zero_p.index('--bddcProlongateZeroPreserving')
        self.assertEqual(zero_p[i:i+2],['--bddcProlongateZeroPreserving','1'])
        self.assertEqual(zero_p[:i]+zero_p[i+2:],eight)
        for mode,changed_flag in [('flex-zero-r8-p16','--bddcProlongateCompressionBits'),
                                  ('flex-zero-r16-p8','--bddcRestrictCompressionBits')]:
            mixed=make_command(seed,mode,Path('/eight'))
            differences={eight[i-1] for i,(a,b) in enumerate(zip(eight,mixed)) if a!=b}
            self.assertEqual(differences,{changed_flag})
            self.assertEqual(mixed[mixed.index('--bddcRestrictZeroPreserving')+1],'1')
        with self.assertRaises(ValueError):make_command(seed,'flex-zero-r1-p8',Path('/bad'))
        with self.assertRaises(ValueError):make_command(new,'flex-zero-r12-p12',Path('/again'))

    def test_controls_parity_and_failure_continuation(self):
        script=Path(__file__).with_name('run_emi_flexible_cg_check.py')
        with tempfile.TemporaryDirectory() as tmp:
            repo=Path(tmp).resolve();ref=repo/'saved';ref.mkdir()
            model=repo/'emiModel';model.write_text('new executable')
            fake=repo/'mpi'
            fake.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
a=sys.argv[1:]
def v(k):return a[a.index(k)+1]
out=Path(v('--dir'));flex=v('--bddcFlexibleCg')=='1';bits=v('--bddcRestrictCompressionBits')
with open('calls','a') as f:f.write(out.name+'\\n')
if flex:print('''+repr(MARKER)+''')
if '--bddcProlongateZeroPreserving' in a:print('''+repr(ZERO_P_MARKER)+''')
if flex and bits=='16':print('Synthetic flexible-control failure');sys.exit(134)
print('cells: 100\\ndofs: 2\\nBDDC subdomains: 2, shared dofs: 1')
print('BDDC coarse interfaces: corners=0, edges=1, faces=0')
for k in [1,2]:
 print(f'BDDC global residual step {k}: abs=1e-10 rel=7e-8 target=1e-7 status=converged iterations=11')
 (out/f'bddc-convergence-step{k}.csv').write_text('same history\\n')
 if v('--bddcConstantCorrection')=='1':
  print(f'BDDC constant correction step {k}: delta=1e-8 sum_r=1e-14 one_A_one=1e-6 pre_rel=8e-8 post_rel=7e-8')
if v('--bddcConstantCorrection')=='1':print('BDDC constant correction additional broadcast payload bytes: 448')
print('raw sent bytes: 10000\\ntransmitted sent bytes: 600')
print('Profile (wall seconds, EMI end-to-end phases; 8 rank(s)):')
print('  solver wrapper (including setup) 1.0 1.0 1.0')
print('total cpu-time: 2s wall, 1s CPU\\nEMI model completed')
for name in ['emiBDDCLast.vtu','emiBDDCLastCompression.vtu']:(out/name).touch()
(out/'potential-state-step2.csv').write_text('dof,u\\n0,1.00000001\\n1,2\\n')
''')
            fake.chmod(0o755)
            (repo/'compare_emi_vtu.py').write_text("print('same coordinates: yes')\nprint('relative L2 error: 1e-8')\nprint('max absolute error: 1e-8')\nprint('mean signed offset (first - reference): 0')\nprint('relative L2 after removing mean offset: 1e-8')\n")
            runs={}
            for mode in ['direct','uncompressed-cg','full-r16-p16','full-zero-r14-p16']:
                out=ref/'10cells'/mode;out.mkdir(parents=True)
                cmd=command(model,[str(fake)],Path('/mesh'),2,2,8,mode,out,
                    matched_original_rho=True,matched_cg_diagnostics=True,candidate_bits=14,
                    potential_diagnostics=True,matched_constant_correction=True)
                (out/'command.sh').write_text(shlex.join(cmd))
                for step in [1,2]:(out/f'bddc-convergence-step{step}.csv').write_text('same history\n')
                runs[mode]=dict(status='PASS',metrics=dict(transmitted_sent_bytes=600,
                    constant_correction_payload_bytes=448 if mode.startswith('full-') else 0))
            (ref/'10cells/direct/emiLast.vtu').touch()
            (ref/'10cells/direct/potential-state-step2.csv').write_text('dof,u\n0,1\n1,2\n')
            (ref/'summary.json').write_text(json.dumps({'10cells':{'runs':runs}}))
            (ref/'environment.json').write_text(json.dumps(dict(steps=2,bddc_ranks=8,candidate_restriction_bits=14)))
            old=(ref/'summary.json').read_bytes()
            args=[sys.executable,str(script),'--repo',str(repo),'--reference-campaign',str(ref),
                  '--cases','10cells','--output',str(repo/'new')]
            result=subprocess.run(args,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            rows=json.loads((repo/'new/summary.json').read_text())['10cells']
            self.assertEqual([r['status'] for r in rows.values()],['PASS','PASS','FAIL','PASS'],rows)
            self.assertEqual(rows['flex-zero-r12-p12']['saved_vs_standard16_percent'],0)
            self.assertIsNone(rows['flex-zero-r12-p12']['saved_vs_flexible16_percent'])
            self.assertEqual((ref/'summary.json').read_bytes(),old)
            # A passing archive can supply controls to an 8/8 run without repeating them.
            fake.write_text(fake.read_text().replace("if flex and bits=='16':",'if False:'))
            result=subprocess.run(args[:-1]+[str(repo/'passing')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            before=(repo/'calls').read_text()
            reuse_args=args[:-1]+[str(repo/'eight'),'--bits','8',
                                 '--reuse-controls-from',str(repo/'passing')]
            result=subprocess.run(reuse_args,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertEqual((repo/'calls').read_text(),before+'flex-zero-r8-p8\n')
            eight=json.loads((repo/'eight/summary.json').read_text())['10cells']
            self.assertIn('reused_from',eight['flex-r16-p16'])
            self.assertEqual(eight['flex-zero-r8-p8']['status'],'PASS')
            # Changed executable must reject reuse before creating output or launching MPI.
            model.write_text('changed binary')
            rejected=repo/'rejected'
            reject_args=[str(rejected) if x==str(repo/'eight') else x for x in reuse_args]
            result=subprocess.run(reject_args,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse(rejected.exists())
            self.assertEqual((repo/'calls').read_text(),before+'flex-zero-r8-p8\n')
            model.write_text('new executable')
            # Directional mode launches BOTH candidates and retains the first failure.
            fake.write_text(fake.read_text().replace('if False:',"if flex and bits=='8':"))
            directional_args=args[:-1]+[str(repo/'directional'),'--directional-check',
                                       '--reuse-controls-from',str(repo/'passing')]
            before=(repo/'calls').read_text()
            result=subprocess.run(directional_args,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            directional=json.loads((repo/'directional/summary.json').read_text())['10cells']
            self.assertEqual(directional['flex-zero-r8-p16']['status'],'FAIL')
            self.assertEqual(directional['flex-zero-r16-p8']['status'],'PASS')
            self.assertEqual((repo/'calls').read_text(),before+'flex-zero-r8-p16\nflex-zero-r16-p8\n')
            denv=json.loads((repo/'directional/environment.json').read_text())
            self.assertTrue(denv['directional_check'])
            self.assertIsNone(denv['candidate_bits'])
            self.assertEqual(denv['candidate_modes'],['flex-zero-r8-p16','flex-zero-r16-p8'])
            # A failed original 8/8 run must not prevent the new quantizer experiment.
            fake.write_text(fake.read_text().replace("if flex and bits=='8':",
                "if flex and bits=='8' and '--bddcProlongateZeroPreserving' not in a:"))
            before=(repo/'calls').read_text()
            qargs=args[:-1]+[str(repo/'zeroP'),'--prolongation-quantizer-check',
                            '--reuse-controls-from',str(repo/'passing')]
            result=subprocess.run(qargs,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            qrows=json.loads((repo/'zeroP/summary.json').read_text())['10cells']
            self.assertEqual(qrows['flex-zero-r8-p8']['status'],'FAIL')
            self.assertEqual(qrows['flex-zero-r8-p8-zeroP']['status'],'PASS')
            self.assertTrue(qrows['flex-zero-r8-p8-zeroP']['zero_preserving_prolongation'])
            self.assertEqual((repo/'calls').read_text(),before+'flex-zero-r8-p8\nflex-zero-r8-p8-zeroP\n')
            # A standard-CG regression must stop experimental runs, not pass silently.
            (ref/'10cells/full-r16-p16/bddc-convergence-step1.csv').write_text('different history\n')
            result=subprocess.run(args[:-1]+[str(repo/'changed')],capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            changed=json.loads((repo/'changed/summary.json').read_text())['10cells']
            self.assertEqual(list(changed),['uncompressed-cg','full-r16-p16'])
            self.assertEqual(changed['full-r16-p16']['status'],'FAIL')

if __name__=='__main__':unittest.main()
