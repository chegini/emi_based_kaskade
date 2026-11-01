import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

from run_emi_flexible_cg_check import make_command,MARKER
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
            # A standard-CG regression must stop experimental runs, not pass silently.
            (ref/'10cells/full-r16-p16/bddc-convergence-step1.csv').write_text('different history\n')
            result=subprocess.run(args[:-1]+[str(repo/'changed')],capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            changed=json.loads((repo/'changed/summary.json').read_text())['10cells']
            self.assertEqual(list(changed),['uncompressed-cg','full-r16-p16'])
            self.assertEqual(changed['full-r16-p16']['status'],'FAIL')

if __name__=='__main__':unittest.main()
