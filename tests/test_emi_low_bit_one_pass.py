import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from run_emi_low_bit_one_pass import candidate_command
from run_emi_mixed_precision_validation import command,digest


class LowBits(unittest.TestCase):
    def test_change_only_bits_and_output(self):
        seed=command(Path('/model'),['mpi'],Path('/mesh'),2,10,8,'full-zero-r14-p16',
                     Path('/old'),matched_original_rho=True,matched_cg_diagnostics=True,
                     candidate_bits=14,potential_diagnostics=True,matched_constant_correction=True)
        new=candidate_command(seed,8,Path('/new'))
        self.assertEqual([i for i,(a,b) in enumerate(zip(seed,new)) if a!=b],
                         [seed.index('--bddcRestrictCompressionBits')+1,seed.index('--dir')+1])
        r12=candidate_command(seed,12,Path('/old'))
        r12p12=candidate_command(r12,12,Path('/new'),12)
        self.assertEqual([i for i,(a,b) in enumerate(zip(r12,r12p12)) if a!=b],
                         [r12.index('--bddcProlongateCompressionBits')+1,r12.index('--dir')+1])
        seed[seed.index('--bddcRestrictionPasses')+1]='2'
        with self.assertRaises(ValueError):candidate_command(seed,8,Path('/new'))

    def test_reuses_controls_and_continues_after_failure(self):
        script=Path(__file__).with_name('run_emi_low_bit_one_pass.py')
        with tempfile.TemporaryDirectory() as temp:
            repo=Path(temp).resolve(); ref=repo/'saved';ref.mkdir()
            model=repo/'emiModel';model.write_text('mock executable')
            fake=repo/'fake-mpi'
            fake.write_text("""#!/usr/bin/env python3
import sys
from pathlib import Path
a=sys.argv[1:]
def v(k):return a[a.index(k)+1]
bits=int(v('--bddcRestrictCompressionBits'));out=Path(v('--dir'))
assert v('--bddcRestrictionPasses')=='1' and v('--bddcProlongateCompressionBits') in ['16','12']
with open('calls','a') as f:f.write(str(bits)+'\\n')
if bits==8: print('Synthetic iteration-limit failure');sys.exit(134)
print('cells: 100\\ndofs: 2\\nBDDC subdomains: 2, shared dofs: 1')
print('BDDC coarse interfaces: corners=0, edges=1, faces=0')
for k in [1,2]:
 print(f'BDDC global residual step {k}: abs=1e-10 rel=7e-8 target=1e-7 status=converged iterations=11')
 print(f'BDDC constant correction step {k}: delta=1e-8 sum_r=1e-14 one_A_one=1e-6 pre_rel=8e-8 post_rel=7e-8')
print('BDDC constant correction additional broadcast payload bytes: 448')
print('raw sent bytes: 10000\\ntransmitted sent bytes: 600')
print('Profile (wall seconds, EMI end-to-end phases; 8 rank(s)):')
print('  solver wrapper (including setup) 1.0 1.0 1.0')
print('total cpu-time: 2s wall, 1s CPU\\nEMI model completed')
(out/'emiBDDCLastCompression.vtu').touch()
(out/'potential-state-step2.csv').write_text('dof,u\\n0,1.00000001\\n1,2\\n')
""")
            fake.chmod(0o755)
            (repo/'compare_emi_vtu.py').write_text("""print('same coordinates: yes')
print('relative L2 error: 1e-5')
print('max absolute error: 1e-5')
print('mean signed offset (first - reference): 1e-5')
print('relative L2 after removing mean offset: 1e-9')
""")
            modes=['direct','uncompressed-cg','full-r16-p16','full-zero-r14-p16']
            runs={}
            for mode in modes:
                out=ref/'10cells'/mode;out.mkdir(parents=True)
                cmd=command(model,[str(fake)],Path('/mesh'),2,2,8,mode,out,
                    matched_original_rho=True,matched_cg_diagnostics=True,candidate_bits=14,
                    potential_diagnostics=True,matched_constant_correction=True)
                (out/'command.sh').write_text(shlex.join(cmd))
                runs[mode]={'status':'PASS','metrics':{'transmitted_sent_bytes':1000,
                    'constant_correction_payload_bytes':448},
                    'comparison_to_direct':{'relative L2 error':1e-7,'max absolute error':1e-7}}
            (ref/'10cells/direct/emiLast.vtu').touch()
            (ref/'10cells/uncompressed-cg/emiBDDCLast.vtu').touch()
            (ref/'10cells/direct/potential-state-step2.csv').write_text('dof,u\n0,1\n1,2\n')
            (ref/'summary.json').write_text(json.dumps({'10cells':{'runs':runs}}))
            env={'executable_sha256':digest(model),'steps':2,'bddc_ranks':8,'candidate_restriction_bits':14}
            (ref/'environment.json').write_text(json.dumps(env))
            old=(ref/'summary.json').read_bytes()
            args=[sys.executable,str(script),'--repo',str(repo),'--reference-campaign',str(ref),
                  '--cases','10cells','--bits','12','8','10','--output',str(repo/'new')]
            result=subprocess.run(args,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            rows=json.loads((repo/'new/summary.json').read_text())
            self.assertEqual([r['status'] for r in rows],['PASS','FAIL','PASS'])
            self.assertTrue(rows[0]['communication_better_than16'])
            self.assertTrue(rows[0]['legacy_strict_objective_failures']) # explicitly not a veto now
            self.assertEqual((repo/'calls').read_text(),'12\n8\n10\n')
            self.assertEqual((ref/'summary.json').read_bytes(),old)
            paired_args=[sys.executable,str(script),'--repo',str(repo),'--reference-campaign',str(ref),
                         '--comparison-campaign',str(repo/'new'),'--cases','10cells','--bits','12',
                         '--prolongation-bits','12','--output',str(repo/'paired')]
            result=subprocess.run(paired_args,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            paired=json.loads((repo/'paired/summary.json').read_text())[0]
            self.assertEqual(paired['prolongation_bits'],12)
            self.assertEqual(paired['status'],'PASS')
            self.assertEqual(paired['saved_prolongation_control']['saved_percent'],0)
            self.assertEqual(paired['saved_prolongation_control']['prolongation_bits'],16)
            self.assertEqual((repo/'calls').read_text(),'12\n8\n10\n12\n')
            saved_command=repo/'new/10cells/full-zero-r12-p16/command.sh'
            saved_command.write_text(saved_command.read_text().replace('--bddcRestrictionPasses 1',
                                                                       '--bddcRestrictionPasses 2'))
            result=subprocess.run(paired_args[:-1]+[str(repo/'mismatched')],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((repo/'mismatched').exists())
            # A stale binary must never reuse these controls.
            model.write_text('different executable')
            result=subprocess.run(args[:-1]+[str(repo/'stale')],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((repo/'stale').exists())


if __name__=='__main__':unittest.main()
