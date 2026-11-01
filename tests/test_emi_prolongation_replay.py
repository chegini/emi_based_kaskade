import json
import csv
import contextlib
import io
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

from replay_emi_prolongation_diagnostics import replay_command, signature, inspect_grouping
from run_emi_mixed_precision_validation import command
from run_emi_sd_cg_comparison import make_command


class ReplayTest(unittest.TestCase):
    def test_only_observer_changes(self):
        seed=command(Path('/model'),['mpi'],Path('/mesh'),2,10,8,'full-r16-p16',Path('/old'),
                     potential_diagnostics=True,matched_original_rho=True)
        off=replay_command(seed,Path('/same'),False)
        on=replay_command(seed,Path('/same'),True)
        self.assertEqual({on[i-1] for i,(a,b) in enumerate(zip(off,on)) if a!=b},
                         {'--bddcProlongationDiagnostics'})
        self.assertEqual(seed[seed.index('--dir')+1],'/old')
        probe=replay_command(seed,Path('/same'),True,7)
        self.assertEqual(probe[probe.index('--bddcProlongationConstraintCall')+1],'7')
        probeoff=replay_command(seed,Path('/same'),False,7)
        self.assertEqual(probeoff[probeoff.index('--bddcProlongationConstraintCall')+1],'0')
        grouped=replay_command(seed,Path('/same'),True,grouping_call=7)
        self.assertEqual(grouped[grouped.index('--bddcProlongationGroupingCall')+1],'7')
        off=replay_command(seed,Path('/same'),False,grouping_call=7)
        self.assertEqual({grouped[i-1] for i,(a,b) in enumerate(zip(off,grouped)) if a!=b},
                         {'--bddcProlongationDiagnostics','--bddcProlongationGroupingCall'})
        text='[rank trace] BDDC global residual step 1: abs=1e-3 rel=0.5 target=1e-7 status=not-converged iterations=7\n'
        self.assertEqual(signature(text)['residuals'][0][-1],'7')

    def test_grouping_screen_checks_parity_and_cost(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'bddc-prolongation-grouping.csv'
            rows=[]
            for bits,group in [(16,0),(8,0),*((8,g) for g in [8,16,32,64,128,256,512])]:
                rows.append(dict(solver_call=7,bits=bits,group_size=group,current_precision=int(bits==8 and group==0),
                    relative_difference_vs_actual=0 if bits==8 and group==0 else 0.1,
                    rho_candidate=-1 if group==0 and bits==8 else 1,metadata_bytes=16,
                    modeled_p_bytes=200 if bits==16 or group==512 else 150,codec_valid=1))
            def write():
                with path.open('w') as f:
                    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
            write()
            with contextlib.redirect_stdout(io.StringIO()):
                result=inspect_grouping(root,7)
            self.assertFalse(result[1]['screen_eligible']) # Negative rho despite small packets.
            self.assertTrue(result[2]['screen_eligible'])
            self.assertFalse(result[-1]['screen_eligible']) # No byte improvement.
            rows[1]['relative_difference_vs_actual']=0.1;write()
            with self.assertRaisesRegex(ValueError,'parity'):inspect_grouping(root,7)
            rows[1]['relative_difference_vs_actual']=0;rows.pop();write()
            with self.assertRaisesRegex(ValueError,'Incomplete'):inspect_grouping(root,7)

    def test_expected_failure_is_valid_replay(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo=Path(tmp).resolve();source=repo/'source';model=repo/'emiModel'
            model.write_text('mock');fake=repo/'mpi'
            passed='BDDC global residual step 1: abs=1e-10 rel=1e-8 target=1e-7 status=converged iterations=10\nEMI model completed\n'
            failed='BDDC global residual step 1: abs=1e-3 rel=0.5 target=1e-7 status=not-converged iterations=7\nBDDC step 1 failed after 7 solver calls: solver breakdown. State was not advanced.\n'
            header='solver_call,restriction_pass,true_residual_norm,restriction_error_rel,rho_before_p,rho_after_p,rho_change_p,rho_closure_error,z_before_norm,z_after_norm,p_error_norm,p_error_rel,p_error_energy,actual_z_jump,weight_partition_defect,inactive_global_dofs\n'
            fake.write_text('''#!/usr/bin/env python3
from pathlib import Path
import sys
a=sys.argv[1:]
def v(k):return a[a.index(k)+1]
out=Path(v('--dir'));mode=out.parent.name
with open('calls','a') as f:f.write(mode+'/'+out.name+'\\n')
failed=mode=='full-sd8'
print('''+repr(failed)+''' if failed else '''+repr(passed)+''')
(out/'potential-state-step0.csv').write_text('dof,u\\n0,0\\n')
if v('--bddcProlongationDiagnostics')=='1':
 row='7,1,1,0,2,-1,-3,0,2,1,3,1.5,9,0,0,0\\n' if failed else '10,1,1e-10,0,2,2,0,0,2,2,0,0,0,0,0,0\\n'
 (out/'bddc-prolongation-diagnostics-step1.csv').write_text('''+repr(header)+'''+row)
sys.exit(134 if failed else 0)
''')
            fake.chmod(0o755)
            for mode in ['uncompressed-cg','full-sd8']:
                folder=source/'10cells'/mode;folder.mkdir(parents=True)
                seed=command(model,[str(fake)],Path('/mesh'),2,1,8,'full-r16-p16',folder,
                             potential_diagnostics=True,matched_original_rho=True)
                cmd=make_command(seed,mode,folder)
                (folder/'command.sh').write_text(shlex.join(cmd))
                (folder/'run.log').write_text(failed if mode=='full-sd8' else passed)
            script=Path(__file__).with_name('replay_emi_prolongation_diagnostics.py')
            args=[sys.executable,str(script),'--repo',str(repo),'--campaign',str(source),
                  '--modes','uncompressed-cg','full-sd8','--output',str(repo/'new')]
            result=subprocess.run(args,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            records=json.loads((repo/'new/summary.json').read_text())
            for row in records.values():
                self.assertTrue(row['observer_parity']);self.assertTrue(row['source_outcome_parity'])
            self.assertEqual(records['full-sd8']['on']['return_code'],134)
            self.assertIn('changes positive',records['full-sd8']['diagnostics'][0]['finding'])
            # The optional probe must exist at the requested call, even on expected failure.
            probe_summary=('solver_call,groups,covered_dofs,rho_before_p,rho_after_p,'
                'rho_mean_corrected,rho_float32_mean_corrected,attributed_rho_error,'
                'unattributed_rho_error,max_remaining_mean_error,max_remaining_float32_mean_error,'
                'estimated_payload32_per_call,estimated_payload64_per_call\n'
                '7,1,1,2,-1,2,2,-3,0,0,0,8,16\n')
            probe_groups=('group,dofs,owners,rho_error_with_interior,rho_mean_error_with_interior\n'
                          '0,1,2,-3,-3\n')
            inject=("if '--bddcProlongationConstraintCall' in a and v('--bddcProlongationConstraintCall')=='7':\n"
                    " (out/'bddc-prolongation-constraints-summary.csv').write_text("+repr(probe_summary)+")\n"
                    " (out/'bddc-prolongation-constraints-groups.csv').write_text("+repr(probe_groups)+")\n")
            fake.write_text(fake.read_text().replace('sys.exit(134',inject+'sys.exit(134'))
            result=subprocess.run(args[:-1]+[str(repo/'probe'),'--constraint-call','7'],
                                  capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            probe=json.loads((repo/'probe/summary.json').read_text())
            self.assertEqual(probe['full-sd8']['constraints']['summary']['rho_mean_corrected'],2)
            result=subprocess.run(args[:-1]+[str(repo/'missing'),'--constraint-call','6'],
                                  capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0,'Missing target probe was accepted')
            group_fields=['solver_call','bits','group_size','current_precision',
                          'relative_difference_vs_actual','rho_candidate','metadata_bytes',
                          'modeled_p_bytes','codec_valid']
            buffer=io.StringIO();writer=csv.writer(buffer);writer.writerow(group_fields)
            for bits,group in [(16,0),(8,0),*((8,g) for g in [8,16,32,64,128,256,512])]:
                writer.writerow([7,bits,group,int(bits==8 and group==0),0,
                                 -1 if bits==8 and group==0 else 1,16,200 if bits==16 else 150,1])
            inject=("if '--bddcProlongationGroupingCall' in a and v('--bddcProlongationGroupingCall')=='7':\n"
                    " (out/'bddc-prolongation-grouping.csv').write_text("+repr(buffer.getvalue())+")\n")
            fake.write_text(fake.read_text().replace('sys.exit(134',inject+'sys.exit(134'))
            group_args=args[:args.index('--modes')]+['--modes','full-sd8','--output',
                       str(repo/'grouped'),'--grouping-call','7']
            result=subprocess.run(group_args,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            screened=json.loads((repo/'grouped/summary.json').read_text())['full-sd8']['grouping']
            self.assertTrue(screened[2]['screen_eligible'])
            # Deliberate observer side effect must invalidate the diagnosis.
            fake.write_text(fake.read_text().replace("row='7,1", "(out/'potential-state-step0.csv').write_text('changed snapshot'); row='7,1"))
            result=subprocess.run(args[:-1]+[str(repo/'bad')],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            bad=json.loads((repo/'bad/summary.json').read_text())
            self.assertFalse(bad['uncompressed-cg']['observer_parity'])
            self.assertNotIn('full-sd8',bad)


if __name__=='__main__':unittest.main()
