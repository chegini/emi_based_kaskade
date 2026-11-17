"""Offline checks; never run a simulation."""
import contextlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import run_check as check

class MethodAStudy(unittest.TestCase):
    def test_plan_matrix_and_no_execution(self):
        with patch.object(check.pilot.hb.study, 'execute', side_effect=AssertionError('Launch')):
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                self.assertEqual(check.main(['--bits', '12']), 0)
        plan=json.loads(output.getvalue())
        self.assertEqual(plan['status'], 'PREPARED_NOT_EXECUTED')
        self.assertEqual(len(plan['commands']), 9)
        self.assertEqual((plan['nodes'],plan['ranks'],plan['timesteps']), (2,4,3))
        self.assertEqual({(c['depth'],c['split']) for c in plan['configurations'][1:]},
                         {(d,s) for d in range(1,5) for s in ['on','off']})
        self.assertTrue(all(c['variant']=='A' and c['bits']==12 for c in plan['configurations'][1:]))

    def test_requires_allocation(self):
        with patch.dict(os.environ, {}, clear=True), patch.object(check.pilot.hb.study, 'execute', side_effect=AssertionError('Launch')):
            with self.assertRaisesRegex(ValueError,'user-submitted'):
                check.main(['--bits','12','--execute'])

    def test_failed_raw_skips_all_compressed_cases(self):
        with tempfile.TemporaryDirectory() as tmp, \
             patch.dict(os.environ, dict(SLURM_JOB_ID='offline-test', SLURM_JOB_NUM_NODES='2', SLURM_NTASKS='4')), \
             patch.object(check.base, 'digest', return_value='test-only'), \
             patch.object(check.pilot.hb, 'mesh_inputs', return_value={}), \
             patch.object(check.pilot.hb.study, 'execute', return_value=(1,False,0)) as execute, \
             patch.object(check.base, 'diagnose', side_effect=ValueError('raw failed')), \
             contextlib.redirect_stdout(io.StringIO()):
            out=Path(tmp)/'results'
            self.assertEqual(check.main(['--execute','--bits','12','--output',str(out)]),1)
            rows=json.loads((out/'runs.json').read_text())
            self.assertEqual(execute.call_count,1)
            self.assertEqual([r['run_status'] for r in rows],['FAILED']+['SKIPPED']*8)

    def test_depth_parity_failure_is_not_reported_as_success(self):
        sample=dict(serialized_inter_rank_bytes=100, serialized_inter_node_bytes=40, accepted_iterations=2)
        with tempfile.TemporaryDirectory() as tmp, \
             patch.dict(os.environ, dict(SLURM_JOB_ID='offline-test', SLURM_JOB_NUM_NODES='2', SLURM_NTASKS='4')), \
             patch.object(check.base, 'digest', return_value='test-only'), \
             patch.object(check.pilot.hb, 'mesh_inputs', return_value={}), \
             patch.object(check.pilot.hb.study, 'execute', return_value=(0,False,0)), \
             patch.object(check.base, 'diagnose', return_value=(sample,[])), \
             patch.object(check.base, 'invariants', return_value={'same':True}), \
             patch.object(check.pilot.hb, 'exact_trajectory', side_effect=[False,True,True,True,True,True]), \
             contextlib.redirect_stdout(io.StringIO()):
            out=Path(tmp)/'results'
            self.assertEqual(check.main(['--execute','--bits','12','--output',str(out)]),1)
            rows=json.loads((out/'runs.json').read_text())
            self.assertEqual(rows[2]['run_status'],'FAILED')
            self.assertFalse(rows[2]['depth1_solution_history_parity'])
            self.assertEqual(sum(r['run_status']=='PASS' for r in rows),8)

    def test_production_application_guards_accept_generated_commands(self):
        from test_emi_compression_report import extract, kaskade_include
        repo = Path(__file__).resolve().parents[2]
        source = (repo / 'EMI_model.cpp').read_text()
        options = extract(source, 'struct EmiOptions').replace('Kaskade::BDDC::EDGE', '0') + ';'
        start = source.index('  if (options.quiet != 0')
        end = source.index('  OptionalMpiSession mpiSession', start)
        fields = dict(re.findall(r'\("(\w+)",\s*options\.(\w+)\s*,', source))
        strings = set(re.findall(r'std::string\s+(\w+)', options))
        prefix = ('#include <cmath>\n#include <iostream>\n#include <string>\n'
                  '#include "EMI_partitioning.hh"\n#include "EMI_bddc_precision.hh"\n'
                  '#include "EMI_bddc_packet_packing.hh"\n#include "EMI_bddc_step_schedule.hh"\n' +
                  options + '\nvoid validate(EmiOptions& options){\n' + source[start:end] +
                  'EmiBddc::validateStepSchedule(options,options.maximumNumberOfTimeSteps);\n}\n')
        cfg = check.base.configuration(check.base.arguments([]))
        checks = []
        for case in check.cases(16) + check.cases(12):
            cmd = check.command(cfg, case, Path('/tmp') / case['name'])[1]
            payload = cmd[cmd.index('--mpi'):]
            assignments = []
            for flag, value in zip(payload[::2], payload[1::2]):
                field = fields[flag.removeprefix('--')]
                literal = json.dumps(value)
                if field not in strings:
                    literal = 'std::stod(' + literal + ')'
                assignments.append('options.' + field + '=' + literal + ';')
            checks.append('{EmiOptions options;\n' + '\n'.join(assignments) +
                          '\ntry{validate(options);}catch(std::exception const& e){'
                          'std::cerr<<e.what();return 1;}\n}\n')
        with tempfile.TemporaryDirectory(prefix='pepe-weight-aware-guards-') as tmp:
            cpp = Path(tmp) / 'guards.cpp'
            exe = Path(tmp) / 'guards'
            cpp.write_text(prefix + 'int main(){\n' + ''.join(checks) + '}\n')
            subprocess.run(shlex.split(os.environ.get('CXX', 'c++')) + [
                '-std=c++20', '-I', str(repo), '-I', kaskade_include(repo),
                str(cpp), '-o', str(exe)], check=True, timeout=180)
            subprocess.run([str(exe)], check=True, timeout=30)


if __name__ == '__main__':
    unittest.main()
