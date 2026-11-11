#!/usr/bin/env python3
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

from run_emi_histogram_start_check import arguments, command, collect_bddc_run, snapshots_equal, summarize_histograms
from test_emi_compression_report import extract
from test_emi_scaling_smoke import log as solver_log


class HistogramScheduleTests(unittest.TestCase):
    def test_solver_rebuild_at_compression_switch(self):
        a = arguments([])
        header = 'BDDC coarse interfaces: corners=1, edges=42, faces=32\n'
        raw = 'BDDC transfer phase: steps 1..2; uncompressed; state retained; physical time=0\n'
        compressed = 'BDDC transfer phase: steps 3..3; compressed; state retained; physical time=0.02\n'
        text = raw+solver_log(steps=3,ranks=2)+compressed+header
        for name in ('scheduled','scheduled-histograms'):
            self.assertEqual(collect_bddc_run(a,name,text,0)['errors'], [])
        for bad in (text.replace(compressed,''), text+header,
                    text.replace('steps 3..3','steps 2..3'),
                    text.replace(compressed+header,compressed+header.replace('edges=42','edges=43'))):
            self.assertEqual(collect_bddc_run(a,'scheduled',bad,0)['status'], 'failed')
        plain = solver_log(steps=3,ranks=2)
        self.assertEqual(collect_bddc_run(a,'uncompressed',plain,0)['errors'], [])
        all_raw = raw.replace('1..2','1..3')+plain
        self.assertEqual(collect_bddc_run(a,'scheduled-all-raw',all_raw,0)['errors'], [])
        a.uncompressed_steps=0
        self.assertEqual(collect_bddc_run(a,'scheduled',plain,0)['errors'], [])

    def test_generated_commands_pass_application_validation(self):
        # Compile the actual early option guards, without the FE/MPI stack.
        # A planned command alone cannot establish that emiModel accepts it.
        repo = Path(__file__).resolve().parents[1]
        source = (repo/'EMI_model.cpp').read_text()
        option_struct = extract(source, 'struct EmiOptions')+';'
        # This unused threshold default is the only dependency on a Kaskade enum.
        option_struct = option_struct.replace('Kaskade::BDDC::EDGE', '0')
        guard_start = source.index('  if (options.quiet != 0')
        guard_end = source.index('  OptionalMpiSession mpiSession', guard_start)
        flags = dict(re.findall(r'\("(\w+)",\s*options\.(\w+)\s*,', source))
        string_fields = set(re.findall(r'std::string\s+(\w+)', option_struct))
        preamble = ('#include <cmath>\n#include <iostream>\n#include <string>\n'
                    '#include "EMI_bddc_precision.hh"\n#include "EMI_bddc_step_schedule.hh"\n'+
                    option_struct+'\nvoid validate(EmiOptions& options) {\n'+
                    source[guard_start:guard_end]+
                    'EmiBddc::validateStepSchedule(options, options.maximumNumberOfTimeSteps);\n}\n')
        checks = []
        def add(values, label, expected=True):
            assignments = []
            for flag, value in values.items():
                field = flags[flag.removeprefix('--')]
                literal = json.dumps(value)
                assignments.append(f'options.{field} = '+
                    (literal if field in string_fields else f'std::stod({literal})')+';')
            checks.append('{ EmiOptions options;\n'+'\n'.join(assignments)+
                '\nbool accepted=true; try { validate(options); } catch (std::exception const& e) {'+
                'accepted=false; '+('std::cerr << e.what() << "\\n";' if expected else '')+'}\n'+
                f'if (accepted != {str(expected).lower()}) {{ std::cerr << '+
                json.dumps(label)+' << "\\n"; return 1; }\n}')
        for solver, split in [('cg','off'), ('fcg','off'), ('fcg','on'), ('fgmres','off'), ('fgmres','on')]:
            for order in ('before', 'after', 'both'):
                a = arguments(['--solver',solver,'--split',split,'--transform-order',order])
                for name in ('direct','uncompressed','scheduled-all-raw','scheduled','scheduled-histograms'):
                    cmd = command(a,name,Path('/tmp/out'),['mpirun'])
                    start = cmd.index(str(a.repo/'emiModel'))+1
                    values = dict(zip(cmd[start::2],cmd[start+1::2]))
                    add(values, f'{solver}/{split}/{order}/{name}')
                    if solver=='fgmres' and split=='on' and order=='before':
                        # Reproduce each invalid direct/BDDC flag from the failed runner.
                        bad_flags = ({'--bddcRelativeTolerance':'1e-7', '--bddcMaterialAssembly':'1',
                                      '--bddcInterfaceAssembly':'2'} if name=='direct' else
                                     {'--directBackend':'mumps-spd'})
                        for flag, value in bad_flags.items():
                            add(dict(values, **{flag:value}), name+' rejects '+flag, expected=False)
        with tempfile.TemporaryDirectory() as tmp:
            cpp, exe = Path(tmp)/'guards.cpp', Path(tmp)/'guards'
            cpp.write_text(preamble+'\nint main() {\n'+'\n'.join(checks)+'\n}\n')
            subprocess.run(shlex.split(os.environ.get('CXX','c++'))+
                           ['-std=c++20','-I',str(repo),str(cpp),'-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)

    def test_command_modes(self):
        a = arguments(['--steps', '6', '--uncompressed-steps', '5', '--histogram-steps', '6'])
        a.repo = Path('/BDDC/example')
        for name, compressed, warmup, capture in [
                ('direct', '0', '0', False), ('uncompressed', '0', '0', False),
                ('scheduled-all-raw', '1', '6', False), ('scheduled', '1', '5', False),
                ('scheduled-histograms', '1', '5', True)]:
            cmd = command(a, name, Path('/tmp/out'), ['mpirun'])
            start = cmd.index(str(a.repo/'emiModel'))+1
            values = dict(zip(cmd[start::2], cmd[start+1::2]))
            self.assertEqual(values['--bddcCompression'], compressed)
            self.assertEqual(values['--bddcUncompressedSteps'], warmup)
            self.assertEqual(values['--maximumNumberOfTimeSteps'], '6')
            self.assertEqual(values['--finalTime'], '0.06')
            self.assertEqual(values.get('--bddcHistogramSteps'), '6' if capture else None)
            self.assertEqual(values['--potentialDiagnostics'], '1')
            self.assertEqual(values['--vtk'], '0')

    def test_explicit_transform_orders(self):
        for order, floating, integer, both in [('before','1','0','0'), ('after','0','1','0'), ('both','1','1','1')]:
            a = arguments(['--transform-order', order])
            cmd = command(a, 'scheduled-histograms', Path('/tmp/out'), ['mpirun'])
            for flag, expected in [('bddcGraphLifting',floating), ('bddcGraphLiftingAfterQuantization',integer),
                                   ('bddcGraphLiftingBothOrders',both)]:
                self.assertEqual(cmd[cmd.index('--'+flag)+1], expected)

    def test_invalid_schedules(self):
        for args in [ ['--steps','0'], ['--uncompressed-steps','-1'], ['--uncompressed-steps','3'],
                      ['--histogram-steps','2'], ['--histogram-steps','4'], ['--histogram-steps','3','3'],
                      ['--solver','cg'], ['--solver','fcg','--steps','101'], ['--steps','11'], ['--ranks','0'] ]:
            with self.subTest(args=args), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                arguments(args)

    def test_fcg_split_commands_at_both_onsets(self):
        for steps in (1,6):
            a=arguments(['--solver','fcg','--split','on','--steps',str(steps),
                         '--uncompressed-steps',str(steps-1),'--physical-values'])
            for name in ('direct','uncompressed','scheduled-all-raw','scheduled','scheduled-histograms'):
                cmd=command(a,name,Path('/tmp/out'),['mpirun'])
                def value(flag): return cmd[cmd.index('--'+flag)+1]
                experiment='0' if name in ('direct','uncompressed') else '1'
                self.assertEqual(value('bddcCgConstantSplitExperiment'),experiment)
                self.assertEqual(value('bddcProlongateConstantSplit'),experiment)
                self.assertEqual(value('bddcFgmres'),'0')
                self.assertEqual(value('bddcFlexibleCg'),'0' if name=='direct' else '1')

    def test_later_steps_without_split(self):
        a = arguments(['--steps','20','--uncompressed-steps','5','--histogram-steps','20','6','10',
                       '--solver','fcg','--split','off'])
        self.assertEqual(a.histogram_steps, [6,10,20])

    def test_observer_changes_are_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            left, right = Path(tmp)/'left', Path(tmp)/'right'
            left.mkdir(); right.mkdir()
            for folder in (left, right):
                (folder/'potential-state-step0.csv').write_text('dof,u\n0,1\n')
            snapshots_equal(left, right, [0])
            (right/'potential-state-step0.csv').write_text('dof,u\n0,1.0000000000000002\n')
            with self.assertRaises(ValueError):
                snapshots_equal(left, right, [0])

    @staticmethod
    def packets(root):
        step = root/'symbol-histograms/step3'
        subdomain = step/'subdomain0'
        subdomain.mkdir(parents=True)
        (step/'metadata.txt').write_text('physical_step=3\nstart_time=0.02\n')
        for phase in ('restrict', 'prolongate'):
            (subdomain/f'{phase}-exchange1-neighbor1.csv').write_text(
                'stage,role,symbol,count,bits,huffman_encoded,transform_order\n'+
                ''.join(f'{s},coefficient,0,2,16,1,before\n{s},coefficient,1,2,16,1,before\n'
                        for s in ('quantized_input','raw','unsigned','bitlength')))

    def test_histogram_population_entropy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.packets(root)
            metrics = summarize_histograms(root, [3])
            self.assertEqual(len(metrics), 8)
            for row in metrics:
                self.assertEqual(row['symbols'], 4)
                self.assertEqual(row['empirical_entropy_bits_per_symbol'], 1)
                self.assertEqual(row['fallback_packet_count'], 0)
            self.assertTrue((root/'histogram-counts.csv').is_file())

    def test_wrong_step_population_and_physical_time(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.packets(root)
            with self.assertRaises(ValueError):
                summarize_histograms(root, [2])
            metadata = root/'symbol-histograms/step3/metadata.txt'
            metadata.write_text('physical_step=3\nstart_time=0\n')
            with self.assertRaisesRegex(ValueError, 'time reset'):
                summarize_histograms(root, [3])

    def test_missing_symbols(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.packets(root)
            packet = root/'symbol-histograms/step3/subdomain0/restrict-exchange1-neighbor1.csv'
            packet.write_text(packet.read_text().replace('raw,coefficient,1,2,', 'raw,coefficient,1,3,'))
            with self.assertRaisesRegex(ValueError, 'populations disagree'):
                summarize_histograms(root, [3])


if __name__ == '__main__':
    unittest.main()
