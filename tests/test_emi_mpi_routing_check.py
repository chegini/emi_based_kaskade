#!/usr/bin/env python3
"""Verify routing-validation commands and fail-closed topology/accounting parsing."""
import unittest
import json
import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path
from run_emi_mpi_routing_check import arguments, make_command, planned_modes, routing_counts, slurm_container_command
from test_emi_compression_report import extract


class RoutingCheck(unittest.TestCase):
    def test_single_node_raw_backend(self):
        repo = Path(__file__).resolve().parents[1]
        source = (repo/'EMI_model.cpp').read_text()
        options = extract(source, 'struct EmiOptions').replace('Kaskade::BDDC::EDGE', '0')+';'
        helper = extract(source, 'bool selectSingleNodeRawBypass(')
        checks = r'''
int main() {
  for (int nodes: {1,2,4}) for (int bypass: {0,1}) for (int solver: {0,1,2}) {
    EmiOptions o;
    o.bddcCompression=1; o.bddcLocalCompressionBypass=bypass; o.bddcMpiRouting=1;
    o.bddcUseCg=solver!=2; o.bddcFlexibleCg=solver==1; o.bddcFgmres=solver==2;
    o.bddcProlongateConstantSplit=1; o.bddcCgConstantSplitExperiment=solver!=2;
    o.bddcGraphLifting=1; o.bddcGraphLiftingAfterQuantization=1; o.bddcGraphLiftingBothOrders=1;
    o.bddcZigzag=1; o.bddcBitlength=1; o.bddcHuffman=1;
    o.bddcUncompressedSteps=5; o.bddcRestrictionPasses=2; o.dt=0.013;
    bool selected=selectSingleNodeRawBypass(o,nodes);
    if (selected!=(bypass && nodes==1)) return 1;
    if (o.bddcUseCg!=(solver!=2) || o.bddcFlexibleCg!=(solver==1) ||
        o.bddcFgmres!=(solver==2) || o.bddcMpiRouting!=1 ||
        o.bddcRestrictionPasses!=2 || o.dt!=0.013) return 2;
    if (selected) {
      if (o.bddcCompression || o.bddcLocalCompressionBypass || o.bddcProlongateConstantSplit ||
          o.bddcCgConstantSplitExperiment || o.bddcGraphLifting || o.bddcGraphLiftingAfterQuantization ||
          o.bddcGraphLiftingBothOrders || o.bddcZigzag || o.bddcBitlength || o.bddcHuffman ||
          o.bddcUncompressedSteps) return 3;
    } else {
      if (!o.bddcCompression || o.bddcLocalCompressionBypass!=bypass ||
          !o.bddcProlongateConstantSplit || !o.bddcGraphLiftingBothOrders ||
          !o.bddcHuffman || o.bddcUncompressedSteps!=5) return 4;
    }
  }
}
'''
        with tempfile.TemporaryDirectory(prefix='emi-bypass-backend-') as tmp:
            cpp, binary = Path(tmp)/'backend.cpp', Path(tmp)/'backend'
            cpp.write_text('#include <string>\n#include "EMI_bddc_precision.hh"\n'+options+'\n'+helper+'\n'+checks)
            subprocess.run(shlex.split(os.environ.get('CXX', 'c++'))+
                           ['-std=c++20', '-I', str(repo), str(cpp), '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)

    def test_local_bypass_commands(self):
        a = arguments(['--check-local-bypass'])
        self.assertEqual(len(planned_modes(a)), 16)
        for solver in ('cg', 'fcg', 'fgmres'):
            cmd = make_command(a, Path('/tmp/bypass'), ['mpirun'], solver, True, True, bypass=True)
            self.assertEqual(cmd[cmd.index('--bddcLocalCompressionBypass')+1], '1')
            self.assertEqual(cmd[cmd.index('--bddcMpiRouting')+1], '1')
            self.assertEqual(cmd[cmd.index('--bddcCompression')+1], '1')

    def test_tcp_only_container(self):
        a = arguments(['--slurm-container', '/shared/image.sif',
                       '--slurm-btl', 'self,tcp', '--setup-diagnostic-only'])
        cmd = make_command(a, Path('/tmp/diagnostic'), ['mpirun'], 'cg', False, False)
        script = cmd[cmd.index('-c')+1]
        self.assertIn('OMPI_MCA_btl=self,tcp\n', script)
        self.assertNotIn('vader', script)
        self.assertIn('OMPI_MCA_pml=ob1', script)
        self.assertIn('export EMI_MATERIAL_TRACE=1\n', script)
        a.slurm_btl = 'self,vader,tcp'
        cmd = make_command(a, Path('/tmp/diagnostic'), ['mpirun'], 'cg', False, False)
        self.assertIn('OMPI_MCA_btl=self,vader,tcp\n', cmd[cmd.index('-c')+1])

    def test_setup_diagnostic(self):
        a = arguments(['--setup-diagnostic-only'])
        self.assertEqual(planned_modes(a), [('cg-uncompressed-legacy', 'cg', False, False)])
        cmd = make_command(a, Path('/tmp/diagnostic'), ['mpirun'], 'cg', False, False)
        i = cmd.index(str(a.repo/'emiModel'))
        self.assertEqual(cmd[i-2:i], ['env', 'EMI_MATERIAL_TRACE=1'])
        self.assertEqual(cmd[cmd.index('--bddcCompression')+1], '0')
        self.assertEqual(cmd[cmd.index('--bddcMpiRouting')+1], '0')
        a.slurm_container = Path('/shared/image.sif')
        cmd = make_command(a, Path('/tmp/diagnostic'), ['mpirun'], 'cg', False, False)
        self.assertIn('export EMI_MATERIAL_TRACE=1\n', cmd[cmd.index('-c')+1])
        self.assertIn('--time=00:10:00', cmd)
        normal = arguments([])
        self.assertEqual(len(planned_modes(normal)), 13)
        normal.slurm_container = a.slurm_container
        cmd = make_command(normal, Path('/tmp/normal'), ['mpirun'], 'cg', False, False)
        self.assertNotIn('EMI_MATERIAL_TRACE', cmd[cmd.index('-c')+1])

    def test_slurm_container_launch(self):
        a = arguments(['--ranks', '4', '--expected-nodes', '2',
                       '--slurm-container', '/shared/image.sif', '--kaskade', '/shared/kaskade'])
        for direct in (False, True):
            cmd = make_command(a, Path('/shared/results'), ['unused-mpirun'],
                               'fgmres', not direct, not direct, direct)
            self.assertEqual(cmd[0], 'srun')
            self.assertIn('--mpi=pmix_v5', cmd)
            self.assertIn('--nodes='+('1' if direct else '2'), cmd)
            self.assertIn('--ntasks='+('1' if direct else '4'), cmd)
            self.assertIn('--ntasks-per-node='+('1' if direct else '2'), cmd)
            self.assertIn('--distribution=block:block', cmd)
            self.assertIn('--export=ALL', cmd)
            for flag in ('--cleanenv', '-np', '-x', '--allow-run-as-root', 'unused-mpirun'):
                self.assertNotIn(flag, cmd)
            self.assertIn(f'{a.repo.parent}:{a.repo.parent}', cmd)
            self.assertIn('/shared/kaskade:/kaskade7', cmd)
            self.assertIn(str(a.repo/'emiModel'), cmd)
            self.assertEqual(cmd[cmd.index('--bddcMpiRouting')+1], '0' if direct else '1')
        probe = slurm_container_command(a, ['/bin/true'], 2, 2)
        self.assertIn('--mpi=none', probe)
        self.assertIn('--ntasks-per-node=1', probe)

    def test_application_guards(self):
        repo = Path(__file__).resolve().parents[1]
        source = (repo/'EMI_model.cpp').read_text()
        options = extract(source, 'struct EmiOptions').replace('Kaskade::BDDC::EDGE', '0')+';'
        guards = source[source.index('  if (options.quiet != 0'):
                        source.index('  OptionalMpiSession mpiSession')]
        fields = dict(re.findall(r'\("(\w+)",\s*options\.(\w+)\s*,', source))
        strings = set(re.findall(r'std::string\s+(\w+)', options))
        checks = []
        a = arguments([])

        def add(cmd, expected):
            start = cmd.index(str(a.repo/'emiModel'))+1
            values = dict(zip(cmd[start::2], cmd[start+1::2]))
            assignments = []
            for flag, value in values.items():
                field = fields[flag.removeprefix('--')]
                literal = json.dumps(value)
                assignments.append(f'options.{field} = '+
                                   (literal if field in strings else f'std::stod({literal})')+';')
            checks.append('{ EmiOptions options;\n'+'\n'.join(assignments)+
                          '\nbool accepted=true; try { validate(options); } '
                          'catch (std::exception const&) {accepted=false;}\n'+
                          f'if (accepted != {str(expected).lower()}) return 1; }}\n')

        for order in ('before', 'after', 'both'):
            a.transform_order = order
            for solver in ('cg', 'fcg', 'fgmres'):
                for compressed in (False, True):
                    for routed in (False, True):
                        cmd = make_command(a, Path('/tmp/out'), ['mpirun'], solver, compressed, routed)
                        add(cmd, True)
        add(make_command(a, Path('/tmp/out'), ['mpirun'], 'fgmres', False, False, True), True)
        for solver in ('cg', 'fcg', 'fgmres'):
            cmd = make_command(a, Path('/tmp/out'), ['mpirun'], solver, True, True, bypass=True)
            add(cmd, True)
            for flag in ('--bddcCompression', '--bddcMpiRouting'):
                invalid = list(cmd)
                invalid[invalid.index(flag)+1] = '0'
                add(invalid, False)
            invalid = list(cmd)
            invalid[invalid.index('--bddcHistogram')+1] = '1'
            add(invalid, False)
        cmd = make_command(a, Path('/tmp/out'), ['mpirun'], 'fgmres', False, True)
        for flag, value in [('--mpi', '0'), ('--sdc', '1'), ('--bddcMpiRouting', '2')]:
            invalid = list(cmd)
            invalid[invalid.index(flag)+1] = value
            add(invalid, False)
        preamble = ('#include <cmath>\n#include <iostream>\n#include <string>\n'
                    '#include "EMI_bddc_precision.hh"\n#include "EMI_bddc_step_schedule.hh"\n'+options+
                    '\nvoid validate(EmiOptions& options) {\n'+guards+'}\n')
        with tempfile.TemporaryDirectory(prefix='emi-routing-guards-') as tmp:
            cpp, binary = Path(tmp)/'guards.cpp', Path(tmp)/'guards'
            cpp.write_text(preamble+'int main() {\n'+''.join(checks)+'}\n')
            subprocess.run(shlex.split(os.environ.get('CXX', 'c++'))+
                           ['-std=c++20', '-I', str(repo), str(cpp), '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)

    def test_controls(self):
        a = arguments([])
        for solver in ('cg', 'fcg', 'fgmres'):
            for compressed in (False, True):
                for routed in (False, True):
                    cmd = make_command(a, Path('/tmp/test'), ['mpirun'], solver, compressed, routed)
                    def value(flag):
                        self.assertEqual(cmd.count(flag), 1)
                        return cmd[cmd.index(flag)+1]
                    self.assertEqual(value('--bddcMpiRouting'), str(int(routed)))
                    self.assertEqual(value('--bddcHistogram'), '0')
                    self.assertNotIn('--bddcHistogramSteps', cmd)
                    self.assertEqual(value('--maximumNumberOfTimeSteps'), '1')
                    self.assertEqual(value('--potentialDiagnostics'), '1')
                    self.assertEqual(value('--bddcCompression'), str(int(compressed)))
                    self.assertEqual(value('--bddcProlongateConstantSplit'), str(int(compressed)))
                    self.assertEqual(value('--bddcCgConstantSplitExperiment'),
                                     str(int(compressed and solver != 'fgmres')))
        direct = make_command(a, Path('/tmp/test'), ['mpirun'], 'cg', False, True, True)
        self.assertEqual(direct[direct.index('--bddcMpiRouting')+1], '0')
        self.assertEqual(direct[direct.index('-np')+1], '1')

    def test_counter_validation(self):
        header = 'BDDC MPI routing: destination; shared-memory domains=1; rank leaders= 0 0\n'
        lines = [f'BDDC routed payload: phase={p} route={r} bytes={b} packets={n}\n'
                 for p in ('restrict', 'prolongate')
                 for r, b, n in [('same-rank', 64, 2), ('same-node', 32, 1), ('inter-node', 0, 0)]]
        log = header+''.join(lines)
        self.assertEqual(sum(r['bytes'] for r in routing_counts(log, 1)), 192)
        for invalid in (log+lines[0], header+''.join(lines[:-1]),
                        log.replace('inter-node bytes=0', 'inter-node bytes=8'),
                        log.replace('domains=1', 'domains=2')):
            with self.assertRaises(ValueError):
                routing_counts(invalid, 1)
        two_nodes = log.replace('domains=1', 'domains=2')
        with self.assertRaises(ValueError):
            routing_counts(two_nodes, 2)
        two_nodes = two_nodes.replace('inter-node bytes=0 packets=0', 'inter-node bytes=8 packets=1')
        self.assertEqual(sum(r['bytes'] for r in routing_counts(two_nodes, 2)), 208)


if __name__ == '__main__':
    unittest.main()
