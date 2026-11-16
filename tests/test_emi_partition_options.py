"""Compile actual EMI option guards against generated smoke/pilot commands.

No MPI launch or PDE simulation. Negative regressions reproduce incompatible
diagnostic flags rejected before mesh loading by split and fixed-HB guards.
"""
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

import run_emi_partition_pilot as pilot
import run_emi_partition_smoke as smoke
from test_emi_compression_report import extract, kaskade_include


class ApplicationOptions(unittest.TestCase):
    def test_complete_generated_options(self):
        repo=Path(__file__).resolve().parents[1]
        source=(repo/'EMI_model.cpp').read_text()
        options=extract(source,'struct EmiOptions').replace('Kaskade::BDDC::EDGE','0')+';'
        start=source.index('  if (options.quiet != 0')
        end=source.index('  OptionalMpiSession mpiSession',start)
        fields=dict(re.findall(r'\("(\w+)",\s*options\.(\w+)\s*,',source))
        strings=set(re.findall(r'std::string\s+(\w+)',options))
        prefix='''#include <cmath>
#include <iostream>
#include <string>
#include "EMI_partitioning.hh"
#include "EMI_bddc_precision.hh"
#include "EMI_bddc_packet_packing.hh"
#include "EMI_bddc_step_schedule.hh"
'''+options+'\nvoid validate(EmiOptions& options) {\n'+source[start:end]+\
            'EmiBddc::validateStepSchedule(options,options.maximumNumberOfTimeSteps);\n}\n'
        with tempfile.TemporaryDirectory(prefix='emi-partition-option-guards-') as tmp:
            root=Path(tmp);specifications=[]
            a=smoke.arguments(['--output',str(root/'smoke')])
            for case in smoke.cases():
                cfg,cmd=smoke.command(a,case,root/case['name'])
                specifications.append((case['name'],cmd,True,''))
                bad=list(cmd);bad[bad.index('--bddcCgDiagnostics')+1]='1'
                specifications.append((case['name']+'-old-conflicting-observer',bad,False,
                                       'diagnostics/corrections disabled'))
            a=pilot.arguments(['--output',str(root/'pilot')])
            raw=smoke.arguments(['--uncompressed','--output',str(root/'raw')])
            for case in smoke.cases(True):
                cfg,cmd=smoke.command(raw,case,root/case['name'])
                specifications.append((case['name']+'-one-step-raw',cmd,True,''))
            for case in pilot.cases():
                for validation in [False,True]:
                    cfg,cmd=pilot.command(a,case,root/case['name'],validation)
                    specifications.append((case['name']+'-validation-'+str(validation),cmd,True,''))
            case=next(c for c in pilot.cases() if c['name']=='P2-B-recenter')
            cfg,bad=pilot.command(a,case,root/'unsplit-HB-conflict',True)
            for flag,value in [('bddcCgDiagnostics','1'),('bddcProlongateConstantSplit','0'),
                               ('bddcCgConstantSplitExperiment','0')]:
                bad[bad.index('--'+flag)+1]=value
            specifications.append(('unsplit-HB-conflicting-observer',bad,False,
                                   'Fixed-topology HB study requires'))
            valid=[(label,cmd) for label,cmd,expected,reason in specifications if expected]
            for label,cmd in valid:
                specifications.append((label+'-byte-entropy',cmd+['--bddcSerializedByteEntropy','1'],True,''))
            for value in ['-1','2']:
                specifications.append(('invalid-byte-entropy-'+value,
                    valid[0][1]+['--bddcSerializedByteEntropy',value],False,'Serialized byte entropy requires'))
            checks=[]
            for label,cmd,expected,reason in specifications:
                payload=cmd[cmd.index('--mpi'):]
                self.assertEqual(len(payload)%2,0)
                assignments=[]
                for flag,value in zip(payload[::2],payload[1::2]):
                    field=fields[flag.removeprefix('--')]
                    literal=json.dumps(value)
                    assignments.append('options.'+field+' = '+
                        (literal if field in strings else 'std::stod('+literal+')')+';')
                checks.append('{ EmiOptions options;\n'+'\n'.join(assignments)+\
                    '\nbool accepted=true;std::string error;try{validate(options);}'
                    'catch(std::exception const& e){accepted=false;error=e.what();}\n'+
                    'if(accepted!='+str(expected).lower()+') {'+
                    'std::cerr<<'+json.dumps(label)+"<<\": \"<<error<<'\\n';return 1;}\n"+
                    ('if(error.find('+json.dumps(reason)+')==std::string::npos){std::cerr<<'+
                     json.dumps(label)+"<<\": \"<<error<<'\\n';return 2;}\n" if reason else '')+
                    '}\n')
            cpp=root/'guards.cpp';exe=root/'guards'
            cpp.write_text(prefix+'int main(){\n'+''.join(checks)+'}\n')
            subprocess.run(shlex.split(os.environ.get('CXX','c++'))+[
                '-std=c++20','-I',str(repo),'-I',kaskade_include(repo),str(cpp),'-o',str(exe)],
                check=True,timeout=180)
            subprocess.run([str(exe)],check=True,timeout=30)
            print('PASS: '+str(len(specifications))+' complete generated option sets; '
                  'old diagnostic conflict rejected, corrected smoke/pilot accepted')


if __name__=='__main__':unittest.main()
