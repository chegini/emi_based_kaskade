import csv
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


class Lifecycle(unittest.TestCase):
    def test_production_recorder(self):
        repo=Path(__file__).resolve().parents[1]
        k=Path(os.environ.get('KASKADE7',str(repo.parent/'kaskade7_new')))
        if not (k/'mg/bddcTransferTiming.hh').is_file():
            k=Path('/Users/fchegini/Desktop/kaskade1/kaskade7_new')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cpp=root/'test.cpp';exe=root/'test'
            cpp.write_text(r'''#include "EMI_bddc_time_study.hh"
int main(int argc,char** argv) {
 if(argc!=2)return 1;
 using namespace Kaskade::BDDC::TransferTiming;
 using P=std::tuple<std::uint64_t,std::uint64_t,std::array<double,5>,std::uint64_t>;
 EmiBddc::TimeStudy off(false);off.write(std::string(argv[1])+"/disabled",0);
 EmiBddc::TimeStudy run(true);
 {Scope setup(Stage::PlanSetup,1);}
 run.setupComplete(P{});
 for(int step=1;step<=2;++step) {
  run.beginStep();
  for(int call=0;call<3;++call) {
   run.beginCall();
   {Context phase(true);Scope q(Stage::Quantize,1);Scope net(Stage::PacketExchange,1);}
   run.endCall(step);
  }
  run.endStep(step,3,1e-8,P{100u*step,40u*step,{.01*step,.02*step,.03*step,.04*step,.05*step},20u*step},{60u*step,24u*step,40u*step,16u*step});
 }
 run.write(argv[1],0);
}''')
            subprocess.run(shlex.split(os.environ.get('CXX','c++'))+[
                '-std=c++20','-O2','-Wall','-Wextra','-Werror','-I',str(repo),
                '-I',str(k),str(cpp),'-o',str(exe)],check=True)
            for mode in ['0','2']:
                target=root/('mode'+mode)
                subprocess.run([str(exe),str(target)],check=True,
                    env=dict(os.environ,KASKADE_TRANSFER_TIMING=mode))
                with (target/'hb-time-rank-0-total.csv').open() as f:total=list(csv.DictReader(f))
                self.assertGreater(int(total[0]['process_peak_rss_bytes']),0)
                with (target/'hb-time-rank-0.csv').open() as f:rows=list(csv.DictReader(f))
                self.assertEqual([int(r['step']) for r in rows],[0,1,2])
                self.assertEqual([int(r['packed_bytes']) for r in rows],[0,40,40])
                self.assertEqual([int(r['calls']) for r in rows],[0,3,3])
                self.assertEqual([int(r['restrict_packed_bytes']) for r in rows],[0,24,24])
                self.assertEqual([int(r['prolongate_packed_bytes']) for r in rows],[0,16,16])
                self.assertFalse((target/'disabled').exists())
                with (target/'hb-time-rank-0-stages.csv').open() as f:stages=list(csv.DictReader(f))
                if mode=='0':self.assertEqual(stages,[])
                else:
                    q=[r for r in stages if r['stage']=='quantize']
                    self.assertEqual(sum(int(r['calls']) for r in q if r['call']=='0'),6)
                    self.assertEqual(sum(int(r['calls']) for r in q if r['call']!='0'),4)
                    self.assertEqual(len([r for r in stages if r['stage']=='plan_setup']),1)

if __name__=='__main__':unittest.main()
