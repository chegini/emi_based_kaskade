"""Small counter/export checks; no MPI process or EMI simulation is launched."""
import csv
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


class Report(unittest.TestCase):
    def test_interval_exports(self):
        repo=Path(__file__).resolve().parents[1]
        k=Path(os.environ['KASKADE7'])
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cpp=root/'report.cpp';exe=root/'report'
            cpp.write_text(r'''#include "EMI_partition_report.hh"
int main(int argc,char** argv) {
 if(argc!=2)return 1;
 using namespace Kaskade::BDDC;
 using namespace TransferTiming;
 EmiBddc::Partitioning::TransportReport off(false);off.write(std::string(argv[1])+"/disabled",0);
 EmiBddc::Partitioning::TransportReport report(true);
 for(int step=1;step<=2;++step){
  report.beginStep();
  MpiTransportTraffic traffic;
  auto& r=traffic.phase[0]; r.applications=step;
  if(step==1){r.route[1].add(100,80,40,true);r.route[2].add(200,160,60,true);}
  else r.route[1].add(10,8,4,true);
  traffic.phase[1].applications=step+1;
  {Context c(true);Scope compression(Stage::Compression,1);Scope q(Stage::Quantize,2);}
  report.endStep(step,step+1,step,"PASS",traffic,.01);
 }
 report.write(argv[1],0);
}''')
            subprocess.run(shlex.split(os.environ.get('CXX','c++'))+[
                '-std=c++20','-O2','-Wall','-Wextra','-Werror','-I',str(repo),
                '-I',str(k),str(cpp),'-o',str(exe)],check=True)
            for mode in ['0','2']:
                dest=root/mode
                subprocess.run([str(exe),str(dest)],env=dict(os.environ,KASKADE_TRANSFER_TIMING=mode),check=True)
                def read(name):
                    with (dest/f'partition-{name}-rank-0.csv').open() as f:return list(csv.DictReader(f))
                rows=read('transport');self.assertEqual(len(rows),12)
                r=next(r for r in rows if r['step']=='1' and r['phase']=='restriction' and r['scope']=='inter_rank')
                self.assertEqual((r['sent_bytes'],r['mpi_messages'],r['mean_sent_bytes']),('100','2','50'))
                r=next(r for r in rows if r['step']=='2' and r['phase']=='restriction' and r['scope']=='same_node')
                self.assertEqual(r['max_sent_bytes'],'4')
                empty=next(r for r in rows if r['step']=='2' and r['scope']=='inter_node')
                self.assertEqual(empty['mean_sent_bytes'],'');self.assertEqual(empty['max_sent_bytes'],'')
                self.assertEqual(read('timing')[0]['quantize_s']=='',mode=='0')
                self.assertEqual(read('timing')[0]['compression_s']=='',mode=='0')
                self.assertEqual([r['accepted_updates'] for r in read('steps')],['1','2'])
                self.assertFalse((dest/'disabled').exists())


if __name__=='__main__':unittest.main()
