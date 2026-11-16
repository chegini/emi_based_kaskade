#!/usr/bin/env python3
"""Native collector, bounded spills, unsigned32 and conditional/pooled H0 tests."""
import csv,math,os,random,shlex,struct,subprocess,sys,tempfile,unittest
from collections import Counter
from pathlib import Path
from summarize_emi_message_entropy import summarize
HERE=Path(__file__).resolve().parent
CPP=r'''
#include "emi_message_entropy.hh"
using namespace Kaskade::BDDC::MessageEntropy;
int main(){
 topology(0,1,{0});
 auto& c=collector();
 c.observe(0,0,4,16,std::vector<unsigned>{0,0,0,0});
 c.observe(0,0,4,16,std::vector<unsigned>{1,1,1,1});
 c.observe(1,0,4,32,std::vector<std::uint32_t>{0,0xffffffffu,0,0xffffffffu});
 c.observe(1,1,4,32,std::vector<std::uint32_t>{0,0xffffffffu,0,0xffffffffu});
 coefficients("x/subdomain0","restrict",0,16,std::vector<unsigned>{100,100,101,101},std::vector<unsigned>{0,0,1,1},std::vector<unsigned>{0,0,1,1},true,true,false);
 coefficients("x/subdomain0","restrict",0,16,std::vector<unsigned>{5,5,5},std::vector<unsigned>{5,5,5},std::vector<unsigned>{3,3,3},false,false,false);
 topology(0,3,{0,1,2});
 coefficients("x/subdomain0","restrict",1,8,std::vector<unsigned>{10,20},std::vector<unsigned>{0,2},std::vector<unsigned>{0,2},true,false,true);
 routed("restrict",1,{0,1,2,3});routed("prolongate",2,{255,255,255});
 for(unsigned j=0;j<10;++j){std::vector<unsigned> v;for(unsigned i=0;i<20;++i)v.push_back(j*20+i);c.observe(1,0,4,16,v);}
 topology(0,3,{0,1,2},{1,0,2,0});
 if(peerRoute("x/subdomain1",0)!=1 || peerRoute("x/subdomain3",2)!=2)return 1;
 bool rejected=false;
 try{peerRoute("x/subdomain0",1);}catch(std::runtime_error const&){rejected=true;}
 if(!rejected)return 2;
}
'''
class EntropyTests(unittest.TestCase):
 def test_native(self):
  with tempfile.TemporaryDirectory(prefix='emi-entropy-test-') as td:
   t=Path(td);src=t/'test.cpp';src.write_text(CPP)
   subprocess.run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++20','-O2','-Wall','-Wextra','-Werror','-I',str(HERE),str(src),'-o',str(t/'test')],check=True)
   env=dict(os.environ,EMI_MESSAGE_ENTROPY_PREFIX=str(t/'counts/entropy'),EMI_ENTROPY_MAX_BINS='16',OMPI_COMM_WORLD_RANK='0')
   subprocess.run([str(t/'test')],env=env,check=True)
   self.assertGreater(len(list((t/'counts').glob('*.bin'))),10)
   rows,ranks=summarize(t/'counts',t/'summary',merge_fanin=2);self.assertEqual(ranks,[0])
   direct,_=summarize(t/'counts',t/'summary-direct')
   self.assertEqual(rows,direct)
   self.assertFalse(list((t/'summary').glob('entropy-merge-*')))
   find=lambda **kw:next(r for r in rows if all(r[k]==v for k,v in kw.items()))
   r=find(phase='restrict',stage='quantized_input',floating=False,encoded=True,source_bits=16)
   self.assertEqual(r['symbols'],8);self.assertEqual(r['messages'],2)
   self.assertAlmostEqual(r['pooled_shannon_bits_per_symbol'],1)
   self.assertAlmostEqual(r['message_weighted_shannon_bits_per_symbol'],0)
   r=find(phase='prolongate',stage='quantized_input',source_bits=32)
   self.assertAlmostEqual(r['message_weighted_shannon_bits_per_symbol'],1)
   r=find(phase='prolongate',stage='quantized_input',source_bits=16)
   self.assertEqual(r['distinct_symbols'],200)
   self.assertAlmostEqual(r['pooled_shannon_bits_per_symbol'],math.log2(200))
   self.assertAlmostEqual(r['message_weighted_shannon_bits_per_symbol'],math.log2(20))
   r=find(phase='restrict',stage='routed_bytes',route='same_node');self.assertAlmostEqual(r['message_weighted_shannon_bits_per_symbol'],2)
   self.assertEqual(r['unit'],'bits/byte')
   r=find(phase='prolongate',stage='routed_bytes');self.assertAlmostEqual(r['message_weighted_shannon_bits_per_symbol'],0)
   for stage in ('quantized_input','encoder_coefficients'):
    self.assertAlmostEqual(find(phase='restrict',stage=stage,floating=True)['message_weighted_shannon_bits_per_symbol'],1)
   self.assertEqual(find(phase='restrict',stage='encoder_coefficients',integer=True)['route'],'same_node')
   self.assertFalse((t/'summary/entropy-counts.sqlite').exists())
 def test_no_sqlite_dependency(self):
  code="import sys; sys.modules['sqlite3']=None; sys.modules['_sqlite3']=None; sys.path.insert(0,sys.argv[1]); import summarize_emi_message_entropy; print('PASS: SQLite unavailable, entropy import succeeds')"
  subprocess.run([sys.executable,'-B','-c',code,str(HERE)],check=True)
 def test_multirank_external_merge(self):
  with tempfile.TemporaryDirectory(prefix='emi-entropy-merge-') as td:
   folder=Path(td)/'counts';folder.mkdir();rng=random.Random(6184527)
   global_counts={0:Counter(),1:Counter()};weighted={0:0.,1:0.};messages={0:0,1:0}
   for part in range(77):
    phase=part%2;counts=Counter(rng.choice([0,0xffffffff,*range(100)]) for _ in range(300))
    global_counts[phase].update(counts);messages[phase]+=1
    weighted[phase]+=300*math.log2(300)-math.fsum(n*math.log2(n) for n in counts.values())
    bucket=(32<<56)|(20<<48)|(phase<<40)|(1<<32)
    with (folder/f'entropy-rank-{part%3}-pid-1-part-{part}.bin').open('wb') as f:
     f.write(b'EMIENT01')
     for symbol,n in sorted(counts.items()):f.write(struct.pack('<QQ',bucket|symbol,n))
   with (folder/'entropy-rank-0-pid-1-totals.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['phase','stage','flags','bits','symbols','messages','weighted_entropy_bits'])
    for phase,c in global_counts.items():w.writerow([phase,1,20,32,sum(c.values()),messages[phase],weighted[phase]])
   for rank in (1,2):
    (folder/f'entropy-rank-{rank}-pid-1-totals.csv').write_text('phase,stage,flags,bits,symbols,messages,weighted_entropy_bits\n')
   rows,ranks=summarize(folder,Path(td)/'out',merge_fanin=4);self.assertEqual(ranks,[0,1,2])
   for row in rows:
    phase=['restrict','prolongate'].index(row['phase']);c=global_counts[phase];n=sum(c.values())
    expected=math.log2(n)-math.fsum(v*math.log2(v) for v in c.values())/n
    self.assertEqual(row['symbols'],n);self.assertEqual(row['distinct_symbols'],len(c))
    self.assertAlmostEqual(row['pooled_shannon_bits_per_symbol'],expected,places=12)
    self.assertAlmostEqual(row['message_weighted_shannon_bits_per_symbol'],weighted[phase]/n,places=12)
   self.assertEqual(len(list(folder.glob('*.bin'))),77)
   self.assertFalse(list((Path(td)/'out').glob('entropy-merge-*')))
   bad=folder/'entropy-rank-0-pid-1-part-bad.bin';bad.write_bytes(b'EMIENT01x')
   with self.assertRaises(ValueError):summarize(folder,Path(td)/'bad')
if __name__=='__main__':unittest.main()
