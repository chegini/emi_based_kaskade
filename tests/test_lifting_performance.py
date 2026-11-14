#!/usr/bin/env python3
"""Compile actual predictor/cache/transform methods and compare all optimization masks."""
from pathlib import Path
import argparse,subprocess,os,hashlib
p=argparse.ArgumentParser();p.add_argument('--kaskade',type=Path,required=True);p.add_argument('--fixture',type=Path);p.add_argument('--test',type=Path,default=Path(__file__).with_name('test_multilevel_hb_transfer.py'));p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
# Reuse a generated production-method fixture; its synthetic main is discarded.
from lifting_performance_fixture import fixture
s=a.fixture.read_text().split('#include <chrono>')[0] if a.fixture else fixture(a.kaskade,a.test)
s+=r'''
#include <fstream>
#include <random>
void bytes(std::ofstream& out,auto const& v){
 auto n=v.size();out.write(reinterpret_cast<char const*>(&n),sizeof(n));
 out.write(reinterpret_cast<char const*>(v.data()),v.size()*sizeof(v[0]));
}
int main(int argc,char** argv){
 if(argc!=2)return 2;
 std::ofstream out(argv[1],std::ios::binary);std::mt19937 rng(61114);
 for(int n:{0,1,2,3,7,31,64,127,257})for(int graph=0;graph<5;++graph)
 for(unsigned depth:{0u,1u,2u,4u}){
  T t;t.globalIndex.resize(n);std::iota(t.globalIndex.begin(),t.globalIndex.end(),0);
  t.interface_graph.resize(n);t.interface_graph_weights.resize(n);
  t.localIfs[1]=t.globalIndex;
  t.localIfs[2]=std::vector<int>(t.globalIndex.begin(),t.globalIndex.begin()+n/2);
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
   bool edge=graph==0?j==i+1:graph==1?(j==i+1 || j==i+8):graph==2?true:
             graph==3?(rng()%19==0):(j==i+1 && i%7!=0);
   if(!edge)continue;
   double w=.125+double(rng()%31)/8;
   t.interface_graph[i].push_back(j);t.interface_graph[j].push_back(i);
   t.interface_graph_weights[i].push_back(w);t.interface_graph_weights[j].push_back(w);
  }
  t.graph_lifting_levels=depth;t.reuse_lifting_plans=true;t.prepareIntegerLiftingPlans();
  for(auto const& idx:{t.globalIndex,t.localIfs[2]}){
   auto old=t.integerLiftingPlan(idx);auto const& ref=t.cachedIntegerLiftingPlan(idx);
   check(old.pack()==ref.pack());bytes(out,ref.pack());
   std::vector<double> u(idx.size());std::vector<std::uint32_t> q(idx.size());
   for(std::size_t i=0;i<u.size();++i){u[i]=double(rng()%9999)/32;q[i]=rng();}
   auto coeff=t.floatLift(u,idx,false);bytes(out,coeff);
   auto inverse=t.floatLift(coeff,idx,true);bytes(out,inverse);
   for(std::size_t i=0;i<u.size();++i)check(std::abs(inverse[i]-u[i])<1e-9);
   auto v=t.applyIntegerGraphLiftingTransform(q,idx,true);bytes(out,v);
   check(t.inverseIntegerGraphLiftingTransform(v,idx,true)==q);
   check(v==old.forward<1>(q,true));
   std::vector<std::uint8_t> q8(idx.size(),255);std::vector<std::uint16_t> q16(idx.size(),65535);
   check(ref.inverse<1>(ref.forward<1>(q8,true),true)==q8);
   check(ref.inverse<1>(ref.forward<1>(q16,true),true)==q16);
   std::vector<std::uint32_t> q2(idx.size()*2,0xffffffffu);
   check(ref.inverse<2>(ref.forward<2>(q2,true),true)==q2);
  }
  t.setGraphLiftingLevels(depth==1?4:1);
  check(t.local_integer_plans.empty() && t.local_float_plans.empty() && t.local_group_plans.empty());
  bool rejected=false;try{t.cachedIntegerLiftingPlan(t.globalIndex);}catch(std::runtime_error const&){rejected=true;}
  check(rejected);t.prepareIntegerLiftingPlans();
  check(t.cachedIntegerLiftingPlan(t.globalIndex).nodeCount==static_cast<std::size_t>(n));
 }
 std::vector<std::vector<int>> tooLarge(24000);
 Kaskade::BDDC::LiftingPerformance::FillBitmap fallback(tooLarge);check(!fallback.active());
}
'''
cpp=a.output/'regression.cpp';exe=a.output/'regression';cpp.write_text(s)
subprocess.run(['c++','-std=c++20','-O2','-Wall','-Wextra','-Werror','-I',str(a.kaskade),str(cpp),'-o',str(exe)],check=True)
hashes=[]
for mask in range(8):
 path=a.output/f'mask-{mask}.bin'
 subprocess.run([str(exe),str(path)],env=dict(os.environ,KASKADE_LIFTING_OPTIMIZATION=str(mask),KASKADE_TRANSFER_TIMING='0',KASKADE_LIFTING_PROFILE_PREFIX=''),check=True)
 hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
assert len(set(hashes))==1,hashes
print('PASS: 180 graph/depth cases x 8 masks; exact serialized plans, floating coefficients/inverses and integer symbols; 8/16/32-bit inverses, vector components, invalidation, bounded-memory fallback.')
print(hashes[0])
