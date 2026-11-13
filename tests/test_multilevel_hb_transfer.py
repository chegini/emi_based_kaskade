#!/usr/bin/env python3
"""Compile actual production plan construction and setup methods without Dune.

The fixture replaces only interface storage; predictor construction, forward,
inverse and serialized remote plans are extracted verbatim from production.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--kaskade-root',type=Path,required=True)
a=p.parse_args()
source=(a.kaskade_root/'mg/bddcSpaceTransfer.hpp').read_text()
names=['buildTwoLevelSplit','graphLiftingPlan','floatLiftingPlan','floatLift',
       'integerLiftingPlan','exportIntegerLiftingPlan','importIntegerLiftingPlan']
methods=[]
for name in names:
    import re
    match=re.search(r'(?:::|\n  )'+name+r'\(',source)
    assert match,name
    start=source.rfind('  template <int m,',0,match.start())
    brace=source.index('{',match.end())
    pos=brace+1;depth=1
    while depth:
        if source[pos]=='{':depth+=1
        elif source[pos]=='}':depth-=1
        pos+=1
    methods.append(source[start:pos])
fixture=r"""
#include "mg/bddcLiftingPlan.hh"
#include "mg/bddcRecenter.hh"
#include <map>
#include <set>
#include <queue>
#include <numeric>
#include <iostream>
namespace Kaskade::BDDC {
template<int m,class Scalar,class TransmissionScalar,class quantizedScalar,class StorageType>
struct SpaceTransferDataCompression {
  using Plan=IntegerLiftingPlan<Scalar>;
  std::vector<int> globalIndex;
  std::map<int,std::vector<int>> localIfs;
  std::vector<std::vector<int>> interface_graph;
  std::vector<std::vector<Scalar>> interface_graph_weights;
  unsigned graph_lifting_levels=1;
  bool graph_lifting_path_hierarchy=false;
  bool transform_first_recentering=false;
  bool producer_lifting_setup_stale=false, require_producer_lifting_plans=true;
  bool transform_after_quantization_restrict=true,transform_after_quantization_prolongate=true;
  mutable std::map<std::vector<int>,Plan> local_float_plans,local_integer_plans;
  std::map<int,Plan> producer_float_plans,producer_lifting_plans;
  void ensureInterfaceGraph() const {}
  bool floatLifting(bool) const {return true;}
  std::vector<int> filteredLocalIdx(int peer) const {return localIfs.at(peer);}
  void buildTwoLevelSplit(std::vector<int> const&,std::vector<int>&,std::vector<int>&,
       std::vector<std::vector<int>>&,std::vector<std::vector<Scalar>>&) const;
  Plan graphLiftingPlan(std::vector<int> const&) const;
  Plan const& floatLiftingPlan(std::vector<int> const&,int=-1) const;
  std::vector<Scalar> floatLift(std::vector<Scalar> const&,std::vector<int> const&,bool,int=-1) const;
  Plan integerLiftingPlan(std::vector<int> const&,int=-1) const;
  std::vector<unsigned char> exportIntegerLiftingPlan(int) const;
  void importIntegerLiftingPlan(int,std::vector<unsigned char> const&);
};
"""+'\n'.join(methods)+r"""
}
using T=Kaskade::BDDC::SpaceTransferDataCompression<1,double,double,std::uint32_t,std::uint8_t>;
void check(bool ok) {if(!ok) throw std::runtime_error("Production HB plan test failed");}
T make(unsigned depth) {
 T t;t.globalIndex.resize(31);std::iota(t.globalIndex.begin(),t.globalIndex.end(),0);
 t.interface_graph.resize(31);t.interface_graph_weights.resize(31);
 for(int i=1;i<31;++i) {
  t.interface_graph[i-1].push_back(i);t.interface_graph[i].push_back(i-1);
  t.interface_graph_weights[i-1].push_back(1);t.interface_graph_weights[i].push_back(1);
 }
 t.graph_lifting_levels=depth;t.localIfs[1]=t.globalIndex;
 t.localIfs[2]=std::vector<int>(t.globalIndex.begin(),t.globalIndex.begin()+10);
 return t;
}
int main() {
 for(unsigned depth:{0u,1u,2u,4u}) for(bool irregular:{false,true}) {
  auto sender=make(depth),receiver=make(1);
  if(irregular) for(int i=2;i<31;i+=3) {
   sender.interface_graph[i-2].push_back(i);sender.interface_graph[i].push_back(i-2);
   sender.interface_graph_weights[i-2].push_back(1);sender.interface_graph_weights[i].push_back(1);
  }
  // Deliberately different receiver topology: it must use imported producer plans.
  for(auto& row:receiver.interface_graph) row.clear();
  std::vector<double> values(31);for(int i=0;i<31;++i) values[i]=0.125*i*i-3;
  auto coeff=sender.floatLift(values,sender.globalIndex,false);
  auto own=sender.floatLift(coeff,sender.globalIndex,true);
  for(int i=0;i<31;++i) check(std::abs(own[i]-values[i])<1e-12);
  for(int peer:{1,2}) {
   receiver.localIfs[0]=sender.localIfs[peer];
   receiver.importIntegerLiftingPlan(0,sender.exportIntegerLiftingPlan(peer));
   auto idx=sender.localIfs[peer];std::vector<double> part,original;
   for(int i:idx) {part.push_back(coeff[i]);original.push_back(values[i]);}
   // Owner grouping means one union analysis restricts exactly to each stream.
   check(part==sender.floatLift(original,idx,false));
   auto restored=receiver.floatLift(part,idx,true,0);
   for(std::size_t i=0;i<restored.size();++i) check(std::abs(restored[i]-original[i])<1e-12);
   std::vector<std::uint32_t> u(idx.size());
   for(std::size_t i=0;i<u.size();++i) u[i]=i%2 ? 0xffffffffu : i;
   auto plan=sender.integerLiftingPlan(idx);
   check(receiver.integerLiftingPlan(idx,0).inverse<1>(plan.forward<1>(u,true),true)==u);
  }
  if(depth==0 && !irregular) check(sender.integerLiftingPlan(sender.globalIndex).coarse.size()==2);
  sender.producer_lifting_setup_stale=true;
  bool rejected=false;try {sender.floatLiftingPlan(sender.globalIndex);} catch(std::runtime_error const&) {rejected=true;}
  check(rejected);
 }
 // Exercise the production selector and owner-group mapping for corrected paths.
 for(unsigned depth:{0u,1u,2u,3u}) {
  auto t=make(depth);t.graph_lifting_path_hierarchy=true;
  std::vector<double> x(31);for(int i=0;i<31;++i)x[i]=2+.125*i;
  auto c=t.floatLift(x,t.globalIndex,false);
  auto p=t.floatLiftingPlan(t.globalIndex);
  for(int f:p.fine)check(std::abs(c[f])<1e-12);
  check(t.floatLift(c,t.globalIndex,true)==x);
  auto idx=t.localIfs[1];auto plan=t.integerLiftingPlan(idx);
  std::vector<std::uint32_t> u(idx.size());for(std::size_t i=0;i<u.size();++i)u[i]=i%2?0xffffffffu:i;
  check(plan.inverse<1>(plan.forward<1>(u,true),true)==u);
 }
 std::cout<<"PASS: production hierarchy, owner-group restriction, remote plans, integer/floating inverse and stale setup rejection.\n";
}
"""
header=(a.kaskade_root/'mg/bddcSpaceTransfer.hh').read_text()
start=header.index('    void setGraphLiftingLevels(')
brace=header.index('{',start);pos=brace+1;depth=1
while depth:
    if header[pos]=='{':depth+=1
    elif header[pos]=='}':depth-=1
    pos+=1
fixture=fixture.replace('  void ensureInterfaceGraph() const {}',header[start:pos]+'\n  void ensureInterfaceGraph() const {}')
fixture=fixture.replace('  sender.producer_lifting_setup_stale=true;',
    '  sender.setGraphLiftingLevels(depth==1 ? 0 : 1); check(sender.local_float_plans.empty() && sender.local_integer_plans.empty());')

with tempfile.TemporaryDirectory(prefix='hb-transfer-') as folder:
    src=Path(folder)/'test.cpp';exe=Path(folder)/'test'
    src.write_text(fixture)
    subprocess.run(shlex.split(os.environ.get('CXX','c++'))+['-std=c++20','-O2','-Wall','-Wextra','-Werror','-I',str(a.kaskade_root),str(src),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
