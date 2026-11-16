#!/usr/bin/env python3
"""Compile production owned-storage setup/R/P exchanges; no EMI mesh or simulation.

The fixture supplies small algebra/storage only. Constructor, averaging-weight
exchange, packet packing, routed exchange and receiver dispatch are production.
Optional --baseline-root builds the saved original and compares P2 signatures.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

import test_bddc_mpi_exchange as base


def constructor(source):
    start=source.index('  template <class Subdomain, class Storage>\n  SharedMemoryDomain<Subdomain,Storage>::SharedMemoryDomain(')
    brace=source.index('{',start)
    pos,depth=brace+1,1
    while depth:
        depth+=(source[pos]=='{')-(source[pos]=='}');pos+=1
    return source[start:pos]


def build(root, fallback, repo, output, original=False):
    h=(root/'mg/bddc.hh').read_text();cpp=(root/'mg/bddc.hpp').read_text()
    fields=h[h.index('    std::uint64_t packetOriginalBytes='):h.index('    MpiNodeTopology nodeTopology;',h.index('    std::uint64_t packetOriginalBytes='))]
    pre=base.PREAMBLE
    if not original:
        # Independently capture precisely the bytes passed to MPI_Alltoallv,
        # including all framing/packing. This observes both observer-OFF/ON
        # fixtures and does not alter buffers, counts or the real MPI call.
        capture=r'''
#include <mpi.h>
#include <array>
#include <cstdint>
#include <string>
#include <vector>
std::array<std::uint64_t,1028> capturedCounts{};
std::array<std::uint64_t,2> capturedWireHash{};
std::vector<int> capturedLeaders;
unsigned capturedPhase=0;
std::string entropyOutput;
int captureAlltoallv(void const* send,int const* counts,int const* offsets,MPI_Datatype type,
                    void* recv,int const* receiveCounts,int const* receiveOffsets,
                    MPI_Datatype receiveType,MPI_Comm communicator) {
 int rank=0,ranks=1;MPI_Comm_rank(communicator,&rank);MPI_Comm_size(communicator,&ranks);
 for(int peer=0;peer<ranks;++peer)if(peer!=rank && counts[peer]) {
  auto route=capturedLeaders.at(peer)==capturedLeaders.at(rank)?0u:1u;
  auto base=(capturedPhase*2+route)*257;
  auto bytes=static_cast<unsigned char const*>(send)+offsets[peer];
  auto& hash=capturedWireHash[capturedPhase];
  hash^=peer;hash*=1099511628211ULL;hash^=counts[peer];hash*=1099511628211ULL;
  for(int i=0;i<counts[peer];++i){
   ++capturedCounts[base+bytes[i]];hash^=bytes[i];hash*=1099511628211ULL;
  }
  ++capturedCounts[base+256];
 }
 return MPI_Alltoallv(send,counts,offsets,type,recv,receiveCounts,receiveOffsets,receiveType,communicator);
}
#define MPI_Alltoallv captureAlltoallv
'''
        pre=capture+pre.replace('#include "mg/bddcMpiRouting.hh"',
            '#include "mg/bddcMpiRouting.hh"\n#undef MPI_Alltoallv')
    pre=pre.replace('#include <cstring>', '#include <cstring>\n#include <bit>\n#include <memory>\n#include <mutex>\n#include <iostream>\n#include "mg/bddcPacketPacking.hh"')
    pre=pre.replace('  int id;', '''  struct XVector:std::vector<double>{
    using std::vector<double>::vector;using std::vector<double>::operator=;
    using block_type=double;
  };
  int id;
  std::map<int,XVector> weights;
  void getInterfaceAveragingWeights(int peer,XVector& v) const {v={double(id),double(peer),.125};}
  void setInterfaceAveragingWeights(int peer,XVector const& v) {
    if(!weights.emplace(peer,v).second)throw std::runtime_error("Duplicate averaging weights");
  }''')
    # Keep R/P cardinalities equal, like scalar interface DOFs in production.
    pre=pre.replace('getInterfaceResidualUpdate(to,out); out.push_back(0.125);', 'getInterfaceResidualUpdate(to,out);')
    storage=r'''
using OwnerMap=std::shared_ptr<std::vector<int> const>;
struct OwnedStorage {
 std::vector<std::unique_ptr<SmallSubdomain>> data;
 OwnedStorage(int n,std::vector<int> const& owners,int rank):data(n) {
  for(int i=0;i<n;++i)if(owners[i]==rank){
   data[i]=std::make_unique<SmallSubdomain>();data[i]->id=i;
   for(int j=0;j<n;++j)if(j!=i && (i+j)%3!=0)data[i]->peers.push_back(j);
  }
 }
 int size()const{return static_cast<int>(data.size());}
 bool has(int i)const{return bool(data.at(i));}
 SmallSubdomain& operator[](int i){if(!has(i))throw std::runtime_error("Access to non-owned slot");return *data[i];}
 SmallSubdomain const& operator[](int i)const{if(!has(i))throw std::runtime_error("Const non-owned slot");return *data[i];}
};
template<class F>void parallelFor(int a,int b,F&& f){for(int i=a;i<b;++i)f(i);}
'''
    anchor='template<class Subdomain,class Storage>\nstruct SharedMemoryDomain {'
    pre=pre.replace(anchor,storage+'\n'+anchor)
    a='  SharedMemoryDomain(Storage& s,bool routed,bool master):subdomains(&s),'
    start=pre.index(a);end=pre.index('  bool mpiEnabled()',start)
    has_owner_map='MpiOwnerMap owners' in constructor(cpp)
    declaration='  SharedMemoryDomain(Storage&,bool,bool'+(',MpiOwnerMap' if has_owner_map else '')+');\n'
    pre=pre[:start]+declaration+pre[end:]
    pre=pre.replace('  int restrictionTraffic=0,prolongationTraffic=0,transfers=0;',
        '''  int restrictionTraffic=0,prolongationTraffic=0,transfers=0,initTraffic=0;
  std::mutex mutex;
  std::vector<std::pair<int,int>> neighbors;
  bool losslessPacking=false;
'''+fields)
    pre=pre.replace('  bool mpiEnabled() const', '''  bool storageHas(int i) const{return subdomains->has(i);}
  int mpiRank() const{return mpiLayout.rank();}
  void mpiInitializeAveragingWeights();
  bool mpiEnabled() const''')
    pre=pre.replace('// PRODUCTION_ROUTING_SETTER',base.routing_setter(h))
    # All setup and exchanges below are compiled verbatim from the selected tree.
    start=cpp.index('    template <class T>\n    std::vector<T> mpiAllgatherVector')
    brace=cpp.index('{',start);pos,depth=brace+1,1
    while depth:
        depth+=(cpp[pos]=='{')-(cpp[pos]=='}');pos+=1
    gather='namespace detail {\n'+cpp[start:pos]+'\n}\n'
    methods=gather+constructor(cpp)+'\n'+base.method(cpp,'mpiInitializeAveragingWeights')+'\n'+base.method(cpp,'mpiRestrictExchange')+'\n'+base.method(cpp,'mpiProlongateExchange')
    configure='D d(storage,true,false'+(',map' if has_owner_map else '')+');'
    cases='for(int kind=0;kind<1;++kind)' if original else 'for(int kind=0;kind<4;++kind)'
    observer='' if original else '''d.partitionTransport.enabled=true;
 d.serializedByteEntropy.configure(!entropyOutput.empty());
 capturedCounts={};capturedLeaders=d.nodeTopology.leaders();
 capturedWireHash={14695981039346656037ULL,14695981039346656037ULL};'''
    main=r'''
}
using namespace Kaskade::BDDC;
void require(bool ok,char const* why){if(!ok)throw std::runtime_error(why);}
void check(int kind,bool packed,int rank,int ranks){
 int n=3*ranks+1;
 std::vector<int> owners(n);
 for(int i=0;i<n;++i)owners[i]=kind==0?i%ranks:kind==1?(i+1)%ranks:kind==2?ranks-1:(n-1-i)%ranks;
 OwnerMap map=std::make_shared<std::vector<int> const>(owners);
 OwnedStorage storage(n,owners,rank);
 using D=SharedMemoryDomain<SmallSubdomain,OwnedStorage>;
 CONFIGURE
 d.activeSubdomains.assign(n,true);
 d.setMpiDestinationRouting(true);d.losslessPacking=packed;
 OBSERVER
 for(int i=0;i<n;++i)if(storage.has(i)){
  require(d.mpiOwner(i)==owners[i],"Storage and domain disagree on owner");
  for(int peer:storage[i].peers){
   auto expected=std::vector<double>{double(peer),double(i),.125};
   require(storage[i].weights.at(peer)==expected,"Producer averaging weights were misrouted");
  }
 }
 EXCHANGES
 std::vector<std::uint64_t> hashes(n);
 std::uint64_t active=0;
 for(int i=0;i<n;++i)if(storage.has(i)){
  ++active;
  auto& sub=storage[i];require(sub.restricted.size()==3*sub.peers.size(),"Missing restriction records");
  require(sub.prolongated.size()==3*sub.peers.size(),"Missing prolongation records");
  std::uint64_t hash=14695981039346656037ULL;
  // Gather records in canonical logical endpoint order, not MPI-owner order.
  for(auto records:{sub.restricted,sub.prolongated}){
   std::sort(records.begin(),records.end());
   for(auto const& [from,v]:records){
    SmallSubdomain producer;producer.id=from;std::vector<double> expected;
    producer.getInterfaceResidualUpdate(i,expected);require(v==expected,"Receiver payload changed");
    hash^=from;hash*=1099511628211ULL;
    for(double x:v){hash^=std::bit_cast<std::uint64_t>(x);hash*=1099511628211ULL;}
   }
  }
  hashes[i]=hash;
 }
 MPI_Allreduce(MPI_IN_PLACE,&active,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
 require(active==std::uint64_t(n),"Subdomain is missing or multiply owned");
 MPI_Allreduce(MPI_IN_PLACE,hashes.data(),n,MPI_UINT64_T,MPI_BOR,MPI_COMM_WORLD);
 std::uint64_t all=0;for(auto x:hashes)all^=x;
 if(rank==0)std::cout<<"SIGNATURE "<<kind<<' '<<packed<<' '<<all<<' '
                   <<d.restrictionTraffic<<' '<<d.prolongationTraffic<<'\n';
 CHECK_COUNTERS
 CHECK_BYTE_ENTROPY
}
int main(int argc,char** argv){
 MPI_Init(&argc,&argv);int rank=0,ranks=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
 ENTROPY_ARGUMENT
 try{
  CASES for(bool packed:{false,true}){ CHECK_RUN }
  if(rank==0)std::cout<<"PASS: production owner-only setup, R/P routing, raw/packed parity, three steps and idle ranks\n";
  MPI_Finalize();return 0;
 }catch(std::exception const& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
'''
    counters='' if original else r'''
 require(d.partitionTransport.phase[0].applications==3 && d.partitionTransport.phase[1].applications==3,"Application counters differ");
 std::uint64_t sent=0;
 for(auto const& phase:d.partitionTransport.phase){
  if(d.nodeTopology.nodes()==1)require(phase.route[2].sentBytes==0,"Local-host MPI counted as inter-node");
  require(phase.route[0].messages==0,"Same-rank transfer counted as MPI message");
  sent+=phase.route[1].sentBytes+phase.route[2].sentBytes;
  require(phase.mpi.counts>=0 && phase.mpi.payload>=0,"Negative collective duration");
 }
 require(sent==d.packetSentBytes,"Scoped transport bytes differ from actual serialized bytes");
'''
    byte_entropy='' if original else r'''
 if(d.serializedByteEntropy.enabled()) {
  require(*d.serializedByteEntropy.localCounts()==capturedCounts,"Observer differs from actual MPI send buffer bytes/messages");
  std::array<std::uint64_t,1028> expected{};
  MPI_Reduce(capturedCounts.data(),expected.data(),1028,MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);
  auto directory=entropyOutput+"/map"+std::to_string(kind)+"-packed"+std::to_string(packed);
  d.serializedByteEntropy.writeGlobal(directory);
  if(rank==0) {
   std::ifstream input(directory+"/serialized-message-byte-counts.csv");
   std::string line;std::getline(input,line);unsigned rows=0;
   while(std::getline(input,line)) {
    std::istringstream fields(line);std::string phase,scope,byte,count;
    std::getline(fields,phase,',');std::getline(fields,scope,',');
    std::getline(fields,byte,',');std::getline(fields,count,',');
    unsigned p=phase=="restriction"?0:1,b=std::stoul(byte);
    auto value=scope=="same_node"?expected[(p*2)*257+b]:scope=="inter_node"?
       expected[(p*2+1)*257+b]:expected[(p*2)*257+b]+expected[(p*2+1)*257+b];
    require(std::stoull(count)==value,"Globally pooled histogram differs from MPI_Alltoallv bytes");++rows;
   }
   require(rows==1536,"Incomplete pooled byte histogram export");
  }
 }
'''
    exchanges='for(int step=0;step<3;++step){d.mpiRestrictExchange();d.mpiProlongateExchange();}' if original else '''
 for(int step=0;step<3;++step){capturedPhase=0;d.mpiRestrictExchange();capturedPhase=1;d.mpiProlongateExchange();}'''
    main=main.replace('CONFIGURE',configure).replace('OBSERVER',observer).replace('CASES',cases).replace('CHECK_COUNTERS',counters)
    main=main.replace('CHECK_BYTE_ENTROPY',byte_entropy).replace('EXCHANGES',exchanges)
    main=main.replace('ENTROPY_ARGUMENT','' if original else 'std::string requestedEntropy=argc==2?argv[1]:"";')
    main=main.replace('CHECK_RUN','check(kind,packed,rank,ranks);' if original else '''
   entropyOutput.clear();check(kind,packed,rank,ranks);
   auto off=capturedWireHash;
   if(!requestedEntropy.empty()) {
    entropyOutput=requestedEntropy;check(kind,packed,rank,ranks);
    require(capturedWireHash==off,"Observer changed actual serialized MPI bytes/order");
    if(rank==0)std::cout<<"PASS: observer OFF/ON exact serialized byte/order hash, map="<<kind<<" packed="<<packed<<'\\n';
   }''')
    output.parent.mkdir(parents=True,exist_ok=True)
    source=output.with_suffix('.cpp');source.write_text(pre+methods+main)
    subprocess.run(shlex.split(os.environ.get('MPICXX','mpicxx'))+[
        '-std=c++20','-O2','-Wall','-Wextra','-Werror',
        '-DKASKADE_HAVE_MPI','-I',str(root),'-I',str(fallback),'-I',str(repo),str(source),'-o',str(output)],check=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--kaskade',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--baseline-root',type=Path)
    p.add_argument('--build-only',action='store_true');p.add_argument('--ranks',type=int,default=4)
    a=p.parse_args();repo=Path(__file__).resolve().parents[1]
    build(a.kaskade,a.kaskade,repo,a.output/'ownership-mpi')
    if a.baseline_root:build(a.baseline_root,a.kaskade,repo,a.output/'original-mpi',original=True)
    if a.build_only:return
    launcher=shlex.split(os.environ.get('MPIEXEC','mpiexec'))
    def run(name):
        result=subprocess.run(launcher+['-n',str(a.ranks),str(a.output/name)],capture_output=True,text=True,check=True,timeout=120)
        (a.output/(name+'.log')).write_text(result.stdout+result.stderr)
        print(result.stdout)
        return re.findall(r'^SIGNATURE (.*)$',result.stdout,re.M)
    current=run('ownership-mpi')
    if a.baseline_root:
        original=run('original-mpi')
        assert current[:2]==original,'Default P2 differs from the saved original implementation'
        print('PASS: P2 payload signatures and logical byte counts match saved original')
    values=[x.split() for x in current]
    assert len({x[2] for x in values})==1,'Ownership/raw-packing changed reconstructed fixture values'
    (a.output/'result.json').write_text(json.dumps(dict(status='PASS',ranks=a.ranks,fixture='small MPI algebra/storage, actual production setup and routing',P2_matches_original=bool(a.baseline_root),signatures=current),indent=2)+'\n')


if __name__=='__main__':main()
