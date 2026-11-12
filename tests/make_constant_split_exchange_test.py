#!/usr/bin/env python3
"""Extract production MPI exchange methods for byte/order regression testing."""
from pathlib import Path
import argparse
p=argparse.ArgumentParser()
p.add_argument('--kaskade-root',type=Path,required=True)
p.add_argument('--output-source',type=Path,required=True)
a=p.parse_args()
k=a.kaskade_root
s=(k/'mg/bddc.hpp').read_text()
parts=[]
for phase in ['Restrict','Prolongate']:
 start=s.index('  template <class Subdomain, class Storage>\n  void SharedMemoryDomain<Subdomain,Storage>::mpi'+phase+'Exchange()')
 end=s.index('\n  }',start)+4;parts.append(s[start:end])
head=r'''
#include "mg/bddcMpiRouting.hh"
#include "mg/bddcPacketPacking.hh"
#include <cstring>
#include <type_traits>
#include <utility>
using namespace Kaskade::BDDC;
struct Mock
{
 using TransmissionType=std::vector<double>;
 int id=0,n=0;
 std::vector<std::vector<unsigned char>> received;
 Mock& transfer(){return *this;}
 void setLosslessPacketInspection(bool){}
 bool lastPacketEncoded(int,bool)const{return true;}
 bool usesProlongateConstantSplit()const{return true;}
 std::vector<int> neighbors()const{std::vector<int> v;for(int j=0;j<n;++j)if(j!=id)v.push_back(j);return v;}
 static TransmissionType packet(int from,int to,bool residual)
 {
  PacketPacking::Bytes b;PacketPacking::put(b,1.0f);PacketPacking::put(b,-2.0f);
  PacketPacking::put(b,std::uint32_t(3));PacketPacking::put(b,std::uint32_t(0));
  b.push_back(static_cast<unsigned char>(from));b.push_back(static_cast<unsigned char>(to));b.push_back(residual);
  b.resize(24,0);PacketPacking::put(b,3.14159);
  if(!residual) { PacketPacking::Bytes prefix;PacketPacking::put(prefix,std::uint32_t(0x43535031));PacketPacking::put(prefix,-17.25f);prefix.insert(prefix.end(),b.begin(),b.end());b=prefix; }
  TransmissionType v(b.size()/8);std::memcpy(v.data(),b.data(),b.size());return v;
 }
 void getInterfaceResidualUpdate(int to,TransmissionType& v){v=packet(id,to,true);}
 void getInterfaceValues(int to,TransmissionType& v){v=packet(id,to,false);}
 void receive(int from,TransmissionType const& v,bool residual)
 {
  auto expected=packet(from,id,residual);
  if(v.size()!=expected.size() || std::memcmp(v.data(),expected.data(),v.size()*8))throw std::runtime_error("MPI record changed");
  auto bytes=reinterpret_cast<unsigned char const*>(v.data());received.emplace_back(bytes,bytes+v.size()*8);
 }
 void setInterfaceResidualUpdate(int from,TransmissionType const& v){receive(from,v,true);}
 void setInterfaceValues(int from,TransmissionType const& v){receive(from,v,false);}
};
struct OwnedStorage
{
 std::vector<Mock>& values;
 std::size_t size()const{return values.size();}
 Mock& operator[](int id) {
  int rank,n;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&n);
  if(id%n!=rank)throw std::runtime_error("access to non-owned BDDC subdomain");
  return values.at(id);
 }
};
namespace Kaskade::BDDC
{
 template<class Subdomain,class Storage> class SharedMemoryDomain
 {
 public:
  using Subdom=Subdomain;Storage* subdomains;
  std::vector<bool> activeSubdomains;
  bool destinationRouting=true,losslessPacking=false;
  MpiNodeTopology nodeTopology;MpiRoutingTraffic routedTraffic;
  int restrictionTraffic=0,prolongationTraffic=0,transfers=0;
  std::uint64_t packetOriginalBytes=0,packetSentBytes=0,packetWideBytes=0,packetTwoLengthBytes=0;
  std::uint64_t packetExplicitBytes=0,packetPlanIdBytes=0,packetOmittedIdBytes=0,packetPlanMessages=0,packetPlanReuseMessages=0,packetWithCountBytes=0;
  std::array<std::vector<PacketPacking::EndpointPlan>,2> packetSentPlans,packetReceivedPlans;
  std::array<std::vector<PacketPacking::ExactCountPlan>,2> packetSentCounts,packetReceivedCounts;
  PacketPacking::RecordLengthStats packetLengthStats;
  PacketPacking::DictionaryStats packetDictionaryStats;
  bool packetTimingEnabled=true;
  std::array<double,5> packetWallSeconds{};
  explicit SharedMemoryDomain(Storage& s):subdomains(&s),activeSubdomains(s.size(),true){nodeTopology.configure();}
  int mpiRank()const{int x;MPI_Comm_rank(MPI_COMM_WORLD,&x);return x;}
  int mpiSize()const{int x;MPI_Comm_size(MPI_COMM_WORLD,&x);return x;}
  int size()const{return static_cast<int>(subdomains->size());}
  int mpiOwner(int id)const{return id%mpiSize();}
  bool mpiOwns(int id)const{return mpiOwner(id)==mpiRank();}
  bool activeOnThisRank(int id)const{return mpiOwns(id);}
  void mpiRestrictExchange();void mpiProlongateExchange();
 };
'''
foot=r'''
}
int main(int argc,char** argv)
{
 MPI_Init(&argc,&argv);int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 try
 {
  std::vector<Mock> subs(8);for(int i=0;i<8;++i){subs[i].id=i;subs[i].n=8;}
  OwnedStorage owned{subs};
  SharedMemoryDomain<Mock,OwnedStorage> d(owned);
  d.mpiRestrictExchange();d.mpiProlongateExchange();
  auto before=subs;
  for(auto& sub:subs)sub.received.clear();
  d.losslessPacking=true;d.packetOriginalBytes=d.packetSentBytes=0;
  for(int repeat=0;repeat<4;++repeat) {
  for(auto& sub:subs)sub.received.clear();
  d.mpiRestrictExchange();d.mpiProlongateExchange();
  for(std::size_t i=0;i<subs.size();++i)
   if(subs[i].received!=before[i].received)throw std::runtime_error("Receive order changed");
  }
  if(!d.packetPlanMessages || !d.packetPlanReuseMessages || d.packetExplicitBytes-d.packetWithCountBytes!=d.packetOmittedIdBytes || d.packetWithCountBytes<=d.packetLengthStats.wideBytes || d.packetLengthStats.wideBytes-d.packetDictionaryStats.wideBytes!=2*d.packetLengthStats.shortRecords || d.packetDictionaryStats.wideBytes-d.packetDictionaryStats.referenceBytes!=3*d.packetDictionaryStats.byteReferences || d.packetDictionaryStats.referenceBytes-d.packetDictionaryStats.beforeSmallBytes!=3*d.packetDictionaryStats.shortCounts || d.packetDictionaryStats.beforeSmallBytes-d.packetSentBytes!=d.packetDictionaryStats.smallSaved)throw std::runtime_error("Endpoint plans not exercised");
  if(d.packetWallSeconds[0]<=0 || d.packetWallSeconds[1]<=0 || d.packetWallSeconds[2]<=0 || d.packetWallSeconds[3]<=0)throw std::runtime_error("Missing timing samples");
  if(d.packetSentBytes>d.packetTwoLengthBytes || d.packetTwoLengthBytes>d.packetWideBytes || d.packetWideBytes>d.packetOriginalBytes)throw std::runtime_error("Expected packing saving");
  if(rank==0)std::cout<<"PASS: extracted production R/P destination exchanges preserve bytes/order with packing on/off.\n";
  MPI_Finalize();return 0;
 }catch(std::exception const& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
'''
header=(k/'mg/bddc.hh').read_text()
start=header.index('    void printPacketTransport() const')
end=header.index('\n    MpiRoutingTraffic const& routingTraffic()',start)
head=head.replace('  void mpiRestrictExchange();void mpiProlongateExchange();',
    '  bool mpiEnabled()const{return true;}\n  bool mpiRankZero()const{return mpiRank()==0;}\n'+
    header[start:end]+'\n  void mpiRestrictExchange();void mpiProlongateExchange();')
foot=foot.replace('  MPI_Finalize();return 0;', '  d.printPacketTransport();\n  MPI_Finalize();return 0;')
a.output_source.write_text(head+'\n'.join(parts)+foot)
