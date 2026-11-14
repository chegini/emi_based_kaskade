// Full production transfer/codec check. Compile inside the Dune/MPI container.
#include "mg/bddcSpaceTransfer.hpp"
#include "mg/bddcPacketPacking.hh"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

using T=Kaskade::BDDC::SpaceTransferDataCompression<1,double,double,std::uint32_t,std::uint8_t>;
namespace P=Kaskade::BDDC::PacketPacking;
struct Row {
 std::vector<int> cols;
 struct It {
  std::vector<int>::const_iterator p;
  int index()const{return *p;}
  It& operator++(){++p;return *this;}
  bool operator!=(It b)const{return p!=b.p;}
 };
 It begin()const{return {cols.begin()};}
 It end()const{return {cols.end()};}
};
void require(bool ok,char const* why){if(!ok)throw std::runtime_error(why);}
int main(int argc,char** argv)
{
 try {
  if(argc!=7)throw std::runtime_error("bits depth split recenter input-kind output-dir");
  int bits=std::stoi(argv[1]),depth=std::stoi(argv[2]),kind=std::stoi(argv[5]);
  bool split=std::stoi(argv[3]),recenter=std::stoi(argv[4]);
  std::filesystem::path output=argv[6];std::filesystem::create_directories(output);
  constexpr int n=67;std::vector<int> ids(n);std::iota(ids.begin(),ids.end(),0);
  T a(0,{{1,ids},{2,ids}}),b(1,{{0,ids},{2,ids}}),c(2,{{0,ids},{1,ids}});
  std::vector<T*> transfers{&a,&b,&c};
  std::vector<T::Vector> u(3,T::Vector(n)),du(3,T::Vector(n)),r(3,T::Vector(n));
  std::vector<T::XVector> weights(3,T::XVector(n));
  std::vector<int> exact=kind==2 ? ids : kind==1 ? std::vector<int>{0,31,66} : std::vector<int>{};
  for(int s=0;s<3;++s){
   auto& t=*transfers[s];t.setGlobalIndex_active();t.setQuantizationBits(bits);
   t.setProlongateConstantSplit(split);t.setTransformFirstRecentering(recenter);
   t.setGraphLiftingLevels(depth);t.setReuseLiftingPlans(true);t.setExcludedLocalDofs(exact);
   t.setLosslessPacketInspection(true);
   t.enableTransform(Kaskade::BDDC::TransformType::GRAPH_LIFTING);
   t.setRestrictTransform(true);t.setProlongateTransform(true);
   t.setRestrictTransformAfterQuantization(false);t.setProlongateTransformAfterQuantization(false);
   t.setRestrictEncoding(true);t.setProlongateEncoding(true);
   t.setRestrictBitlengthEncoding(true);t.setProlongateBitlengthEncoding(true);
   t.setRestrictZigzagEncoding(true);t.setProlongateZigzagEncoding(true);
   if(kind==3)t.setCompressionMinInterfaceSize(n+1); // Exercise raw fallback with recenter configured.
   std::vector<Row> graph(n);
   for(int i=1;i<n;++i){int j=i-1;graph[i].cols.push_back(j);graph[j].cols.push_back(i);}
   if(s==1)for(int i=4;i<n;i+=7){graph[i].cols.push_back(i-3);graph[i-3].cols.push_back(i);}
   t.setGraphLiftingFromMatrix(graph);
   for(int i=0;i<n;++i){
    weights[s][i][0]=(s+1)*(1.0+(i*(s+3)%17));
    u[s][i][0]=(kind==4 ? 0 : 212772.48203188565)+(kind==4 ? 0 : (s+1)*std::sin(i*.7));
   }
  }
  for(int s=0;s<3;++s)transfers[s]->setInterfaceAveragingWeights(s,weights[s]);
  for(int s=0;s<3;++s)for(int peer=0;peer<3;++peer)if(s!=peer){
   T::XVector w;transfers[peer]->getInterfaceAveragingWeights(s,w);
   transfers[s]->setInterfaceAveragingWeights(peer,w);
  }
  for(auto* t:transfers)t->prepareIntegerLiftingPlans();
  for(int s=0;s<3;++s)for(int peer=0;peer<3;++peer)if(s!=peer)
   transfers[peer]->importIntegerLiftingPlan(s,transfers[s]->exportIntegerLiftingPlan(peer));
  std::ofstream state(output/"state.bin",std::ios::binary);
  std::ofstream packets(output/"packets.bin",std::ios::binary);
  for(bool residual:{true,false}){
   for(int s=0;s<3;++s){if(residual)transfers[s]->preRestrict(u[s]);else transfers[s]->preProlongate(u[s]);}
   std::vector<std::uint32_t> training;
   for(auto* t:transfers){if(residual)t->appendRestrictSymbols(training);else t->appendProlongateSymbols(training);}
   if(residual)T::buildSharedRestrictCodebookFromSymbols(training);else T::buildSharedProlongateCodebookFromSymbols(training);
   for(auto* t:transfers){if(residual)t->finalizeRestrictEncoding();else t->finalizeProlongateEncoding();}
   for(int s=0;s<3;++s)for(int peer=0;peer<3;++peer)if(s!=peer){
    T::TransmissionType wire;
    if(residual)transfers[s]->getInterfaceResidualUpdate(peer,wire);else transfers[s]->getInterfaceValues(peer,wire);
    auto original=wire;
    bool const encoded=transfers[s]->lastPacketEncoded(peer,residual);
    if(kind==3)require(!encoded,"Size threshold did not force raw fallback");
    if(split && !residual){
     std::uint32_t tag;std::memcpy(&tag,wire.data(),4);require(tag==0x43535031u,"CSP1 constant prefix changed");
    }
    auto bptr=reinterpret_cast<unsigned char const*>(wire.data());
    P::Record record{s,peer,P::Bytes(bptr,bptr+wire.size()*sizeof(T::TEntry)),transfers[s]->lastPacketEncoded(peer,residual)};
    bool prefix=split && !residual;
    if(prefix)P::moveConstantPrefix(record.bytes,false);
    record.splitPrefix=prefix;
    // Exercise real compact encoding in both stateless and plan-reuse modes.
    auto best=P::best({record});auto back=P::decode(best);
    if(prefix)P::moveConstantPrefix(back[0].bytes,true);
    require(back[0].bytes.size()==wire.size()*sizeof(T::TEntry) &&
            std::memcmp(back[0].bytes.data(),wire.data(),back[0].bytes.size())==0,"Packed packet not exact");
    P::EndpointPlan sp,rp;P::ExactCountPlan sc,rc;
    for(int repeat=0;repeat<3;++repeat){
     auto encoded=P::encodePlanned({record},best,sp,&sc);auto decoded=P::decodePlanned(encoded,rp,&rc);
     if(prefix)P::moveConstantPrefix(decoded[0].bytes,true);
     require(decoded[0].bytes==back[0].bytes,"Planned compact packet not exact");
    }
    if(residual)transfers[peer]->setInterfaceResidualUpdate(s,wire);else transfers[peer]->setInterfaceValues(s,wire);
    require(wire.size()==original.size() && std::memcmp(wire.data(),original.data(),wire.size()*sizeof(T::TEntry))==0,"Receiver modified packet");
    auto size=wire.size();packets.write(reinterpret_cast<char const*>(&size),sizeof(size));
    packets.write(reinterpret_cast<char const*>(wire.data()),wire.size()*sizeof(T::TEntry));
   }
   for(int s=0;s<3;++s){
    if(residual){r[s]=0;transfers[s]->postRestrict(r[s]);}
    else transfers[s]->postProlongate(u[s],du[s]);
    for(int i=0;i<n;++i){double value=residual ? r[s][i][0] : u[s][i][0]+du[s][i][0];
     require(std::isfinite(value),"Nonfinite reconstructed value");state.write(reinterpret_cast<char const*>(&value),sizeof(value));}
   }
  }
  if(recenter && kind!=2){
   a.setProlongateConstantSplit(!split);
   bool rejected=false;
   try{a.preProlongate(u[0]);}catch(std::runtime_error const&){rejected=true;}
   require(rejected,"Changing split did not invalidate recentered producer setup");
  }
  require(bool(state) && bool(packets),"Cannot write native test output");
 }catch(std::exception const& e){std::cerr<<e.what()<<'\n';return 1;}
}
