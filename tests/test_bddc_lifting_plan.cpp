#include "mg/bddcLiftingPlan.hh"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using Plan=Kaskade::BDDC::IntegerLiftingPlan<double>;
using Symbols=std::vector<std::uint32_t>;
namespace fs=std::filesystem;
void require(bool ok, char const* why)
{
  if (!ok) throw std::runtime_error(why);
}
template<class F> void rejects(F action)
{
  bool caught=false;
  try { action(); } catch (std::runtime_error const&) { caught=true; }
  require(caught,"Invalid plan accepted");
}
template<class T> std::vector<T> read(fs::path const& file)
{
  std::ifstream in(file);
  std::size_t n;
  if (!(in>>n)) throw std::runtime_error("Cannot read "+file.string());
  std::vector<T> data(n);
  for (auto& x: data)
    if (!(in>>x)) throw std::runtime_error("Truncated "+file.string());
  std::string extra;
  require(!(in>>extra),"Trailing capture entries");
  return data;
}
Plan captured(fs::path const& folder, std::string const& prefix)
{
  Plan p;
  p.coarse=read<int>(folder/(prefix+"coarse.txt"));
  p.fine=read<int>(folder/(prefix+"fine.txt"));
  p.nodeCount=p.coarse.size()+p.fine.size();
  auto offsets=read<std::size_t>(folder/(prefix+"offsets.txt"));
  auto neighbors=read<int>(folder/(prefix+"neighbors.txt"));
  auto weights=read<double>(folder/(prefix+"weights.txt"));
  require(offsets.size()==p.fine.size()+1 && offsets.front()==0 &&
          offsets.back()==neighbors.size() && neighbors.size()==weights.size(),"Bad capture offsets");
  for (std::size_t k=0; k<p.fine.size(); ++k)
  {
    require(offsets[k]<=offsets[k+1],"Unsorted capture offsets");
    p.neighbors.emplace_back(neighbors.begin()+offsets[k],neighbors.begin()+offsets[k+1]);
    p.weights.emplace_back(weights.begin()+offsets[k],weights.begin()+offsets[k+1]);
  }
  p.validate();
  return p;
}
int main(int argc, char** argv)
{
  try
  {
    Plan a{4,{0,2},{1,3},{{0,2},{2}},{{1,1},{1}}};
    Plan b{4,{0,3},{1,2},{{0,3},{3}},{{1,1},{1}}};
    a.validate(); b.validate();
    Symbols example{14,75,92,3};
    auto symbols=a.forward<1>(example,true);
    require(b.inverse<1>(symbols,true)!=example,"Test must expose differing endpoint plans");
    auto bytes=a.pack();
    auto received=Plan::unpack(bytes,4);
    require(received.inverse<1>(symbols,true)==example,"Producer plan must restore nodal codes");
    std::mt19937 random(713);
    for (bool zigzag: {false,true})
      for (int trial=0; trial<500; ++trial)
      {
        Symbols x(8);
        for (auto& v: x) v=random();
        if (!trial) x={0,0xffffffffu,0x80000000u,0x7fffffffu,1,2,3,4};
        require(received.inverse<2>(a.forward<2>(x,zigzag),zigzag)==x,
                "Multi-component modular integer round trip failed");
      }
    for (std::size_t n=0; n<bytes.size(); ++n)
      rejects([&]{ Plan::unpack({bytes.begin(),bytes.begin()+n},4); });
    rejects([&]{ Plan::unpack(bytes,5); });
    auto corrupt=bytes; corrupt[0]^=1;
    rejects([&]{ Plan::unpack(corrupt,4); });
    corrupt=bytes; corrupt.push_back(0);
    rejects([&]{ Plan::unpack(corrupt,4); });
    auto invalid=a; invalid.fine[0]=0;
    rejects([&]{ invalid.pack(); });
    invalid=a; invalid.neighbors[0][0]=1;
    rejects([&]{ invalid.pack(); });
    invalid=a; invalid.weights[0][0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{ invalid.pack(); });
    require(Plan::unpack(Plan{}.pack(),0).inverse<1>(Symbols{},true).empty(),"Empty plan");
    std::cout<<"Lifting-plan tests passed: endpoint mismatch, setup serialization, modular codes, malformed input.\n";
    if (argc==2)
    {
      fs::path folder=argv[1];
      int packets=0, oldMismatches=0, changedPlans=0;
      for (auto const& entry: fs::directory_iterator(folder))
      {
        std::string name=entry.path().filename().string();
        auto at=name.find("_send_");
        std::string suffix="_symbols.txt";
        if (at==std::string::npos || !name.ends_with(suffix)) continue;
        auto phase=name.substr(0,at);
        auto pair=name.substr(at+6,name.size()-(at+6)-suffix.size());
        auto sep=pair.find('_');
        std::string sender=phase+"_send_"+pair+"_";
        std::string receiver=phase+"_recv_"+pair.substr(sep+1)+"_"+pair.substr(0,sep)+"_";
        auto source=read<std::uint32_t>(folder/(sender+"nodal.txt"));
        auto encoded=read<std::uint32_t>(folder/(sender+"symbols.txt"));
        auto oldDecoded=read<std::uint32_t>(folder/(receiver+"nodal.txt"));
        auto sendPlan=captured(folder,sender), recvPlan=captured(folder,receiver);
        require(sendPlan.forward<1>(source,true)==encoded,"Forward transform changed saved symbols");
        require(recvPlan.inverse<1>(encoded,true)==oldDecoded,"Cannot reproduce old receiver inverse");
        require(sendPlan.inverse<1>(encoded,true)==source,"Sender inverse failed");
        auto imported=Plan::unpack(sendPlan.pack(),recvPlan.nodeCount);
        require(imported.inverse<1>(encoded,true)==source,"Imported producer inverse failed");
        ++packets;
        oldMismatches+=oldDecoded!=source;
        changedPlans+=sendPlan.pack()!=recvPlan.pack();
      }
      require(packets>0,"No captured packets tested");
      std::cout<<"Captured packets="<<packets<<"; unequal endpoint plans="<<changedPlans
               <<"; reproduced old nodal mismatches="<<oldMismatches
               <<"; mismatches with imported producer plan=0\n";
    }
  }
  catch (std::exception const& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
