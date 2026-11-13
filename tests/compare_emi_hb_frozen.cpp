// Audit copy of Poisson frozen replay: same production transfer/codec, selectable HB depth.
// Replay captured scalar BDDC inputs through the production transfer and codec.
// This program runs OUTSIDE the live solver: its shared codebooks are independent.
#include "mg/bddcSpaceTransfer.hpp"
#include <filesystem>
#include <tuple>
#include <fstream>
#include <memory>
#include <numeric>

namespace fs = std::filesystem;
unsigned auditLevels=1;
using T = Kaskade::BDDC::SpaceTransferDataCompression<1,double,double,std::uint32_t,std::uint8_t>;
struct Row
{
  std::vector<int> cols;
  struct It
  {
    std::vector<int>::const_iterator p;
    int index() const { return *p; }
    It& operator++() { ++p; return *this; }
    bool operator!=(It b) const { return p != b.p; }
  };
  It begin() const { return {cols.begin()}; }
  It end() const { return {cols.end()}; }
};
struct Input
{
  std::string phase;
  int call=0, id=0;
  std::vector<double> values;
  std::vector<int> exact;
  Kaskade::BDDC::Interfaces interfaces;
  std::vector<Row> graph;
};
void require(bool ok, std::string const& message)
{
  if (!ok) throw std::runtime_error(message);
}
int count(std::istream& in)
{
  int n=-1;
  in >> n;
  require(bool(in) && n>=0 && n<=10000000,"Invalid frozen input count");
  return n;
}
Input read(fs::path const& path)
{
  std::ifstream in(path);
  Input f;
  std::string magic;
  in >> magic >> f.phase >> f.call >> f.id;
  require(magic=="KASKADE_FROZEN_V1" && f.call>0 && f.id>=0 &&
          (f.phase=="restrict" || f.phase=="prolongate"),"Invalid frozen header");
  int n=count(in);
  f.values.resize(n);
  for (int i=0; i<n; ++i)
  {
    int exact=-1;
    in >> exact >> f.values[i];
    require(bool(in) && (exact==0 || exact==1) && std::isfinite(f.values[i]),"Invalid frozen value");
    if (exact) f.exact.push_back(i);
  }
  int peers=count(in);
  std::set<int> seen, covered;
  for (int j=0; j<peers; ++j)
  {
    int peer=count(in), size=count(in);
    require(peer!=f.id && seen.insert(peer).second,"Duplicate/self neighbor");
    std::vector<int> idx(size);
    std::set<int> unique;
    for (int& i : idx)
    {
      i=count(in);
      require(i<n && unique.insert(i).second,"Invalid interface index");
      covered.insert(i);
    }
    f.interfaces.emplace_back(peer,idx);
  }
  require(covered.size()==f.values.size(),"Frozen union is incomplete");
  f.graph.resize(n);
  for (auto& row : f.graph)
  {
    row.cols.resize(count(in));
    for (int& i : row.cols) { i=count(in); require(i<n,"Invalid graph index"); }
  }
  in >> magic;
  require(bool(in) && magic=="END","Incomplete frozen file: "+path.string());
  in >> magic;
  require(in.eof(),"Trailing frozen input");
  return f;
}
void configure(T& t, Input const& f, std::string const& mode, int bits)
{
  bool coded=mode!="quantization", before=mode=="before", after=mode=="after";
  t.setGlobalIndex_active();
  t.setGraphLiftingLevels(auditLevels);
  t.setQuantizationBits(bits);
  t.setExcludedLocalDofs(f.exact);
  t.setRestrictEncoding(coded); t.setProlongateEncoding(coded);
  t.requireRestrictHuffman(coded); t.requireProlongateHuffman(coded);
  t.setRestrictBitlengthEncoding(coded); t.setProlongateBitlengthEncoding(coded);
  t.setRestrictZigzagEncoding(coded); t.setProlongateZigzagEncoding(coded);
  t.enableTransform(Kaskade::BDDC::TransformType::GRAPH_LIFTING);
  t.setRestrictTransform(before); t.setProlongateTransform(before);
  t.setRestrictTransformAfterQuantization(after);
  t.setProlongateTransformAfterQuantization(after);
  t.setCompressionMinInterfaceSize(0);
  t.setGraphLiftingFromMatrix(f.graph);
  T::XVector w(f.values.size()); w=1;
  t.setInterfaceAveragingWeights(f.id,w);
  t.prepareIntegerLiftingPlans();
}
T::Vector vector(std::vector<double> const& values)
{
  T::Vector v(values.size());
  for (std::size_t i=0; i<values.size(); ++i) v[i][0]=values[i];
  return v;
}
void pre(T& t, Input const& f)
{
  auto v=vector(f.values);
  if (f.phase=="restrict") t.preRestrict(v); else t.preProlongate(v);
}
T::TransmissionType packet(T const& t, Input const& f, int peer)
{
  T::TransmissionType wire;
  if (f.phase=="restrict") t.getInterfaceResidualUpdate(peer,wire);
  else t.getInterfaceValues(peer,wire);
  return wire;
}
std::vector<double> own(T& t, Input const& f)
{
  T::Vector zero(f.values.size()), result(f.values.size()); zero=0; result=0;
  if (f.phase=="restrict") t.postRestrict(result);
  else t.postProlongate(zero,result);
  std::vector<double> out(f.values.size());
  for (std::size_t i=0; i<out.size(); ++i) out[i]=result[i][0];
  return out;
}
std::ofstream output(fs::path const& path)
{
  std::ofstream out(path);
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out << std::setprecision(17);
  return out;
}
// Decode only the producer setup plans; no topology is reconstructed here.
auto exportedPlan(T const& t,int peer,std::size_t count) {
 auto bytes=t.exportIntegerLiftingPlan(peer);
 std::uint32_t tag=0;if(bytes.size()>=4)std::memcpy(&tag,bytes.data(),4);
 if(tag==0x464c5001){std::uint64_t n=0;require(bytes.size()>=12,"Short setup");std::memcpy(&n,bytes.data()+4,8);require(n<=bytes.size()-12,"Bad setup");bytes.erase(bytes.begin(),bytes.begin()+12+n);}
 return Kaskade::BDDC::IntegerLiftingPlan<double>::unpack(bytes,count);
}
void replay(fs::path const& root, fs::path const& dest, int bits)
{
  std::map<std::pair<std::string,int>,std::map<int,Input>> groups;
  for (auto const& entry : fs::directory_iterator(root))
    if (entry.path().extension()==".buffer")
    {
      auto f=read(entry.path());
      auto& group=groups[{f.phase,f.call}];
      require(group.emplace(f.id,f).second,"Duplicate producer snapshot");
    }
  require(!groups.empty(),"No frozen inputs");
  auto const& first=groups.begin()->second;
  std::set<int> producers;
  for (auto const& [id,f] : first) producers.insert(id);
  std::map<std::string,int> last;
  for (auto const& [key,group] : groups)
  {
    std::set<int> ids;
    for (auto const& [id,f] : group) ids.insert(id);
    require(ids==producers,"Missing producer in a captured exchange");
    require(key.second==++last[key.first],"Missing captured invocation");
    for (auto const& [id,f] : group)
      for (auto const& [peer,idx] : f.interfaces)
        require(producers.count(peer),"Missing neighbor producer");
  }
  require(last.size()==2,"Capture must contain both R and P");
  require(fs::create_directory(dest),"Output already exists");
  auto packets=output(dest/"packets.csv");
  packets << "mode,phase,call,producer,peer,raw_bytes,wire_bytes,grid_bytes,framing_bytes,"
             "huffman_bytes,tail_bytes,padding_bytes,exact_bytes,slot_bytes,plan_bytes,"
             "squared_error,squared_norm,max_absolute_error\n";
  // Emit bounded per-packet/per-level aggregates, not one row per coefficient.
  auto levels=output(dest/"level-stats.csv");
  levels << "mode,phase,call,producer,peer,level,count,minimum,maximum,sum_squares,tail_bits\n";
  using CacheKey=std::tuple<std::string,int,int>;
  std::map<CacheKey,std::vector<int>> blockCache;
  std::map<std::tuple<std::string,int,unsigned>,std::unique_ptr<T>> probeCache;
  // Plans are reusable only if all captured topologies/corner masks agree.
  std::map<int,Input const*> topology;
  for(auto const& [key,group]:groups) for(auto const& [id,f]:group) {
    auto [it,inserted]=topology.emplace(id,&f);
    if(!inserted){auto const& a=*it->second;
      require(a.exact==f.exact && a.interfaces==f.interfaces && a.graph.size()==f.graph.size(),"Capture topology changed; cannot cache plans");
      for(std::size_t j=0;j<f.graph.size();++j)require(a.graph[j].cols==f.graph[j].cols,"Capture graph changed; cannot cache plans");
    }
  }
  auto books=output(dest/"codebooks.csv");
  books << "mode,phase,call,training_symbols,estimated_table_bits\n";
  for (auto const& [key,group] : groups)
  {
    std::cout<<"Processing "<<key.first<<" call "<<key.second<<" at depth "<<auditLevels<<std::endl;
    bool residual=key.first=="restrict";
    // Quantization-only reconstructed packet values are the exact reference for
    // reversible integer lifting. They need not match transform-first values.
    std::map<std::pair<int,int>,std::vector<double>> nodalReference;
    for (std::string mode : {"quantization","before","after"})
    {
      std::cout<<"  mode "<<mode<<std::endl;
      std::map<int,std::unique_ptr<T>> senders;
      std::vector<std::uint32_t> training;
      for (auto const& [id,f] : group)
      {
        auto t=std::make_unique<T>(id,f.interfaces);
        configure(*t,f,mode,bits);
        pre(*t,f);
        if (residual) t->appendRestrictSymbols(training);
        else t->appendProlongateSymbols(training);
        senders.emplace(id,std::move(t));
      }
      if (mode!="quantization")
      {
        if (residual) T::buildSharedRestrictCodebookFromSymbols(training);
        else T::buildSharedProlongateCodebookFromSymbols(training);
        for (auto& [id,t] : senders)
          if (residual) t->finalizeRestrictEncoding(); else t->finalizeProlongateEncoding();
      }
      books << mode << ',' << key.first << ',' << key.second << ',' << training.size() << ','
            << (mode=="quantization" ? 0 : residual ? T::sharedRestrictCodebookBits()
                                                   : T::sharedProlongateCodebookBits()) << '\n';
      for (auto const& [id,f] : group)
      {
        auto& sender=*senders.at(id);
        auto local=own(sender,f);
        for (auto const& [peer,idx] : f.interfaces)
        {
          auto wire=packet(sender,f,peer);
          auto plan=sender.exportIntegerLiftingPlan(peer);
          Input r;
          r.id=peer; r.phase=f.phase; r.values.resize(idx.size(),0); r.graph.resize(idx.size());
          std::vector<int> sequential(idx.size());
          std::iota(sequential.begin(),sequential.end(),0);
          r.interfaces={{id,sequential}};
          for (std::size_t i=0; i<idx.size(); ++i)
            if (std::binary_search(f.exact.begin(),f.exact.end(),idx[i])) r.exact.push_back(i);
          T receiver(peer,r.interfaces);
          configure(receiver,r,mode,bits);
          receiver.importIntegerLiftingPlan(id,plan);
          pre(receiver,r);
          if (residual) receiver.setInterfaceResidualUpdate(id,wire);
          else receiver.setInterfaceValues(id,wire);
          auto decoded=own(receiver,r);
          double error=0, norm=0, maximum=0;
          for (std::size_t i=0; i<idx.size(); ++i)
          {
            double v=f.values[idx[i]], d=decoded[i]-v;
            require(std::isfinite(decoded[i]),"Nonfinite reconstruction");
            require(decoded[i]==local[idx[i]],"Sender/receiver reconstruction mismatch");
            if (std::binary_search(r.exact.begin(),r.exact.end(),int(i)))
              require(decoded[i]==v,"Exact corner changed");
            error+=d*d; norm+=v*v; maximum=std::max(maximum,std::abs(d));
          }
          if (mode=="quantization") nodalReference[{id,peer}]=decoded;
          if (mode=="after")
            require(decoded==nodalReference.at({id,peer}),"Integer lifting changed nodal reconstruction");
          auto stage=sender.inspectPacketStages(peer,residual);
          auto filtered=idx;
          filtered.erase(std::remove_if(filtered.begin(),filtered.end(),[&](int i){return std::binary_search(f.exact.begin(),f.exact.end(),i);}),filtered.end());
          auto cacheKey=CacheKey{mode,id,peer};
          auto cached=blockCache.find(cacheKey);
          if(cached==blockCache.end()) {
            std::vector<int> block(filtered.size(),0);
            if(mode!="quantization") {
              unsigned saved=auditLevels;
              auto final=exportedPlan(sender,peer,filtered.size());
              for(unsigned cap=1;cap<=final.levels;++cap){
                auto probeKey=std::make_tuple(mode,id,cap);
                auto probe=probeCache.find(probeKey);
                if(probe==probeCache.end()) {
                  auditLevels=cap;
                  auto t=std::make_unique<T>(id,f.interfaces);configure(*t,f,mode,bits);
                  probe=probeCache.emplace(probeKey,std::move(t)).first;
                }
                auto plan=exportedPlan(*probe->second,peer,filtered.size());
                for(int j:plan.fine)if(!block[j])block[j]=cap;
              }
              auditLevels=saved;
              if(mode=="after"){
                auto nodal=block;block.clear();
                for(int j:final.coarse)block.push_back(nodal[j]);
                for(int j:final.fine)block.push_back(nodal[j]);
              }
            }
            cached=blockCache.emplace(cacheKey,std::move(block)).first;
          }
          auto const& block=cached->second;
          require(block.size()==stage.symbols.size(),"Level/symbol alignment mismatch");
          struct Stats {std::size_t count=0;double lo=0,hi=0,squares=0;std::uint64_t tails=0;};
          std::map<int,Stats> stats;
          for(std::size_t j=0;j<block.size();++j) {
            double value=stage.quantizationInput[j];
            if(mode=="after") {
              auto v=std::uint64_t(stage.lifted[j]);
              value=v<(1ULL<<31) ? double(v) : double(std::int64_t(v)-(1LL<<32));
            }
            auto& v=stats[block[j]];
            if(!v.count)v.lo=v.hi=value;
            ++v.count;v.lo=std::min(v.lo,value);v.hi=std::max(v.hi,value);v.squares+=value*value;
            auto symbol=stage.symbols[j];unsigned length=0;for(;symbol;symbol>>=1)++length;
            v.tails+=length ? length-1 : 0;
          }
          for(auto const& [level,v]:stats)
            levels<<mode<<','<<f.phase<<','<<f.call<<','<<id<<','<<peer<<','<<level<<','
                  <<v.count<<','<<v.lo<<','<<v.hi<<','<<v.squares<<','<<v.tails<<'\n';
          std::size_t bytes=wire.size()*sizeof(T::TransmissionType::value_type);
          std::size_t exact=8*r.exact.size(), framing=0, huffman=0, tails=0, padding=0, slots=0;
          require(bytes>=8+exact,"Short packet");
          if (mode=="quantization") slots=8*stage.symbols.size();
          else
          {
            require(bytes>=16+exact,"Short encoded packet");
            std::uint32_t h=0,t=0;
            auto data=reinterpret_cast<unsigned char const*>(wire.data());
            std::memcpy(&h,data+8,4); std::memcpy(&t,data+12,4);
            framing=8; huffman=h; tails=t;
            padding=(8-(huffman+tails)%8)%8;
          }
          require(bytes==8+framing+huffman+tails+padding+exact+slots,"Byte accounting mismatch");
          packets << mode << ',' << f.phase << ',' << f.call << ',' << id << ',' << peer << ','
                  << 8*idx.size() << ',' << bytes << ",8," << framing << ',' << huffman << ','
                  << tails << ',' << padding << ',' << exact << ',' << slots << ',' << plan.size()
                  << ',' << error << ',' << norm << ',' << maximum << '\n';
        }
      }
    }
  }
  std::cout << "Frozen replay passed: production packet roundtrips, exact corners, own/copy equality,"
               " and integer-lifting/nodal equality.\n";
}

// Synthetic fixtures exercise overlap, nonidentical producer graphs, mixed/all
// exact corners, zeros and nonzero data. Capture-on/off packets must be identical.
void fixture(fs::path const& root, int bits)
{
  require(fs::create_directory(root),"Fixture directory exists");
  for (int id=0; id<3; ++id)
  {
    Input f;
    f.id=id; f.values={0.1,-0.4,0,0.7,0.9};
    f.interfaces={{(id+1)%3,{0,1,2,3}},{(id+2)%3,{1,2,3,4}}};
    f.graph.resize(5);
    for (int i=1; i<5; ++i)
    {
      int j=id==1 ? 0 : i-1;
      f.graph[i].cols.push_back(j); f.graph[j].cols.push_back(i);
    }
    if (id==1) f.exact={1};
    if (id==2) f.exact={0,1,2,3,4};
    T observed(id,f.interfaces), control(id,f.interfaces);
    configure(observed,f,"quantization",bits); configure(control,f,"quantization",bits);
    for (int call=0; call<2; ++call)
      for (std::string phase : {"restrict","prolongate"})
      {
        f.phase=phase;
        if (call==1) std::fill(f.values.begin(),f.values.end(),0);
        unsetenv("KASKADE_FROZEN_TRANSFER_DIR");
        pre(control,f);
        setenv("KASKADE_FROZEN_TRANSFER_DIR",root.c_str(),1);
        pre(observed,f);
        unsetenv("KASKADE_FROZEN_TRANSFER_DIR");
        for (auto const& [peer,idx] : f.interfaces)
        {
          auto a=packet(observed,f,peer), b=packet(control,f,peer);
          require(a.size()==b.size() && std::memcmp(a.data(),b.data(),
                  a.size()*sizeof(T::TransmissionType::value_type))==0,"Capture changed packet");
        }
        require(own(observed,f)==own(control,f),"Capture changed own contribution");
      }
  }
}
int main(int argc, char** argv)
{
  try
  {
    require(argc==5,"Usage: replay INPUT_DIR NEW_OUTPUT_DIR BITS LEVELS");
    int bits=std::stoi(argv[3]), levels=std::stoi(argv[4]);
    require(bits>=1 && bits<=32 && levels>=0 && levels<=64,"Invalid bits/depth");
    auditLevels=levels;
    unsetenv("KASKADE_FROZEN_TRANSFER_DIR");
    replay(argv[1],argv[2],bits);
    return 0;
  }
  catch (std::exception const& e)
  {
    std::cerr << "Frozen replay FAILED: " << e.what() << '\n';
    return 1;
  }
}
