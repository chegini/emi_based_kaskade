#include <chrono>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <atomic>
#include <cstdlib>
#include <new>
#include <sys/resource.h>
static std::atomic<std::uint64_t> allocations{0},allocatedBytes{0};
void* operator new(std::size_t size) {
 if(auto p=std::malloc(size ? size:1)){++allocations;allocatedBytes+=size;return p;}
 throw std::bad_alloc();
}
void operator delete(void* p) noexcept {std::free(p);}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete[](void* p) noexcept {::operator delete(p);}
#ifdef __cpp_sized_deallocation
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
void operator delete[](void* p,std::size_t) noexcept {std::free(p);}
#endif
using Clock=std::chrono::steady_clock;
volatile std::uint64_t sink=0;
T readGraph(char const* file,unsigned depth) {
 std::ifstream in(file);T t;std::string magic,phase;int call,id,n;
 in>>magic>>phase>>call>>id>>n;
 if(magic!="KASKADE_FROZEN_V1" || n<0)throw std::runtime_error("Frozen header");
 t.globalIndex.resize(n);std::iota(t.globalIndex.begin(),t.globalIndex.end(),0);
 t.excludedGlobal.resize(n);
 for(int i=0;i<n;++i){int exact;double value;in>>exact>>value;t.excludedGlobal[i]=exact;}
 int peers;in>>peers;
 for(int p=0;p<peers;++p){int peer,size;in>>peer>>size;auto& idx=t.localIfs[peer];
  for(int i=0;i<size;++i){int j;in>>j;if(!t.excludedGlobal.at(j))idx.push_back(j);}}
 t.interface_graph.resize(n);t.interface_graph_weights.resize(n);
 for(int i=0;i<n;++i){int size;in>>size;t.interface_graph[i].resize(size);
  for(auto& j:t.interface_graph[i])in>>j;t.interface_graph_weights[i].assign(size,1);}
 in>>magic;if(!in || magic!="END")throw std::runtime_error("Frozen graph");
 t.graph_lifting_levels=depth;t.reuse_lifting_plans=true;return t;
}
template<class Fn>void measure(char const* name,int repeat,unsigned depth,int count,Fn fn){
 auto a=allocations.load(),b=allocatedBytes.load();auto start=Clock::now();
 for(int i=0;i<count;++i)fn();
 auto seconds=std::chrono::duration<double>(Clock::now()-start).count();
 auto na=allocations.load()-a,nb=allocatedBytes.load()-b;
 rusage r{};getrusage(RUSAGE_SELF,&r);
 auto peak=static_cast<std::uint64_t>(r.ru_maxrss);
#ifndef __APPLE__
 peak*=1024;
#endif
 std::cout<<name<<','<<repeat<<','<<depth<<','<<count<<','<<std::setprecision(12)<<seconds<<','<<na<<','<<nb<<','<<peak<<'\n';
}
void savePlans(T const& t,std::string const& name) {
 std::ofstream out(name,std::ios::binary);
 for(auto const* plans:{&t.local_integer_plans,&t.local_float_plans})
  for(auto const& [key,plan]:*plans){auto bytes=plan.pack();auto n=bytes.size();
   out.write(reinterpret_cast<char const*>(&n),sizeof(n));
   out.write(reinterpret_cast<char const*>(bytes.data()),bytes.size());}
}
int main(int argc,char** argv){
 if(argc!=4)return 2;int repeats=std::stoi(argv[2]),copies=std::stoi(argv[3]);
 std::cout<<"operation,repeat,depth,calls,seconds,allocations,allocated_bytes,process_peak_rss_bytes\n";
 for(unsigned depth:{1u,4u})for(int repeat=0;repeat<repeats;++repeat){
  auto t=readGraph(argv[1],depth);
  measure("prepare",repeat,depth,1,[&]{t.prepareIntegerLiftingPlans();});
  if(repeat==0)savePlans(t,"plans-d"+std::to_string(depth)+".bin");
  auto const& idx=t.localIfs.begin()->second;auto const& expected=t.local_integer_plans.at(idx);
  check(t.integerLiftingPlan(idx).pack()==expected.pack());
  measure("accessor",repeat,depth,copies,[&]{auto p=t.integerLiftingPlan(idx);sink=p.fine.size();});
  measure("reference_control",repeat,depth,copies,[&]{auto const& p=t.local_integer_plans.at(idx);sink=p.fine.size();});
  std::vector<std::uint32_t> u(idx.size(),0xffffffffu);
  measure("forward_copy",repeat,depth,100,[&]{auto p=t.integerLiftingPlan(idx);auto v=p.forward<1>(u,true);sink=v.size();});
  measure("forward_reference",repeat,depth,100,[&]{auto const& p=t.local_integer_plans.at(idx);auto v=p.forward<1>(u,true);sink=v.size();});
  measure("production_forward",repeat,depth,100,[&]{auto v=t.applyIntegerGraphLiftingTransform(u,idx,true);sink=v.size();});
  auto lifted=t.applyIntegerGraphLiftingTransform(u,idx,true);
  check(lifted==expected.forward<1>(u,true));
  check(t.inverseIntegerGraphLiftingTransform(lifted,idx,true)==u);
  check(expected.inverse<1>(expected.forward<1>(u,true),true)==u);
 }
}
