// Diagnostic-only, read-only message entropy. No production headers are edited.
#ifndef EMI_MESSAGE_ENTROPY_HH
#define EMI_MESSAGE_ENTROPY_HH
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>
namespace Kaskade::BDDC::MessageEntropy {
inline std::string prefix(){auto p=std::getenv("EMI_MESSAGE_ENTROPY_PREFIX");return p?p:"";}
inline int rank(){auto p=std::getenv("OMPI_COMM_WORLD_RANK");if(!p)p=std::getenv("PMI_RANK");if(!p)p=std::getenv("SLURM_PROCID");return p?std::stoi(p):0;}
inline unsigned phaseId(std::string const& p){if(p=="restrict")return 0;if(p=="prolongate")return 1;throw std::runtime_error("Unknown entropy phase");}
class Collector;
inline Collector& collector();
inline std::vector<int> routes{0};
inline std::vector<int> ownerBySubdomain;
inline void topology(int r,int n,std::vector<int> const& p,std::vector<int> const& owners={}){
 if(prefix().empty())return;
 if(n<1 || r!=rank() || p.size()!=static_cast<std::size_t>(n))throw std::runtime_error("Entropy topology mismatch");
 for(int owner:owners)if(owner<0 || owner>=n)throw std::runtime_error("Invalid entropy owner map");
 routes=p;ownerBySubdomain=owners;(void)collector();
}
inline unsigned peerRoute(std::string const& directory,int peer){
 auto name=std::filesystem::path(directory).filename().string();
 if(!name.starts_with("subdomain"))throw std::runtime_error("Missing entropy producer identity");
 int id=std::stoi(name.substr(9));
 auto owner=[](int s){
  if(s<0)throw std::runtime_error("Negative entropy subdomain");
  return ownerBySubdomain.empty()?s%static_cast<int>(routes.size()):ownerBySubdomain.at(s);
 };
 if(owner(id)!=rank())throw std::runtime_error("Entropy producer differs from actual subdomain owner");
 return static_cast<unsigned>(routes.at(owner(peer)));
}
struct Total {std::uint64_t n=0,packets=0;long double weighted=0;};
class Collector {
 std::mutex mutex;
 std::unordered_map<std::uint64_t,std::uint64_t> counts;
 std::map<std::uint64_t,Total> totals;
 std::size_t limit=500000,part=0;
 std::string stem;
 static void word(std::ofstream& out,std::uint64_t x){for(unsigned i=0;i<8;++i)out.put(static_cast<char>((x>>(8*i))&255));}
 void flush(){
  if(counts.empty())return;
  std::vector<std::pair<std::uint64_t,std::uint64_t>> sorted(counts.begin(),counts.end());std::sort(sorted.begin(),sorted.end());
  std::ofstream out(stem+"-part-"+std::to_string(part++)+".bin",std::ios::binary);
  out.exceptions(std::ios::badbit|std::ios::failbit);out.write("EMIENT01",8);
  for(auto const& [key,n]:sorted){word(out,key);word(out,n);}out.close();counts.clear();
 }
public:
 Collector(){
  auto p=prefix();if(p.empty())throw std::runtime_error("Entropy prefix missing");
  std::filesystem::create_directories(std::filesystem::path(p).parent_path());
  stem=p+"-rank-"+std::to_string(rank())+"-pid-"+std::to_string(getpid());
  if(auto l=std::getenv("EMI_ENTROPY_MAX_BINS"))limit=std::stoull(l);
  if(limit<16)throw std::runtime_error("Entropy bin limit must be >=16");
 }
 ~Collector() noexcept {
  try {
   flush();std::ofstream out(stem+"-totals.csv");out.exceptions(std::ios::badbit|std::ios::failbit);
   out<<"phase,stage,flags,bits,symbols,messages,weighted_entropy_bits\n"<<std::setprecision(21);
   for(auto const& [key,t]:totals)out<<((key>>40)&255)<<','<<((key>>32)&255)<<','<<((key>>48)&255)<<','<<((key>>56)&255)<<','<<t.n<<','<<t.packets<<','<<t.weighted<<'\n';
   out.close();
  }catch(std::exception const& e){std::fprintf(stderr,"Entropy flush failed: %s\n",e.what());std::_Exit(86);}
 }
 template<class U> long double observe(unsigned phase,unsigned stage,unsigned flags,unsigned bits,std::vector<U> const& v){
  if(v.empty())return 0;
  if(phase>1 || stage>3 || bits>32 || flags>255)throw std::runtime_error("Invalid entropy bucket");
  std::unordered_map<std::uint32_t,std::uint64_t> local;
  for(auto x:v){if(static_cast<std::uint64_t>(x)>std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("Entropy word exceeds uint32");++local[static_cast<std::uint32_t>(x)];}
  long double weighted=static_cast<long double>(v.size())*std::log2(static_cast<long double>(v.size()));
  for(auto const& [x,n]:local){(void)x;weighted-=static_cast<long double>(n)*std::log2(static_cast<long double>(n));}
  if(weighted<0 && weighted>-1e-9L)weighted=0;
  std::uint64_t bucket=(std::uint64_t(stage)<<32)|(std::uint64_t(phase)<<40)|(std::uint64_t(flags)<<48)|(std::uint64_t(bits)<<56);
  std::lock_guard<std::mutex> lock(mutex);
  auto& t=totals[bucket];t.n+=v.size();++t.packets;t.weighted+=weighted;
  for(auto const& [x,n]:local){counts[bucket|x]+=n;if(counts.size()>=limit)flush();}
  return weighted;
 }
};
inline Collector& collector(){static Collector value;return value;}
template<class U> void coefficients(std::string const& directory,std::string const& phase,int peer,int bits,
 std::vector<U> const& input,std::vector<U> const& mapped,std::vector<U> const& lengths,bool encoded,bool floating,bool integer){
 if(prefix().empty())return;
 if(bits<2 || bits>32 || input.size()!=mapped.size() || mapped.size()!=lengths.size())throw std::runtime_error("Entropy coefficient population mismatch");
 for(std::size_t i=0;i<mapped.size();++i){std::uint32_t len=0;for(auto x=mapped[i];x;x>>=1)++len;if(len!=lengths[i])throw std::runtime_error("Entropy bit lengths differ from encoder");}
 unsigned flags=(floating?1u:0u)|(integer?2u:0u)|(encoded?4u:0u)|(peerRoute(directory,peer)<<3);
 auto& c=collector();auto p=phaseId(phase);
 auto original=c.observe(p,0,flags,bits,input);
 auto recoded=c.observe(p,1,flags,bits,encoded?mapped:input);
 if(encoded)c.observe(p,2,flags,bits,lengths);
 if(floating && !integer && std::abs(original-recoded)>1e-10L*std::max<std::size_t>(1,input.size()))throw std::runtime_error("Scalar recentering changed per-message Shannon entropy");
}
inline void routed(std::string const& phase,int route,std::vector<unsigned char> const& bytes){
 if(prefix().empty())return;
 if(route<0 || route>2)throw std::runtime_error("Invalid entropy route");
 collector().observe(phaseId(phase),3,static_cast<unsigned>(route)<<3,8,bytes);
}
}
#endif
