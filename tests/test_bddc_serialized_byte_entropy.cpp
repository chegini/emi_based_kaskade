// Standalone observer mathematics/export; no PDE solve or MPI initialization.
#include "mg/bddcSerializedByteEntropy.hh"

#include <numeric>

using Observer=Kaskade::BDDC::SerializedByteEntropy;
void require(bool ok,char const* why)
{ if (!ok) throw std::runtime_error(why); }

int main(int argc,char** argv)
{
  try
  {
    require(argc==2,"Expected fresh temporary output directory");
    Observer disabled;
    require(!disabled.enabled() && !disabled.localCounts(),"Disabled observer allocated counts");
    disabled.observe(0,1,{0,1});disabled.writeGlobal(std::string(argv[1])+"/disabled");
    Observer empty;empty.configure(true);
    require(!Observer::summarize(*empty.localCounts(),0,2).bitsPerByte,"Empty entropy is defined");
    Observer uniform;uniform.configure(true);
    std::vector<unsigned char> all(256);std::iota(all.begin(),all.end(),0);
    uniform.observe(0,1,all);
    require(Observer::summarize(*uniform.localCounts(),0,2).bitsPerByte==8,"Uniform H !=8");
    Observer constant;constant.configure(true);constant.observe(0,2,{7,7,7,7});
    require(Observer::summarize(*constant.localCounts(),0,2).bitsPerByte==0,"Constant H !=0");
    Observer pooled;pooled.configure(true);
    pooled.observe(0,1,{0,0,0,0,0});pooled.observe(0,2,{255,255,255,255,255});
    require(Observer::summarize(*pooled.localCounts(),0,0).bitsPerByte==0,"Local stream H !=0");
    require(Observer::summarize(*pooled.localCounts(),0,2).bitsPerByte==1,"Pooled H !=1");
    // Treat a whole synthetic serialized record, including its header/padding,
    // as bytes. The observer must not strip metadata or change input data.
    std::vector<unsigned char> record{0,0,0,1,0,0,0,2,16,0,0,0,0,0,0,0,
                                      0,255,127,128,42,42,42,0};
    auto before=record;pooled.observe(1,1,record);
    require(record==before,"Observer changed serialized bytes");
    pooled.observe(1,2,{}); // Empty destination is not an MPI message.
    require(Observer::summarize(*pooled.localCounts(),1,2).messages==1,"Empty buffer counted");
    auto precision=std::cout.precision();pooled.writeGlobal(argv[1]);
    require(std::cout.precision()==precision,"Observer changed log precision");
    bool rejected=false;
    try { pooled.writeGlobal(argv[1]); } catch (std::runtime_error const&) { rejected=true; }
    require(rejected,"Observer overwrote prior measurements");
    std::cout<<"PASS: disabled/empty, H=0/1/8, pooling, all serialized bytes, unchanged buffer, no overwrite\n";
  }
  catch (std::exception const& error) { std::cerr<<error.what()<<'\n';return 1; }
}
