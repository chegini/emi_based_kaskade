// Exercise byte-width boundaries, exact constants/corners, plans and wide fallbacks.
#include "mg/bddcPacketPacking.hh"
#include <iostream>
namespace P=Kaskade::BDDC::PacketPacking;
P::Record record(unsigned h,unsigned t,unsigned corners,bool encoded,bool split)
{
  P::Record r{1,2,{},encoded,split};
  P::put(r.bytes,std::uint32_t(0x80000000)); // exact negative-zero parameter
  P::put(r.bytes,std::uint32_t(0x7fc00042)); // opaque NaN payload parameter
  if (encoded)
  {
    P::put(r.bytes,std::uint32_t(h)); P::put(r.bytes,std::uint32_t(t));
    r.bytes.insert(r.bytes.end(),h,0x91); r.bytes.insert(r.bytes.end(),t,0x23);
    r.bytes.resize((r.bytes.size()+7)/8*8,0);
  }
  else r.bytes.insert(r.bytes.end(),h/8*8,0x34);
  for(unsigned i=0;i<corners;++i) P::put(r.bytes,std::uint64_t(0x7ff8000000000042ULL+i));
  if (split) { P::put(r.bytes,std::uint32_t(0x43535031)); P::put(r.bytes,std::uint32_t(0x80000000)); }
  return r;
}
int main()
{
  try {
    unsigned exercised=0;
    for(bool split:{false,true}) for(bool encoded:{false,true})
    for(unsigned h:{0u,1u,7u,8u,110u,118u,119u,120u,126u,127u,128u,254u,255u,256u,32750u,32767u,65535u,65536u})
    for(unsigned t:{0u,3u,7u}) for(unsigned corners:{0u,2u})
    {
      auto r=record(h,t,corners,encoded,split);
      auto ids=P::endpointPlan({r}); auto counts=P::exactCountPlan({r});
      for(bool planned:{false,true})
      {
        auto wire=P::encodeSmall({r},planned,planned);
        if (wire.empty()) continue;
        ++exercised;
        auto back=P::decode(wire,&ids,&counts);
        P::check(back.size()==1 && back[0].bytes==r.bytes &&
                 back[0].encoded==encoded,"Small-format round trip changed bits");
        wire.pop_back();
        bool rejected=false;
        try { (void)P::decode(wire,&ids,&counts); } catch(std::runtime_error const&) { rejected=true; }
        P::check(rejected,"Truncated record accepted");
      }
    }
    auto small=record(5,3,2,true,true);
    auto wideH=record(256,3,0,true,true);
    auto mix=P::encodeSmall({small,wideH},false,false);
    P::check(!(mix[0]&0x10) && !(mix[0]&0x08),"Mixed widths not widened");
    auto mixed=P::decode(mix);
    P::check(mixed[0].bytes==small.bytes && mixed[1].bytes==wideH.bytes,"Mixed body changed");
    auto ordinary=record(0,0,0,true,false);
    P::check(P::encodeSmall({small,ordinary},false,false).empty(),"Mixed split state accepted");
    P::EndpointPlan sp,rp; P::ExactCountPlan sc,rc;
    P::RecordLengthStats lengths; P::DictionaryStats stats;
    // Changes in exact counts, endpoints and raw/encoded modes refresh plans.
    for (auto r:{small,small,record(5,0,1,true,true),record(16,0,0,false,true),
                 small,record(65536,0,0,true,true),small})
    {
      auto base=P::best({r});
      auto wire=P::encodePlanned({r},base,sp,&sc,nullptr,&lengths,&stats);
      auto back=P::decodePlanned(wire,rp,&rc);
      P::check(back.size()==1 && back[0].bytes==r.bytes,"Plan transition changed packet");
      P::check(stats.beforeSmallBytes>=stats.smallSaved,"Invalid savings counters");
    }
    P::check(stats.smallMessages>0 && stats.smallSaved>0,"New layout not selected");
    auto empty=P::encodePlanned({}, {},sp,&sc,nullptr,&lengths,&stats);
    P::check(P::decodePlanned(empty,rp,&rc).empty(),"Empty message changed");
    auto changed=small; changed.to=3;
    auto wire=P::encodePlanned({changed},P::best({changed}),sp,&sc,nullptr,&lengths,&stats);
    auto back=P::decodePlanned(wire,rp,&rc);
    P::check(back[0].to==3 && back[0].bytes==changed.bytes,"Endpoint change lost");
    // One constant tag disappears but its 32-bit value is unchanged.
    auto splitSmall=record(0,0,0,true,true);
    auto explicitWire=P::encodeSmall({splitSmall},true,true);
    P::check(explicitWire.size()==15,"Expected 1 message byte + 14 fixed bytes");
    auto noSplit=record(0,0,0,true,false);
    P::check(P::encodeSmall({noSplit},true,true).size()==11,"Expected 10 fixed bytes");
    // Reject impossible count plans rather than reading outside the body.
    auto ids=P::endpointPlan({splitSmall}); P::ExactCountPlan bad{100000};
    bool rejected=false;
    try { (void)P::decode(explicitWire,&ids,&bad); } catch(std::runtime_error const&) { rejected=true; }
    P::check(rejected,"Invalid exact-count plan accepted");
    std::cout<<"PASS: "<<exercised<<" small-format round trips; boundary widths, raw fallbacks, "
             <<"constant/corner bits, truncation, plan changes and 10/14-byte fixed fields.\n";
  } catch(std::exception const& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
