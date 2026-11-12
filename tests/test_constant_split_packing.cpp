// Transport regression: preserve the split prefix, stream bits and exact values.
#include "mg/bddcPacketPacking.hh"
#include <iostream>
namespace P=Kaskade::BDDC::PacketPacking;
int main() { try {
 for(bool encoded:{false,true}) for(unsigned h=0;h<8;++h) for(unsigned corners:{0u,2u}) {
  P::Bytes original;
  P::put(original,std::uint32_t(0x43535031));
  P::put(original,std::uint32_t(0x80000000)); // negative-zero constant: retain bits
  P::put(original,0.03125f);P::put(original,-2.0f);
  P::put(original,std::uint32_t(h));P::put(original,std::uint32_t(3));
  original.insert(original.end(),h,0x9a);original.insert(original.end(),3,0x13);
  original.resize((original.size()+7)/8*8,0);
  for(unsigned c=0;c<corners;++c)P::put(original,std::uint64_t(0x7ff8000000000042ULL+c));
  P::Record r{0,1,original,encoded};P::moveConstantPrefix(r.bytes,false);
  for(int mode=0;mode<=6;++mode) {
   auto wire=P::encode({r},mode);auto back=P::decode(wire);
   P::moveConstantPrefix(back.at(0).bytes,true);
   P::check(back[0].bytes==original,"Split prefix/stream/corner changed");
  }
  P::EndpointPlan sp,rp;P::ExactCountPlan sc,rc;
  for(int call=0;call<3;++call) {
   auto wire=P::encodePlanned({r},P::best({r}),sp,&sc);
   auto back=P::decodePlanned(wire,rp,&rc);
   P::moveConstantPrefix(back.at(0).bytes,true);
   P::check(back[0].bytes==original,"Planned split packet changed");
  }
 }
 // A raw empty interface can contain just prefix plus quantization parameters.
 P::Bytes empty;P::put(empty,std::uint32_t(0x43535031));P::put(empty,0.0f);
 P::put(empty,1.0f);P::put(empty,0.0f);
 auto originalEmpty=empty;P::moveConstantPrefix(empty,false);P::moveConstantPrefix(empty,true);
 P::check(empty==originalEmpty,"Empty raw interface changed");
 for(unsigned n:{0u,8u,16u,25u}) {
  bool rejected=false;P::Bytes b(n);
  try {P::moveConstantPrefix(b,false);}catch(std::runtime_error const&){rejected=true;}
  P::check(rejected,"Malformed split prefix accepted");
 }
 std::cout<<"PASS: split prefix, all padding widths, raw fallback, compact formats, exact bits and repeated plans.\n";
 }catch(std::exception const& e){std::cerr<<e.what()<<'\n';return 1;} }
