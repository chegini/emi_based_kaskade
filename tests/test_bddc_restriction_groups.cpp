#include "mg/bddcRestrictionGroups.hh"
#include "mg/bddcZeroQuantization.hh"
#include <cassert>
#include <iostream>
#include <functional>

int main()
{
  using Kaskade::BDDC::RestrictionGroups;
  using Q=Kaskade::BDDC::ZeroPreservingQuantization;
  auto fails=[](std::function<void()> f) {
    bool rejected=false;
    try { f(); } catch (std::runtime_error const&) { rejected=true; }
    assert(rejected);
  };
  // One producer group appears twice, with a different group in between.
  auto header=RestrictionGroups::pack({0,0,1,0},{2,2,4,2},8);
  assert(header.size()==40); // 4+3*12, already aligned.
  auto [scales,bytes]=RestrictionGroups::unpack(header.data(),header.size(),4,8);
  assert(bytes==40 && scales==std::vector<double>({2,2,4,2}));
  auto empty=RestrictionGroups::pack({}, {}, 8);
  assert(empty.size()==8);
  assert(RestrictionGroups::unpack(empty.data(),8,0,8).first.empty());
  fails([&]{RestrictionGroups::unpack(header.data(),12,4,8);});
  fails([&]{RestrictionGroups::unpack(header.data(),header.size(),3,8);});
  fails([&]{RestrictionGroups::pack({0,0},{2,4},8);});
  fails([&]{RestrictionGroups::pack({0},{0},8);});
  auto corrupt=header;
  float offset=1;
  std::memcpy(corrupt.data()+12,&offset,4);
  fails([&]{RestrictionGroups::unpack(corrupt.data(),corrupt.size(),4,8);});
  // Producer/receiver use the identical float32 grid, including exact zero.
  double s=Q::scaleForRange(0.1234567,4095);
  auto h=RestrictionGroups::pack({1,1},{s,s},8);
  auto decoded=RestrictionGroups::unpack(h.data(),h.size(),2,8).first;
  for (double x : {0.0,-1e-8,0.1234567})
    assert(Q::decode(Q::encode(x,4095,s),4095,s)==
           Q::decode(Q::encode(x,4095,s),4095,decoded[0]));
  std::cout << "Restriction group framing, validation and wire-scale tests passed.\n";
}
