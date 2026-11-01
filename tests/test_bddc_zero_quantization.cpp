#include "mg/bddcZeroQuantization.hh"
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>

using Q = Kaskade::BDDC::ZeroPreservingQuantization;

int main()
{
  std::mt19937_64 generator(731);
  std::uniform_real_distribution<double> random(-1,1);
  for (int bits : {2,4,8,12,16,24,32})
  {
    double const maxCode = static_cast<double>((std::uint64_t(1)<<bits)-1);
    double const zero = (maxCode+1)/2;
    for (double range : {0.0,1e-20,0.000353188,0.31,1.0,1e20})
    {
      double const scale = Q::scaleForRange(range,maxCode);
      // Exactly the float32 value sent on the wire, without loss on serialization.
      assert(scale==static_cast<double>(static_cast<float>(scale)));
      assert(Q::encode(0,maxCode,scale)==zero);
      assert(Q::decode(zero,maxCode,scale)==0);
      assert(Q::decode(Q::encode(-0.0,maxCode,scale),maxCode,scale)==0);
      for (int n=0;n<1000;++n)
      {
        double const x = range*(n==0 ? 1 : n==1 ? -1 : random(generator));
        double const code=Q::encode(x,maxCode,scale);
        assert(code>=1 && code<=maxCode && std::floor(code)==code);
        assert(Q::encode(-x,maxCode,scale)+code==2*zero);
        double const restored = Q::decode(code,maxCode,scale);
        assert(std::abs(restored-x)<=0.5/scale+16*std::numeric_limits<double>::epsilon()*range);
      }
    }
    // Symmetric tie handling; a large offset must not change rounding direction.
    assert(Q::encode(0.5,maxCode,1)==zero+1);
    assert(Q::encode(-0.5,maxCode,1)==zero-1);
  }
  // The feedback state is refined locally using the same grid as the receiver.
  std::array<double,5> x{0,0.000353188, -0.000074, 1e-12,0.000017};
  std::array<double,5> published{};
  double previous=1;
  for (int pass=0;pass<3;++pass)
  {
    double range=0;
    for (size_t i=0;i<x.size();++i) range=std::max(range,std::abs(x[i]-published[i]));
    double const scale=Q::scaleForRange(range,4095);
    for (size_t i=0;i<x.size();++i)
      published[i]+=Q::decode(Q::encode(x[i]-published[i],4095,scale),4095,scale);
    double error=0;
    for (size_t i=0;i<x.size();++i) error=std::max(error,std::abs(x[i]-published[i]));
    assert(published[0]==0 && error<previous);previous=error;
  }
  bool rejected=false;
  try { Q::scaleForRange(1,1); } catch (std::invalid_argument const&) { rejected=true; }
  assert(rejected);
  rejected=false;
  try { Q::scaleForRange(1e300,4095); } catch (std::runtime_error const&) { rejected=true; }
  assert(rejected);
  rejected=false;
  try { Q::encode(std::numeric_limits<double>::infinity(),4095,1); }
  catch (std::runtime_error const&) { rejected=true; }
  assert(rejected);
  std::cout << "Zero-preserving quantizer checks passed (2/4/8/12/16/24/32 bits).\n";
}
