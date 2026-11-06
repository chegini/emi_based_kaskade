#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "mg/bddcSpaceTransfer.hpp"

using T=Kaskade::BDDC::SpaceTransferDataCompression<1,double,double,std::uint32_t,std::uint8_t>;
void require(bool ok, char const* message) {if(!ok) throw std::runtime_error(message);}

std::vector<double> exercise(bool split, bool encoded, int exact, bool constant, bool before=false,
                             bool both=false)
{
  constexpr int n=67;
  std::vector<int> ids(n);std::iota(ids.begin(),ids.end(),0);
  T a(0,{{1,ids},{2,ids}}), b(1,{{0,ids},{2,ids}}), c(2,{{0,ids},{1,ids}});
  std::vector<T*> t{&a,&b,&c};
  std::vector<T::Vector> u(3,T::Vector(n)),du(3,T::Vector(n));
  std::vector<T::XVector> weights(3,T::XVector(n));
  std::vector<int> corners=exact==2 ? ids : exact==1 ? std::vector<int>{0,31,66} : std::vector<int>{};
  for (int s=0;s<3;++s)
  {
    t[s]->setProlongateConstantSplit(split); // Before the existing setup exchange.
    t[s]->setGlobalIndex_active();t[s]->setQuantizationBits(8);
    t[s]->setExcludedLocalDofs(corners);
    t[s]->setProlongateEncoding(encoded);t[s]->requireProlongateHuffman(encoded);
    t[s]->setProlongateBitlengthEncoding(encoded);t[s]->setProlongateZigzagEncoding(encoded);
    t[s]->enableTransform((encoded || before) ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                                 : Kaskade::BDDC::TransformType::NONE);
    t[s]->setProlongateTransform(before);
    t[s]->setTransformBeforeAndAfterQuantization(both);
    t[s]->setProlongateTransformAfterQuantization(encoded && (!before || both));
    for(int i=0;i<n;++i)
    {
      // Strongly varying weights, a large non-float32-exact offset, small variation.
      weights[s][i][0]=(s+1)*(1.0+(i*(s+3)%17));
      u[s][i][0]=212772.48203188565+(constant ? 0 : (s+1)*std::sin(i*.7));
    }
  }
  for(int s=0;s<3;++s)t[s]->setInterfaceAveragingWeights(s,weights[s]);
  for(int s=0;s<3;++s)for(int peer=0;peer<3;++peer)if(s!=peer)
  {
    T::XVector wireWeights;t[peer]->getInterfaceAveragingWeights(s,wireWeights);
    t[s]->setInterfaceAveragingWeights(peer,wireWeights);
  }
  for(auto* x:t) x->prepareIntegerLiftingPlans();
  for(int s=0;s<3;++s)for(int peer=0;peer<3;++peer)if(s!=peer)
    t[peer]->importIntegerLiftingPlan(s,t[s]->exportIntegerLiftingPlan(peer));
  // A transform changes the quantization range and can amplify error on inverse
  // lifting. Bound that error analytically instead of reusing the nodal fixture's
  // empirical .08 limit. The predictor is a convex average, so ||T^-1||_inf <= 2.
  std::vector<double> beforeBound(3,0);
  if(before && split && exact!=2)
    for(int s=0;s<3;++s)
    {
      std::vector<int> idx;
      long double sum=0;
      for(int i=0;i<n;++i)
        if(std::find(corners.begin(),corners.end(),i)==corners.end())
        {idx.push_back(i);sum+=u[s][i][0];}
      float center=static_cast<float>(sum/idx.size());
      std::vector<double> raw;
      for(int i:idx)raw.push_back(weights[s][i][0]*(u[s][i][0]-double(center)));
      auto coefficients=t[s]->applyTransform(raw,idx);
      auto [lo,hi]=std::minmax_element(coefficients.begin(),coefficients.end());
      double scale=255.0/(*hi==*lo ? 1 : *hi-*lo);
      double wireScale=static_cast<float>(scale),wireMin=static_cast<float>(*lo);
      double coefficientError=.5/scale+std::abs(wireMin-*lo)+
                              255*std::abs(1/wireScale-1/scale);
      beforeBound[s]=2*coefficientError;
    }
  for(int s=0;s<3;++s)t[s]->preProlongate(u[s]);
  if(encoded)
  {
    std::vector<std::uint32_t> symbols;
    for(auto* x:t)x->appendProlongateSymbols(symbols);
    T::buildSharedProlongateCodebookFromSymbols(symbols);
    for(auto* x:t)x->finalizeProlongateEncoding();
  }
  std::uint64_t sent=0,received=0;
  for(int s=0;s<3;++s)for(int peer=0;peer<3;++peer)if(s!=peer)
  {
    T::TransmissionType packet;t[s]->getInterfaceValues(peer,packet);
    if(!encoded)
      require(packet.size()==std::size_t(n+1+(split ? 1 : 0)),"Plain packet/header size");
    if(split)
    {
      std::uint32_t tag=0;std::memcpy(&tag,packet.data(),4);
      require(tag==0x43535031,"Missing CSP1 header");
      auto bad=packet;std::uint32_t invalid=0;std::memcpy(bad.data(),&invalid,4);
      bool rejected=false;
      try {t[peer]->setInterfaceValues(s,bad);}catch(std::runtime_error const&) {rejected=true;}
      require(rejected,"Invalid header accepted");
      bad=packet;float nan=std::numeric_limits<float>::quiet_NaN();
      std::memcpy(reinterpret_cast<unsigned char*>(bad.data())+4,&nan,4);
      rejected=false;
      try {t[peer]->setInterfaceValues(s,bad);}catch(std::runtime_error const&) {rejected=true;}
      require(rejected,"Nonfinite constant accepted");
    }
    t[peer]->setInterfaceValues(s,packet);
    sent+=packet.size()*sizeof(T::TEntry);
  }
  for(int s=0;s<3;++s)t[s]->postProlongate(u[s],du[s]);
  std::vector<double> result;
  double maxError=0;
  for(int i=0;i<n;++i)
  {
    double expected=0,denom=0;
    for(int s=0;s<3;++s){expected+=weights[s][i][0]*u[s][i][0];denom+=weights[s][i][0];}
    expected/=denom;
    double actual=u[0][i][0]+du[0][i][0];result.push_back(actual);
    maxError=std::max(maxError,std::abs(actual-expected));
    if(before && split)
    {
      double bound=std::accumulate(beforeBound.begin(),beforeBound.end(),0.0)/denom;
      require(std::abs(actual-expected)<=bound+1e-8,
              "Floating split exceeds quantization/metadata/inverse-lifting error bound");
    }
    for(int s=1;s<3;++s)
      require(std::abs(actual-u[s][i][0]-du[s][i][0])<1e-8,"Shared copies disagree");
    if(std::find(corners.begin(),corners.end(),i)!=corners.end())
      require(std::abs(actual-expected)<1e-8,"Exact trailer corrupted");
  }
  std::uint64_t counted=0;
  for(auto* x:t)
  {
    counted+=x->prolongateStats().total_sent_bytes;
    received+=x->prolongateStats().total_received_bytes;
    require(x->prolongateConstantHeaderBytes()==(split ? 16 : 0),"Extra header counter");
  }
  require(sent==counted && sent==received,"Packet traffic accounting mismatch");
  if(split && !before)require(maxError<(constant ? 1e-4 : .08),"Constant split failed to protect small variation");
  if(!split && !exact)require(maxError>10,"Fixture does not expose weighted-offset damage");
  return result;
}

void checkDifferentLocalNumbering()
{
  std::vector<int> ia{1,4,0},ib{3,0,2};
  T a(0,{{1,ia}}),b(1,{{0,ib}});
  for(auto* t:{&a,&b})
  {t->setProlongateConstantSplit(true);t->setGlobalIndex_active();t->setQuantizationBits(8);}
  a.setExcludedLocalDofs({4});b.setExcludedLocalDofs({0});
  T::XVector wa(5),wb(5);wa=1;wb=1;
  T::Vector ua(5),ub(5),da(5),db(5);ua=0;ub=0;
  for(int i=0;i<3;++i)
  {
    wa[ia[i]][0]=i+1;wb[ib[i]][0]=10-i;
    ua[ia[i]][0]=200000.1+i;ub[ib[i]][0]=200000.1-i;
  }
  a.setInterfaceAveragingWeights(0,wa);b.setInterfaceAveragingWeights(1,wb);
  T::XVector wire;a.getInterfaceAveragingWeights(1,wire);b.setInterfaceAveragingWeights(0,wire);
  b.getInterfaceAveragingWeights(0,wire);a.setInterfaceAveragingWeights(1,wire);
  a.preProlongate(ua);b.preProlongate(ub);
  T::TransmissionType ab,ba;a.getInterfaceValues(1,ab);b.getInterfaceValues(0,ba);
  b.setInterfaceValues(0,ab);a.setInterfaceValues(1,ba);
  a.postProlongate(ua,da);b.postProlongate(ub,db);
  for(int i=0;i<3;++i)
    require(std::abs(ua[ia[i]][0]+da[ia[i]][0]-ub[ib[i]][0]-db[ib[i]][0])<1e-8,
            "Producer weights/constant lost under different local numbering");
}

int main()
{
  checkDifferentLocalNumbering();
  for(int exact: {0,1,2})for(bool constant: {false,true})
  {
    auto plain=exercise(true,false,exact,constant),full=exercise(true,true,exact,constant);
    require(plain==full,"Integer lifting/Huffman changed split reconstruction");
    exercise(false,false,exact,constant);
    auto beforePlain=exercise(true,false,exact,constant,true);
    auto beforeFull=exercise(true,true,exact,constant,true);
    require(beforePlain==beforeFull,"Huffman changed floating split reconstruction");
    auto both=exercise(true,true,exact,constant,true,true);
    require(beforeFull==both,"Combined lifting changed weighted split reconstruction");
  }
  T missing(0,{{1,{0,1}}});missing.setGlobalIndex_active();missing.setQuantizationBits(8);
  missing.setProlongateConstantSplit(true);T::XVector w(2);w=1;
  missing.setInterfaceAveragingWeights(0,w);T::Vector u(2);u=1;missing.preProlongate(u);
  T::TransmissionType packet;missing.getInterfaceValues(1,packet);
  bool rejected=false;
  try {missing.setInterfaceValues(1,packet);}catch(std::runtime_error const&) {rejected=true;}
  require(rejected,"Missing producer weights accepted");
  std::cout<<"Constant-split transfer passed: weighted offsets, float32 center, three owners, "
             "both transform orders, plain/full equivalence, exact trailers, invalid headers, packet counters.\n";
}
