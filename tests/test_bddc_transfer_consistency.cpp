#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "mg/bddcSpaceTransfer.hpp"

using Transfer = Kaskade::BDDC::SpaceTransferDataCompression<1,double,double,
                                                         std::uint32_t,std::uint8_t>;

void checkClose(double a, double b, char const* phase)
{
  double const tolerance = 64 * std::numeric_limits<double>::epsilon()
                           * std::max({1.0,std::abs(a),std::abs(b)});
  if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b) > tolerance)
  {
    std::cerr << phase << ": copies disagree by " << std::abs(a-b) << "\n";
    throw std::runtime_error("Compressed transfer must reconstruct identical interface sums");
  }
}

std::vector<double> checkTransfer(int bits, bool encoded, int exactDofs, bool zeroPreserving,
                                  int groupSize=0, bool zeroProlongate=false, int prolongateGroupSize=0)
{
  std::vector<double> reconstruction;
  Transfer a(0,{{1,{0,1,2}}}), b(1,{{0,{0,1,2}}});
  for (Transfer* t : {&a,&b})
  {
    t->setGlobalIndex_active();
    t->setQuantizationBits(bits);
    t->setRestrictZeroPreserving(zeroPreserving);
    t->setProlongateZeroPreserving(zeroProlongate);
    t->setRestrictGroupSize(groupSize);
    t->setProlongateGroupSize(prolongateGroupSize);
    if (groupSize) t->setProlongateQuantizationBits(16);
    t->setRestrictEncoding(encoded);
    t->setProlongateEncoding(encoded);
    t->requireRestrictHuffman(encoded);
    t->requireProlongateHuffman(encoded);
    t->setRestrictBitlengthEncoding(encoded);
    t->setProlongateBitlengthEncoding(encoded);
    t->setRestrictZigzagEncoding(encoded);
    t->setProlongateZigzagEncoding(encoded);
    t->enableTransform(encoded ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                               : Kaskade::BDDC::TransformType::NONE);
    t->setRestrictTransformAfterQuantization(encoded);
    t->setProlongateTransformAfterQuantization(encoded);
  }
  std::vector<int> exact;
  if (exactDofs == 1) exact = {1};
  if (exactDofs == 3) exact = {0,1,2};
  a.setExcludedLocalDofs(exact); b.setExcludedLocalDofs(exact);
  Transfer::XVector wa(3), wb(3);
  for (int i=0;i<3;++i)
  {
    wa[i][0]=0.2+0.1*i;
    wb[i][0]=1-wa[i][0];
  }
  a.setInterfaceAveragingWeights(0,wa); a.setInterfaceAveragingWeights(1,wb);
  b.setInterfaceAveragingWeights(1,wb); b.setInterfaceAveragingWeights(0,wa);
  Transfer::Vector actualWeights(3);
  a.diagnosticProlongationWeights(actualWeights);
  for (int i=0;i<3;++i) checkClose(actualWeights[i][0],wa[i][0],"diagnostic averaging weights");
  Transfer::Vector ua(3), ub(3), da(3), db(3);
  ua[0][0]=0.123456789123; ua[1][0]=-0.23194872837; ua[2][0]=1.127384983;
  ub[0][0]=0.78549327812; ub[1][0]=0.31415926535; ub[2][0]=-0.52648391;
  if (zeroProlongate) ua[0][0]=ub[0][0]=0;
  // Primal corner values already agree after the coarse solve.
  for (int i : exact) ub[i][0] = ua[i][0];
  a.preProlongate(ua); b.preProlongate(ub);
  if (encoded)
  {
    std::vector<std::uint32_t> symbols;
    a.appendProlongateSymbols(symbols); b.appendProlongateSymbols(symbols);
    Transfer::buildSharedProlongateCodebookFromSymbols(symbols);
    a.finalizeProlongateEncoding(); b.finalizeProlongateEncoding();
  }
  Transfer::TransmissionType ab,ba;
  a.getInterfaceValues(1,ab); b.getInterfaceValues(0,ba);
  if (!prolongateGroupSize && !encoded && exactDofs==0 && (ab.size()!=4 || ba.size()!=4))
    throw std::runtime_error("Prolongation metadata/payload grew beyond one header and three symbols");
  b.setInterfaceValues(0,ab); a.setInterfaceValues(1,ba);
  a.postProlongate(ua,da); b.postProlongate(ub,db);
  if (!zeroProlongate)
  {
    auto const sent=a.prolongateStats().total_sent_bytes;
    auto const received=a.prolongateStats().total_received_bytes;
    int const precision=groupSize ? 16 : bits;
    auto const pa=a.diagnosticGroupedProlongation(ua.size(),precision,prolongateGroupSize);
    auto const pb=b.diagnosticGroupedProlongation(ub.size(),precision,prolongateGroupSize);
    if (!pa.codecValid || pa.totalBytes!=ab.size()*sizeof(Transfer::TEntry) ||
        pb.totalBytes!=ba.size()*sizeof(Transfer::TEntry))
      throw std::runtime_error("Whole-buffer probe byte model differs from production packets");
    for (int i=0;i<3;++i)
      checkClose(pa.contribution[i][0]+pb.contribution[i][0],ua[i][0]+da[i][0],
                 "whole-buffer producer probe");
    auto const ga=a.diagnosticGroupedProlongation(ua.size(),8,1);
    auto const gb=b.diagnosticGroupedProlongation(ub.size(),8,1);
    for (int i=0;i<3;++i)
    {
      bool const isExact=std::find(exact.begin(),exact.end(),i)!=exact.end();
      double va=wa[i][0]*ua[i][0], vb=wb[i][0]*ub[i][0];
      if (!isExact) { va=static_cast<float>(va); vb=static_cast<float>(vb); }
      checkClose(ga.contribution[i][0]+gb.contribution[i][0],va+vb,"singleton float32 grids");
    }
    if (sent!=a.prolongateStats().total_sent_bytes ||
        received!=a.prolongateStats().total_received_bytes)
      throw std::runtime_error("Grouped probe changed production counters");
    Transfer::Vector after(3);
    a.postProlongate(ua,after);
    for (int i=0;i<3;++i) checkClose(after[i][0],da[i][0],"probe preserves transfer state");
  }
  for (int i=0;i<3;++i) reconstruction.push_back(ua[i][0]+da[i][0]);
  for (int i=0;i<3;++i) checkClose(ua[i][0]+da[i][0],ub[i][0]+db[i][0],"prolongation");
  if (zeroProlongate && (ua[0][0]+da[0][0]!=0 || ub[0][0]+db[0][0]!=0))
    throw std::runtime_error("Zero correction changed during prolongation");
  if (zeroProlongate)
  {
    double maxA=0, maxB=0;
    for (int i=0; i<3; ++i)
      if (std::find(exact.begin(),exact.end(),i)==exact.end())
      {
        maxA=std::max(maxA,std::abs(wa[i][0]*ua[i][0]));
        maxB=std::max(maxB,std::abs(wb[i][0]*ub[i][0]));
      }
    using Q=Kaskade::BDDC::ZeroPreservingQuantization;
    double const maxCode=std::ldexp(1.0,bits)-1;
    double const bound=0.5/Q::scaleForRange(maxA,maxCode)+
                       0.5/Q::scaleForRange(maxB,maxCode);
    for (int i=0; i<3; ++i)
      if (std::abs(ua[i][0]+da[i][0]-(wa[i][0]*ua[i][0]+wb[i][0]*ub[i][0]))>
          bound+64*std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Reconstructed correction exceeds sum of producer quantization bounds");
  }

  for (int i : exact)
    checkClose(ua[i][0]+da[i][0],ua[i][0],"unquantized primal value");

  size_t const prolongationWireBytes = ab.size()*sizeof(Transfer::TEntry);
  Transfer::Vector ra(3),rb(3),sa(3),sb(3);
  std::size_t restrictionWireBytes=0;
  sa=0; sb=0;
  for (int step=0;step<3;++step)
  {
    for (int i=0;i<3;++i)
    {
      ra[i][0]=(0.123456789123-0.3*i)/(step+1);
      rb[i][0]=(0.78549327812+0.17*i)/(step+1);
    }
    if (zeroPreserving) ra[0][0]=rb[0][0]=0;
    a.preRestrict(ra); b.preRestrict(rb);
    if (encoded)
    {
      std::vector<std::uint32_t> symbols;
      a.appendRestrictSymbols(symbols); b.appendRestrictSymbols(symbols);
      Transfer::buildSharedRestrictCodebookFromSymbols(symbols);
      a.finalizeRestrictEncoding(); b.finalizeRestrictEncoding();
    }
    a.getInterfaceResidualUpdate(1,ab); b.getInterfaceResidualUpdate(0,ba);
    restrictionWireBytes += ab.size()*sizeof(Transfer::TEntry);
    if (!encoded && !groupSize && exactDofs==0 && (ab.size()!=4 || ba.size()!=4))
      throw std::runtime_error("Restriction metadata/payload grew beyond one header and three symbols");
    b.setInterfaceResidualUpdate(0,ab); a.setInterfaceResidualUpdate(1,ba);
    a.postRestrict(sa); b.postRestrict(sb);
    for (int i=0;i<3;++i) reconstruction.push_back(sa[i][0]);
    if (zeroPreserving && (sa[0][0]!=0 || sb[0][0]!=0))
      throw std::runtime_error("Zero restriction entry changed during transfer");
    for (int i=0;i<3;++i) checkClose(sa[i][0]/wa[i][0],sb[i][0]/wb[i][0],"restriction");
    // Even though corner corrections agree, their local residual contributions
    // do not vanish. The coarse RHS needs their complete, unquantized sum.
    for (int i : exact)
      checkClose(sa[i][0]/wa[i][0],ra[i][0]+rb[i][0],"corner residual sum");
  }
  if (a.prolongateStats().total_original_sent_bytes != 3*sizeof(double) ||
      a.prolongateStats().total_sent_bytes != prolongationWireBytes ||
      b.prolongateStats().total_received_bytes != prolongationWireBytes)
    throw std::runtime_error("Prolongation byte counters must include unquantized corners");
  if (a.restrictStats().total_original_sent_bytes != 9*sizeof(double))
    throw std::runtime_error("Restriction byte counters must include unquantized corners");
  if (a.restrictStats().total_sent_bytes!=restrictionWireBytes ||
      b.restrictStats().total_received_bytes!=restrictionWireBytes)
    throw std::runtime_error("Restriction byte counters missed wire metadata");
  return reconstruction;
}

std::vector<double> checkThreeOwners(bool encoded, int prolongateGroupSize=0)
{
  // Physical DOFs: a=(0,1,2,3), b=(0,1,2), c=(1,2,3).
  // Shared DOFs 1 and 2 appear in different producer group positions.
  Transfer a(0,{{1,{0,1,2}},{2,{1,2,3}}});
  Transfer b(1,{{0,{0,1,2}},{2,{1,2}}});
  Transfer c(2,{{0,{0,1,2}},{1,{0,1}}});
  std::vector<Transfer*> ts{&a,&b,&c};
  for (int owner=0; owner<3; ++owner)
  {
    auto t=ts[owner];
    t->setProlongateGroupSize(prolongateGroupSize);
    t->setGlobalIndex_active();
    t->setQuantizationBits(12);
    t->setProlongateQuantizationBits(16);
    t->setRestrictZeroPreserving(true);
    t->setRestrictGroupSize(2);
    t->setRestrictEncoding(encoded);
    t->requireRestrictHuffman(encoded);
    t->setRestrictBitlengthEncoding(encoded);
    t->setRestrictZigzagEncoding(encoded);
    t->enableTransform(encoded ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                               : Kaskade::BDDC::TransformType::NONE);
    t->setRestrictTransformAfterQuantization(encoded);
    Transfer::XVector weight(owner==0 ? 4 : 3);
    weight=1;
    t->setInterfaceAveragingWeights(owner,weight);
  }
  for (int from=0; from<3; ++from)
    for (int to=0; to<3; ++to)
      if (from!=to)
      {
        Transfer::XVector w;
        ts[from]->getInterfaceAveragingWeights(to,w);
        ts[to]->setInterfaceAveragingWeights(from,w);
      }
  std::vector<Transfer::Vector> r{Transfer::Vector(4),Transfer::Vector(3),Transfer::Vector(3)};
  auto sum=r;
  for (auto& v:sum) v=0;
  std::vector<double> result;
  for (int step=0; step<3; ++step)
  {
    for (int owner=0; owner<3; ++owner)
    {
      for (int j=0; j<r[owner].size(); ++j)
        r[owner][j][0]=(j==0 ? 10.123456*(owner+1) : -0.000123456*(j+owner))/(step+1);
      ts[owner]->preRestrict(r[owner]);
    }
    if (encoded)
    {
      std::vector<std::uint32_t> symbols;
      for (auto t:ts) t->appendRestrictSymbols(symbols);
      Transfer::buildSharedRestrictCodebookFromSymbols(symbols);
      for (auto t:ts) t->finalizeRestrictEncoding();
    }
    for (int from=0; from<3; ++from)
      for (int to=0; to<3; ++to)
        if (from!=to)
        {
          Transfer::TransmissionType wire;
          ts[from]->getInterfaceResidualUpdate(to,wire);
          ts[to]->setInterfaceResidualUpdate(from,wire);
        }
    for (int owner=0; owner<3; ++owner) ts[owner]->postRestrict(sum[owner]);
    checkClose(sum[0][0][0],sum[1][0][0],"two-owner group");
    checkClose(sum[0][1][0],sum[1][1][0],"three-owner group a/b");
    checkClose(sum[0][1][0],sum[2][0][0],"three-owner group a/c");
    checkClose(sum[0][2][0],sum[1][2][0],"three-owner group a/b second");
    checkClose(sum[0][2][0],sum[2][1][0],"three-owner group a/c second");
    checkClose(sum[0][3][0],sum[2][2][0],"two-owner group last");
    for (auto const& v:sum)
      for (int j=0; j<v.size(); ++j) result.push_back(v[j][0]);
  }
  // Prolongation probe on overlapping packets with non-identical producer unions.
  for (int owner=0; owner<3; ++owner)
  {
    auto t=ts[owner];
    t->setProlongateEncoding(encoded);
    t->setProlongateBitlengthEncoding(encoded);
    t->setProlongateZigzagEncoding(encoded);
    t->setProlongateTransformAfterQuantization(encoded);
    t->preProlongate(r[owner]);
  }
  if (encoded)
  {
    std::vector<std::uint32_t> symbols;
    for (auto t:ts) t->appendProlongateSymbols(symbols);
    Transfer::buildSharedProlongateCodebookFromSymbols(symbols);
    for (auto t:ts) t->finalizeProlongateEncoding();
  }
  std::vector<size_t> packetBytes(3,0);
  for (int from=0; from<3; ++from)
    for (int to=0; to<3; ++to)
      if (from!=to)
      {
        Transfer::TransmissionType wire;
        ts[from]->getInterfaceValues(to,wire);
        packetBytes[from]+=wire.size()*sizeof(Transfer::TEntry);
        ts[to]->setInterfaceValues(from,wire);
      }
  std::vector<std::vector<int>> global{{0,1,2,3},{0,1,2},{1,2,3}};
  std::vector<double> averaged(4,0);
  for (int owner=0; owner<3; ++owner)
  {
    auto probe=ts[owner]->diagnosticGroupedProlongation(r[owner].size(),16,prolongateGroupSize);
    if (probe.totalBytes!=packetBytes[owner])
      throw std::runtime_error("Three-owner packet byte model mismatch");
    for (size_t i=0;i<global[owner].size();++i)
      averaged[global[owner][i]]+=probe.contribution[i][0];
    // Candidate buffers must also be independent of the entropy coding mode.
    for (size_t group: {1,2,128})
    {
      auto candidate=ts[owner]->diagnosticGroupedProlongation(r[owner].size(),8,group);
      for (auto const& v:candidate.contribution) result.push_back(v[0]);
    }
  }
  for (int owner=0;owner<3;++owner)
  {
    Transfer::Vector delta(r[owner].size());
    ts[owner]->postProlongate(r[owner],delta);
    for (size_t i=0;i<global[owner].size();++i)
      checkClose(r[owner][i][0]+delta[i][0],averaged[global[owner][i]],"three-owner P probe");
  }
  return result;
}

std::vector<double> checkLargeProlongation(bool encoded)
{
  std::vector<int> indices(99); std::iota(indices.begin(),indices.end(),0);
  Transfer a(0,{{1,indices}}), b(1,{{0,indices}});
  Transfer::Vector ua(99),ub(99),da(99),db(99);
  for (int i=0;i<99;++i)
  {
    // Different dynamic ranges in each group, including a final singleton.
    ua[i][0]=std::pow(10.0,i/32-2)*(0.3+std::sin(i));
    ub[i][0]=std::pow(10.0,i/32-2)*(0.2-std::cos(i));
  }
  for (int i: {3,34}) ub[i][0]=ua[i][0];
  for (auto t: {&a,&b})
  {
    t->setGlobalIndex_active(); t->setQuantizationBits(8);
    t->setExcludedLocalDofs({3,34}); t->setProlongateGroupSize(32);
    t->setProlongateEncoding(encoded); t->requireProlongateHuffman(encoded);
    t->setProlongateBitlengthEncoding(encoded); t->setProlongateZigzagEncoding(encoded);
    t->enableTransform(encoded ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                               : Kaskade::BDDC::TransformType::NONE);
    t->setProlongateTransformAfterQuantization(encoded);
    Transfer::XVector wa(99),wb(99);
    for (int i=0;i<99;++i) {wa[i][0]=0.1+0.008*i; wb[i][0]=1-wa[i][0];}
    t->setInterfaceAveragingWeights(0,wa); t->setInterfaceAveragingWeights(1,wb);
  }
  a.preProlongate(ua); b.preProlongate(ub);
  if (encoded)
  {
    std::vector<std::uint32_t> symbols;
    a.appendProlongateSymbols(symbols); b.appendProlongateSymbols(symbols);
    Transfer::buildSharedProlongateCodebookFromSymbols(symbols);
    a.finalizeProlongateEncoding(); b.finalizeProlongateEncoding();
  }
  Transfer::TransmissionType ab,ba;
  a.getInterfaceValues(1,ab); b.getInterfaceValues(0,ba);
  b.setInterfaceValues(0,ab); a.setInterfaceValues(1,ba);
  a.postProlongate(ua,da); b.postProlongate(ub,db);
  auto pa=a.diagnosticGroupedProlongation(99,8,32);
  auto pb=b.diagnosticGroupedProlongation(99,8,32);
  if (pa.totalBytes!=ab.size()*sizeof(Transfer::TEntry) ||
      pb.totalBytes!=ba.size()*sizeof(Transfer::TEntry) ||
      a.prolongateStats().total_sent_bytes!=pa.totalBytes ||
      b.prolongateStats().total_received_bytes!=pa.totalBytes ||
      b.prolongateStats().total_sent_bytes!=pb.totalBytes ||
      a.prolongateStats().total_received_bytes!=pb.totalBytes)
    throw std::runtime_error("Multi-group packet counters differ from actual wire bytes");
  std::vector<double> result;
  for (int i=0;i<99;++i)
  {
    checkClose(ua[i][0]+da[i][0],ub[i][0]+db[i][0],"multi-group shared copies");
    checkClose(ua[i][0]+da[i][0],pa.contribution[i][0]+pb.contribution[i][0],"multi-group probe");
    result.push_back(ua[i][0]+da[i][0]);
  }
  return result;
}

void checkProlongationHeaders()
{
  using H=Kaskade::BDDC::ProlongationGroups;
  auto packet=H::pack({0,0,2},{3.1,3.1,5.2},{-0.1,-0.1,0.2},8);
  auto [scales,minima,size]=H::unpack(packet.data(),packet.size(),3,8);
  if (size!=32 || scales[0]!=double(float(3.1)) || minima[2]!=double(float(0.2)))
    throw std::runtime_error("Grouped metadata is not float32 with counted padding");
  auto reject=[&](auto bytes, size_t expected) {
    bool rejected=false;
    try {H::unpack(bytes.data(),bytes.size(),expected,8);}
    catch (std::runtime_error const&) {rejected=true;}
    if (!rejected) throw std::runtime_error("Malformed grouped header accepted");
  };
  reject(packet,2);
  auto truncated=packet; truncated.resize(27); reject(truncated,3);
  auto invalid=packet; float nan=std::numeric_limits<float>::quiet_NaN();
  std::memcpy(invalid.data()+8,&nan,4); reject(invalid,3);
  invalid=packet; std::uint32_t count=0;
  std::memcpy(invalid.data()+4,&count,4); reject(invalid,3);
  auto empty=H::pack({},{},{},8);
  if (std::get<2>(H::unpack(empty.data(),empty.size(),0,8))!=8)
    throw std::runtime_error("Empty grouped header differs from wire minimum");
}

int main()
{
  checkProlongationHeaders();
  if (checkLargeProlongation(false)!=checkLargeProlongation(true))
    throw std::runtime_error("Multi-group plain/full reconstruction differs");
  for (int group: {1,2,32})
  {
    for (int exact: {0,1,3})
      if (checkTransfer(8,false,exact,true,0,false,group)!=
          checkTransfer(8,true,exact,true,0,false,group))
        throw std::runtime_error("Grouped prolongation plain/full reconstruction differs");
    if (checkThreeOwners(false,group)!=checkThreeOwners(true,group))
      throw std::runtime_error("Three-owner grouped prolongation differs");
  }
  int failures = 0;
  for (int bits : {8,12,16,24,32})
    for (bool zeroPreserving : {false,true})
      for (int exactDofs : {0,1,3})
      {
        try
        {
          auto plain = checkTransfer(bits,false,exactDofs,zeroPreserving);
          auto encoded = checkTransfer(bits,true,exactDofs,zeroPreserving);
          if (plain!=encoded)
            throw std::runtime_error("Integer lifting/Huffman changed reconstructed values");
        }
        catch (std::exception const& error)
        {
          ++failures;
          std::cerr << "bits=" << bits << " exact_dofs=" << exactDofs
                    << " zero_preserving=" << zeroPreserving << ": " << error.what() << "\n";
        }
      }
  if (failures != 0) return 1;
  for (int bits : {8,12,16,24,32})
    for (int exact : {0,1,3})
      for (bool zeroR : {false,true})
        if (checkTransfer(bits,false,exact,zeroR,0,true)!=
            checkTransfer(bits,true,exact,zeroR,0,true))
          throw std::runtime_error("Exact-zero prolongation differs between plain and full codec");
  for (int group : {1,2,128})
    for (int exact : {0,1,3})
      if (checkTransfer(12,false,exact,true,group)!=checkTransfer(12,true,exact,true,group))
        throw std::runtime_error("Grouped quantization/full reconstruction differs");
  if (checkThreeOwners(false)!=checkThreeOwners(true))
    throw std::runtime_error("Three-owner grouped quantization/full differs");
  std::cout << "Compressed transfer consistency passed (8/12/16/24/32 bits, "
               "legacy/exact-zero restriction and prolongation; quantization equals integer lifting/Huffman; "
               "mixed/all/no exact corners; R/P groups including multi-group P32; three-owner consistency; wire byte counts).\n";
}
