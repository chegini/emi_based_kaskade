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

std::vector<double> checkTransfer(int bits, bool encoded, int exactDofs, bool zeroPreserving, int groupSize=0)
{
  std::vector<double> reconstruction;
  Transfer a(0,{{1,{0,1,2}}}), b(1,{{0,{0,1,2}}});
  for (Transfer* t : {&a,&b})
  {
    t->setGlobalIndex_active();
    t->setQuantizationBits(bits);
    t->setRestrictZeroPreserving(zeroPreserving);
    t->setRestrictGroupSize(groupSize);
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
  Transfer::Vector ua(3), ub(3), da(3), db(3);
  ua[0][0]=0.123456789123; ua[1][0]=-0.23194872837; ua[2][0]=1.127384983;
  ub[0][0]=0.78549327812; ub[1][0]=0.31415926535; ub[2][0]=-0.52648391;
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
  b.setInterfaceValues(0,ab); a.setInterfaceValues(1,ba);
  a.postProlongate(ua,da); b.postProlongate(ub,db);
  for (int i=0;i<3;++i) reconstruction.push_back(ua[i][0]+da[i][0]);
  for (int i=0;i<3;++i) checkClose(ua[i][0]+da[i][0],ub[i][0]+db[i][0],"prolongation");

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

std::vector<double> checkThreeOwners(bool encoded)
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
  return result;
}

int main()
{
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
  for (int group : {1,2,128})
    for (int exact : {0,1,3})
      if (checkTransfer(12,false,exact,true,group)!=checkTransfer(12,true,exact,true,group))
        throw std::runtime_error("Grouped quantization/full reconstruction differs");
  if (checkThreeOwners(false)!=checkThreeOwners(true))
    throw std::runtime_error("Three-owner grouped quantization/full differs");
  std::cout << "Compressed transfer consistency passed (8/12/16/24/32 bits, "
               "legacy/exact-zero restriction; quantization equals integer lifting/Huffman; "
               "mixed/all/no exact corners; groups 1/2/128; three-owner consistency; wire byte counts).\n";
}
