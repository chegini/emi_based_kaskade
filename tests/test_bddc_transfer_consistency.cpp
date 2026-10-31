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

void checkTransfer(int bits, bool encoded, int exactDofs)
{
  Transfer a(0,{{1,{0,1,2}}}), b(1,{{0,{0,1,2}}});
  for (Transfer* t : {&a,&b})
  {
    t->setGlobalIndex_active();
    t->setQuantizationBits(bits);
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
  for (int i=0;i<3;++i) checkClose(ua[i][0]+da[i][0],ub[i][0]+db[i][0],"prolongation");

  for (int i : exact)
    checkClose(ua[i][0]+da[i][0],ua[i][0],"unquantized primal value");

  size_t const prolongationWireBytes = ab.size()*sizeof(Transfer::TEntry);
  Transfer::Vector ra(3),rb(3),sa(3),sb(3);
  sa=0; sb=0;
  for (int step=0;step<3;++step)
  {
    for (int i=0;i<3;++i)
    {
      ra[i][0]=(0.123456789123-0.3*i)/(step+1);
      rb[i][0]=(0.78549327812+0.17*i)/(step+1);
    }
    a.preRestrict(ra); b.preRestrict(rb);
    if (encoded)
    {
      std::vector<std::uint32_t> symbols;
      a.appendRestrictSymbols(symbols); b.appendRestrictSymbols(symbols);
      Transfer::buildSharedRestrictCodebookFromSymbols(symbols);
      a.finalizeRestrictEncoding(); b.finalizeRestrictEncoding();
    }
    a.getInterfaceResidualUpdate(1,ab); b.getInterfaceResidualUpdate(0,ba);
    b.setInterfaceResidualUpdate(0,ab); a.setInterfaceResidualUpdate(1,ba);
    a.postRestrict(sa); b.postRestrict(sb);
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
}

int main()
{
  int failures = 0;
  for (int bits : {12,16,24,32})
    for (bool encoded : {false,true})
      for (int exactDofs : {0,1,3})
    {
      try
      {
        checkTransfer(bits,encoded,exactDofs);
      }
      catch (std::exception const& error)
      {
        ++failures;
        std::cerr << "bits=" << bits << " encoded=" << encoded << " exact_dofs=" << exactDofs << ": "
                  << error.what() << "\n";
      }
    }
  if (failures != 0)
    return 1;
  std::cout << "Compressed transfer consistency passed (12/16/24/32 bits, "
               "quantization and integer lifting/Huffman, mixed/all/no exact corners).\n";
}
