#include <algorithm>
#include <array>
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

int main()
{
  Transfer a(0,{{1,{0,1,2,3,4}}}), b(1,{{0,{0,1,2,3,4}}});
  for (Transfer* t : {&a,&b})
  {
    t->setGlobalIndex_active();
    t->setQuantizationBits(8);
    t->setRestrictEncoding(false);
    t->setProlongateEncoding(false);
  }
  Transfer::XVector wa(5),wb(5);
  Transfer::Vector ra(5),rb(5),sa(5),sb(5);
  std::array<double,5> const ar{0.3,-0.1,0.51,-0.43,0.2};
  std::array<double,5> const br{-0.17,0.24,-0.11,-0.37,-0.07};
  for (int i=0;i<5;++i)
  {
    wa[i][0]=0.1*(i+1); wb[i][0]=1-wa[i][0];
    ra[i][0]=ar[i]; rb[i][0]=br[i];
  }
  a.setInterfaceAveragingWeights(0,wa); a.setInterfaceAveragingWeights(1,wb);
  b.setInterfaceAveragingWeights(1,wb); b.setInterfaceAveragingWeights(0,wa);
  sa=0; sb=0;
  double firstError=0,finalError=0;
  for (int pass=0;pass<4;++pass)
  {
    // Keep the original local residual fixed while quantization feedback is refined.
    a.preRestrict(ra); b.preRestrict(rb);
    Transfer::TransmissionType ab,ba;
    a.getInterfaceResidualUpdate(1,ab); b.getInterfaceResidualUpdate(0,ba);
    b.setInterfaceResidualUpdate(0,ab); a.setInterfaceResidualUpdate(1,ba);
    a.postRestrict(sa); b.postRestrict(sb);
    double error=0, sum=0;
    for (int i=0;i<5;++i)
    {
      error=std::max(error,std::abs(sa[i][0]/wa[i][0]-(ar[i]+br[i])));
      error=std::max(error,std::abs(sb[i][0]/wb[i][0]-(ar[i]+br[i])));
      sum+=sa[i][0]+sb[i][0];
    }
    std::cout << "8-bit restriction pass " << pass+1 << ": max_error=" << error
              << " compatible_sum_defect=" << std::abs(sum) << '\n';
    if (pass==0) firstError=error;
    finalError=error;
  }
  if (!std::isfinite(finalError) || firstError<=1e-6 || finalError>firstError*1e-6)
    throw std::runtime_error("8-bit feedback refinement did not reduce restriction error");
  if (a.restrictStats().total_original_sent_bytes!=4*5*sizeof(double))
    throw std::runtime_error("Refinement communication was not fully counted");
  std::cout << "Restriction refinement checks passed; all four 8-bit exchanges counted.\n";
}
