#include "../EMI_bddc_precision.hh"
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>
struct Options
{
  int bddcCgConstantSplitExperiment=1, bddc=1, sdc=0, direct=0, bddcUseCg=1;
  int bddcFgmres=0, bddcFlexibleCg=0, bddcPreconditionerProbe=0, bddcCgOriginalResidual=1;
  int bddcOperatorAudit=0, bddcConditionEstimate=0, bddcConstantCorrection=0, masterWorker=0;
  int bddcProlongateConstantSplit=1, bddcCompression=1, maximumNumberOfTimeSteps=1;
  int bddcRestrictionPasses=1, bddcRestrictZeroPreserving=0, bddcProlongateZeroPreserving=0;
  int bddcRestrictGroupSize=0, bddcProlongateGroupSize=0, bddcHistogram=0;
  int bddcCgDiagnostics=0, bddcProlongationDiagnostics=0;
  int bddcCompressionBits=16, bddcRestrictCompressionBits=16, bddcProlongateCompressionBits=16;
  int bddcGraphLifting=1, bddcGraphLiftingAfterQuantization=1;
  double bddcRelativeTolerance=1e-7;
};
bool accepts(Options const& o)
{
  try {EmiBddc::validateCgConstantSplitExperiment(o);return true;}
  catch(std::runtime_error const&) {return false;}
}
int main()
{
  for(int bits:{4,6,8,12,16,20,24,32}) for(int flexible:{0,1})
  {
    Options o;o.bddcFlexibleCg=flexible;
    o.bddcRestrictCompressionBits=o.bddcProlongateCompressionBits=bits;
    assert(accepts(o));
    o.bddcRestrictCompressionBits=o.bddcProlongateCompressionBits=-1;
    o.bddcCompressionBits=bits;assert(accepts(o));
  }
  for(int steps:{1,2,10,100}) {Options o;o.bddcFlexibleCg=1;
    o.maximumNumberOfTimeSteps=steps;assert(accepts(o));}
  for(int steps:{0,101}) {Options o;o.bddcFlexibleCg=1;
    o.maximumNumberOfTimeSteps=steps;assert(!accepts(o));}
  std::vector<std::pair<int Options::*,int>> invalid={
    {&Options::bddcCgConstantSplitExperiment,-1},{&Options::bddcCgConstantSplitExperiment,2},
    {&Options::bddc,0},{&Options::sdc,1},{&Options::direct,1},{&Options::bddcUseCg,0},
    {&Options::bddcFgmres,1},{&Options::bddcPreconditionerProbe,1},
    {&Options::bddcCgOriginalResidual,0},{&Options::bddcOperatorAudit,1},
    {&Options::bddcConditionEstimate,1},{&Options::bddcConstantCorrection,1},
    {&Options::masterWorker,1},{&Options::bddcProlongateConstantSplit,0},
    {&Options::bddcCompression,0},{&Options::maximumNumberOfTimeSteps,10},
    {&Options::bddcRestrictionPasses,2},{&Options::bddcRestrictCompressionBits,12},
    {&Options::bddcRestrictZeroPreserving,1},{&Options::bddcProlongateZeroPreserving,1},
    {&Options::bddcRestrictGroupSize,16},{&Options::bddcProlongateGroupSize,16},
    {&Options::bddcHistogram,1},{&Options::bddcCgDiagnostics,1},{&Options::bddcProlongationDiagnostics,1}};
  for(auto const& [field,value]:invalid) {Options o;o.*field=value;assert(!accepts(o));}
  Options o;o.bddcRelativeTolerance=1e-6;assert(!accepts(o));
  o=Options{};o.bddcRestrictCompressionBits=o.bddcProlongateCompressionBits=3;assert(!accepts(o));
  o.bddcCgConstantSplitExperiment=0;assert(accepts(o));
  std::cout << "CG/FCG split opt-in validation passed: 4/6/8/12/16/20/24/32 bits, inherited precision, disabled gate and invalid combinations.\n";
}
