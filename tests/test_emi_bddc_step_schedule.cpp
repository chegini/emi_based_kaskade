#include "../EMI_bddc_step_schedule.hh"
#include <iostream>
#include <stdexcept>

struct Options
{
  int bddcUncompressedSteps=2;
  std::string bddcHistogramSteps="3,5";
  int bddc=1,bddcCompression=1,sdc=0,direct=0,bddcOperatorAudit=0,bddcPreconditionerProbe=0,bddcHistogram=0;
  int bddcRestrictCompressionBits=0,bddcProlongateCompressionBits=0,bddcCompressionBits=16;
  int bddcHuffman=1,bddcBitlength=1,bddcZigzag=1;
  int bddcRestrictZeroPreserving=0,bddcProlongateZeroPreserving=0;
  int bddcRestrictGroupSize=0,bddcProlongateGroupSize=0,bddcCompressionMinInterfaceSize=0;
};
void check(bool yes) { if (!yes) throw std::runtime_error("Schedule assertion failed"); }
template<class F> void rejects(F f)
{
  bool threw=false;
  try { f(); } catch (std::exception const&) { threw=true; }
  check(threw);
}
int main()
{
  using namespace EmiBddc;
  check(parseHistogramSteps("").empty());
  check(parseHistogramSteps("20,6,10")==std::vector<int>({6,10,20}));
  check(parseHistogramSteps("6 10 20")==std::vector<int>({6,10,20}));
  for (auto text: {" ","0","-1","6,6","6,,10","6, ,10",",6","6,","6.0","1e2","2147483648"})
    rejects([&] { parseHistogramSteps(text); });
  for (int total: {0,1,5,10})
    for (int raw: {0,1,5,20})
    {
      int next=0;
      auto phases=bddcStepPhases(total,raw);
      check(phases.size()<=2);
      for (auto const& phase: phases)
      {
        check(phase.first==next && phase.count>0);
        for (int s=phase.first;s<phase.first+phase.count;++s)
          check(phase.compressed==(s>=raw));
        next+=phase.count;
      }
      check(next==total);
    }
  rejects([] { bddcStepPhases(-1,0); });
  rejects([] { bddcStepPhases(1,-1); });
  Options o;
  validateStepSchedule(o,5);
  for (auto member: {&Options::sdc,&Options::direct,&Options::bddcOperatorAudit,
                     &Options::bddcPreconditionerProbe,&Options::bddcHistogram,
                     &Options::bddcRestrictZeroPreserving,&Options::bddcProlongateZeroPreserving,
                     &Options::bddcRestrictGroupSize,&Options::bddcProlongateGroupSize,
                     &Options::bddcCompressionMinInterfaceSize})
  {
    auto bad=o; bad.*member=1;
    rejects([&] { validateStepSchedule(bad,5); });
  }
  for (auto member: {&Options::bddc,&Options::bddcCompression,&Options::bddcHuffman,
                     &Options::bddcBitlength,&Options::bddcZigzag})
  {
    auto bad=o; bad.*member=0;
    rejects([&] { validateStepSchedule(bad,5); });
  }
  for (auto text: {"2","6"})
  { auto bad=o; bad.bddcHistogramSteps=text; rejects([&] { validateStepSchedule(bad,5); }); }
  o.bddcHistogramSteps="";
  o.bddcUncompressedSteps=10;
  validateStepSchedule(o,5); // Entirely raw scheduled run.
  o.bddcUncompressedSteps=0; o.bddcCompression=0; o.direct=1;
  validateStepSchedule(o,5); // Legacy behavior when scheduling is unused.
  std::cout << "Physical-step schedule validation passed\n";
}
