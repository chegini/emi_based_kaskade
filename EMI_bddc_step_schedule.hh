#ifndef EMI_BDDC_STEP_SCHEDULE_HH
#define EMI_BDDC_STEP_SCHEDULE_HH

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace EmiBddc
{
  // Physical steps are one-based in the UI; first is a zero-based time offset.
  // A raw phase uses the ordinary SpaceTransfer, not an uncoded quantizer.
  struct BddcStepPhase
  {
    int first, count;
    bool compressed;
  };

  inline std::vector<int> parseHistogramSteps(std::string text)
  {
    if (text.empty()) return {};
    // The application option parser accepts one string. Commas (or quoted
    // spaces) select distinct physical steps, never a range of iterations.
    std::size_t start=0;
    do
    {
      auto end=text.find(',',start);
      if (text.substr(start,end==std::string::npos ? end : end-start)
              .find_first_not_of(" \t\r\n")==std::string::npos)
        throw std::invalid_argument("Empty histogram step");
      if (end==std::string::npos) break;
      start=end+1;
    } while (true);
    std::replace(text.begin(),text.end(),',',' ');
    std::istringstream input(text);
    std::vector<int> steps;
    for (std::string token; input>>token;)
    {
      if (token.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("Histogram steps must be positive integers");
      int step=std::stoi(token);
      if (step<1) throw std::invalid_argument("Histogram steps are one-based");
      steps.push_back(step);
    }
    if (steps.empty()) throw std::invalid_argument("Empty histogram selection");
    std::sort(steps.begin(),steps.end());
    if (std::adjacent_find(steps.begin(),steps.end())!=steps.end())
      throw std::invalid_argument("Duplicate histogram step");
    return steps;
  }

  inline std::vector<BddcStepPhase> bddcStepPhases(int steps, int uncompressed)
  {
    if (steps<0 || uncompressed<0)
      throw std::invalid_argument("Step counts must be nonnegative");
    int raw=std::min(steps,uncompressed);
    std::vector<BddcStepPhase> result;
    if (raw) result.push_back({0,raw,false});
    if (raw<steps) result.push_back({raw,steps-raw,true});
    return result;
  }

  template<class Options>
  void validateStepSchedule(Options const& o, int steps)
  {
    auto selected=parseHistogramSteps(o.bddcHistogramSteps);
    if (o.bddcUncompressedSteps<0)
      throw std::invalid_argument("bddcUncompressedSteps must be nonnegative");
    if (!o.bddcUncompressedSteps && selected.empty()) return;
    if (!o.bddc || !o.bddcCompression || o.sdc || o.direct ||
        o.bddcOperatorAudit || o.bddcPreconditionerProbe || o.bddcHistogram)
      throw std::invalid_argument("Step scheduling requires compressed non-SDC BDDC; no direct/audit/probe or legacy histogram mode");
    // Packet capture observes full affine coding; no hidden change of codec.
    if (!selected.empty())
    {
      int r=o.bddcRestrictCompressionBits>0 ? o.bddcRestrictCompressionBits : o.bddcCompressionBits;
      int p=o.bddcProlongateCompressionBits>0 ? o.bddcProlongateCompressionBits : o.bddcCompressionBits;
      if (!o.bddcHuffman || !o.bddcBitlength || !o.bddcZigzag ||
          o.bddcRestrictZeroPreserving || o.bddcProlongateZeroPreserving ||
          o.bddcRestrictGroupSize || o.bddcProlongateGroupSize ||
          o.bddcCompressionMinInterfaceSize || r<2 || r>24 || p<2 || p>24)
        throw std::invalid_argument("Selected-step packet histograms require affine 2..24-bit Huffman/bit-length/zigzag coding without grouping or thresholds");
      for (int step:selected)
        if (step>steps || step<=o.bddcUncompressedSteps)
          throw std::invalid_argument("Histogram steps must exist and follow the uncompressed starting steps");
    }
  }
}
#endif
