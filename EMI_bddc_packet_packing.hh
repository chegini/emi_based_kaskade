#ifndef EMI_BDDC_PACKET_PACKING_HH
#define EMI_BDDC_PACKET_PACKING_HH
#include <stdexcept>
namespace EmiBddc
{
  // Transport-only option. Keep both affine quantization parameters unchanged.
  // Raw-prefix phases are excluded until their transport counters are accumulated.
  template<class Options>
  void validateLosslessPacketPacking(Options const& o)
  {
    if (o.bddcLosslessPacking<0 || o.bddcLosslessPacking>1)
      throw std::invalid_argument("bddcLosslessPacking must be 0 or 1");
    if (!o.bddcLosslessPacking) return;
    if (!o.bddc || !o.mpi || !o.bddcMpiRouting || !o.bddcCompression ||
        o.sdc || o.direct || o.bddcOperatorAudit || o.bddcPreconditionerProbe ||
        o.bddcProlongateConstantSplit || o.bddcRestrictGroupSize ||
        o.bddcProlongateGroupSize || o.bddcLocalCompressionBypass || o.bddcUncompressedSteps)
      throw std::invalid_argument("Lossless packing requires compressed non-SDC MPI BDDC "
          "with destination routing, no audit/probe, constant splitting, grouping, "
          "local bypass or raw-prefix schedule; quantization parameters are unchanged");
  }

  template<class Solver,class Options>
  void configureLosslessPacketPacking(Solver& solver,Options const& o)
  {
    if constexpr (requires { solver.setMpiLosslessPacking(false); })
      solver.setMpiLosslessPacking(o.bddcLosslessPacking!=0);
    else if (o.bddcLosslessPacking)
      throw std::runtime_error("Kaskade headers lack lossless MPI packet packing");
  }
}
#endif
