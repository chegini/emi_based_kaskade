#include "EMI_bddc_packet_packing.hh"
#include <cassert>
#include <initializer_list>
struct Options {
  int bddcLosslessPacking=0, bddc=1, mpi=1, bddcMpiRouting=1, bddcCompression=1;
  int sdc=0, direct=0, bddcOperatorAudit=0, bddcPreconditionerProbe=0;
  int bddcProlongateConstantSplit=0, bddcRestrictGroupSize=0, bddcProlongateGroupSize=0;
  int bddcLocalCompressionBypass=0, bddcUncompressedSteps=0;
};
struct NewSolver { bool enabled=false; void setMpiLosslessPacking(bool on) { enabled=on; } };
struct OldSolver {};
int main() {
  Options o;
  EmiBddc::validateLosslessPacketPacking(o);
  o.bddcLosslessPacking=1;
  EmiBddc::validateLosslessPacketPacking(o);
  for (auto field : {&Options::bddc, &Options::mpi, &Options::bddcMpiRouting, &Options::bddcCompression,
      &Options::sdc, &Options::direct, &Options::bddcOperatorAudit, &Options::bddcPreconditionerProbe,
      &Options::bddcRestrictGroupSize,
      &Options::bddcProlongateGroupSize, &Options::bddcLocalCompressionBypass, &Options::bddcUncompressedSteps}) {
    auto bad=o; bad.*field=!(bad.*field);
    bool rejected=false;
    try { EmiBddc::validateLosslessPacketPacking(bad); } catch(std::invalid_argument const&) { rejected=true; }
    assert(rejected);
    bad.bddcLosslessPacking=0;
    EmiBddc::validateLosslessPacketPacking(bad);
  }
  o.bddcProlongateConstantSplit=1;
  EmiBddc::validateLosslessPacketPacking(o);
  NewSolver modern;
  EmiBddc::configureLosslessPacketPacking(modern,o); assert(modern.enabled);
  OldSolver old;
  bool rejected=false;
  try { EmiBddc::configureLosslessPacketPacking(old,o); } catch(std::runtime_error const&) { rejected=true; }
  assert(rejected);
  o.bddcLosslessPacking=0;
  EmiBddc::configureLosslessPacketPacking(old,o);
  EmiBddc::configureLosslessPacketPacking(modern,o); assert(!modern.enabled);
  for (int invalid : {-1,2}) {
    o.bddcLosslessPacking=invalid; rejected=false;
    try { EmiBddc::validateLosslessPacketPacking(o); } catch(std::invalid_argument const&) { rejected=true; }
    assert(rejected);
  }
}
