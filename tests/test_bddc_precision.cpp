#include <stdexcept>
#include <iostream>

#include "EMI_bddc_precision.hh"

struct Options
{
  int bddcCompressionBits = 16;
  int bddcRestrictCompressionBits = -1;
  int bddcProlongateCompressionBits = -1;
};

struct Transfer
{
  int restrictBits = -1;
  int prolongateBits = -1;
  void setQuantizationBits(int bits) { restrictBits = prolongateBits = bits; }
  void setRestrictQuantizationBits(int bits) { restrictBits = bits; }
  void setProlongateQuantizationBits(int bits) { prolongateBits = bits; }
};

void check(Transfer const& transfer, int restriction, int prolongation)
{
  if (transfer.restrictBits != restriction || transfer.prolongateBits != prolongation)
    throw std::runtime_error("Directional precision forwarding is incorrect");
}

int main()
{
  Options options;
  Transfer transfer;
  EmiBddc::configureTransferPrecision(transfer,options);
  check(transfer,16,16);
  options.bddcCompressionBits = 8;
  options.bddcRestrictCompressionBits = 32;
  EmiBddc::configureTransferPrecision(transfer,options);
  check(transfer,32,8);
  options.bddcRestrictCompressionBits = -1;
  options.bddcProlongateCompressionBits = 32;
  EmiBddc::configureTransferPrecision(transfer,options);
  check(transfer,8,32);
  options.bddcRestrictCompressionBits = 12;
  options.bddcProlongateCompressionBits = 24;
  EmiBddc::configureTransferPrecision(transfer,options);
  check(transfer,12,24);
  options.bddcRestrictCompressionBits = options.bddcProlongateCompressionBits = -1;
  EmiBddc::configureTransferPrecision(transfer,options);
  check(transfer,8,8); // Reusing a transfer clears previous overrides.
  struct ZeroOptions : Options { int bddcRestrictZeroPreserving = 1; } zeroOptions;
  struct ZeroTransfer : Transfer
  {
    bool zero = false;
    void setRestrictZeroPreserving(bool enable) { zero = enable; }
  } zeroTransfer;
  EmiBddc::configureTransferPrecision(zeroTransfer,zeroOptions);
  if (!zeroTransfer.zero) return 3;
  zeroOptions.bddcRestrictZeroPreserving = 0;
  EmiBddc::configureTransferPrecision(zeroTransfer,zeroOptions);
  if (zeroTransfer.zero) return 4;
  zeroOptions.bddcRestrictZeroPreserving = 1;
  bool rejected = false;
  try { EmiBddc::configureTransferPrecision(transfer,zeroOptions); }
  catch (std::runtime_error const&) { rejected = true; }
  if (!rejected) return 5;
  struct GroupOptions : Options { int bddcRestrictGroupSize = 128; } groupOptions;
  struct GroupTransfer : Transfer {
    int group = -1;
    void setRestrictGroupSize(int n) { group=n; }
  } groupTransfer;
  EmiBddc::configureTransferPrecision(groupTransfer,groupOptions);
  if (groupTransfer.group!=128) return 6;
  groupOptions.bddcRestrictGroupSize=0;
  EmiBddc::configureTransferPrecision(groupTransfer,groupOptions);
  if (groupTransfer.group!=0) return 7;
  groupOptions.bddcRestrictGroupSize=128;
  rejected=false;
  try { EmiBddc::configureTransferPrecision(transfer,groupOptions); }
  catch (std::runtime_error const&) { rejected=true; }
  if (!rejected) return 8;
  struct LegacyOptions { int bddcCompressionBits = 24; } legacy;
  struct ProlongateOptions : Options { int bddcProlongateZeroPreserving = 1; } pOptions;
  struct ProlongateTransfer : Transfer {
    bool zero = false;
    void setProlongateZeroPreserving(bool b) { zero=b; }
  } pTransfer;
  EmiBddc::configureTransferPrecision(pTransfer,pOptions);
  if (!pTransfer.zero) return 9;
  pOptions.bddcProlongateZeroPreserving=0;
  EmiBddc::configureTransferPrecision(pTransfer,pOptions);
  if (pTransfer.zero) return 10;
  pOptions.bddcProlongateZeroPreserving=1;
  rejected=false;
  try { EmiBddc::configureTransferPrecision(transfer,pOptions); }
  catch (std::runtime_error const&) { rejected=true; }
  if (!rejected) return 11;
  EmiBddc::configureTransferPrecision(transfer,legacy);
  check(transfer,24,24);
  struct PGroupOptions : Options { int bddcProlongateGroupSize=32; } pgOptions;
  struct PGroupTransfer : Transfer {
    int group=-1;
    void setProlongateGroupSize(int n) { group=n; }
  } pgTransfer;
  EmiBddc::configureTransferPrecision(pgTransfer,pgOptions);
  if (pgTransfer.group!=32) return 12;
  rejected=false;
  try { EmiBddc::configureTransferPrecision(transfer,pgOptions); }
  catch (std::runtime_error const&) { rejected=true; }
  if (!rejected) return 13;
  pgOptions.bddcProlongateGroupSize=0;
  EmiBddc::configureTransferPrecision(pgTransfer,pgOptions);
  if (pgTransfer.group!=0) return 14;
  struct UncompressedTransfer {} uncompressed;
  struct SplitOptions : Options { int bddcProlongateConstantSplit=1; } splitOptions;
  struct SplitTransfer : Transfer {
    bool split=false;
    void setProlongateConstantSplit(bool enable) { split=enable; }
  } splitTransfer;
  EmiBddc::configureTransferPrecision(splitTransfer,splitOptions);
  if (!splitTransfer.split) return 15;
  rejected=false;
  try { EmiBddc::configureTransferPrecision(transfer,splitOptions); }
  catch (std::runtime_error const&) { rejected=true; }
  if (!rejected) return 16;
  splitOptions.bddcProlongateConstantSplit=0;
  EmiBddc::configureTransferPrecision(splitTransfer,splitOptions);
  if (splitTransfer.split) return 17;
  EmiBddc::configureTransferPrecision(uncompressed,options);
  for (int bits : {-1,1,8,12,16,24,32})
    if (!EmiBddc::validDirectionalCompressionBits(bits)) return 1;
  for (int bits : {-2,0,33})
    if (EmiBddc::validDirectionalCompressionBits(bits)) return 2;
  std::cout << "Directional precision checks passed (inheritance, overrides, validation).\n";
}
