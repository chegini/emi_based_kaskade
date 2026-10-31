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
  struct LegacyOptions { int bddcCompressionBits = 24; } legacy;
  EmiBddc::configureTransferPrecision(transfer,legacy);
  check(transfer,24,24);
  struct UncompressedTransfer {} uncompressed;
  EmiBddc::configureTransferPrecision(uncompressed,options);
  for (int bits : {-1,1,8,12,16,24,32})
    if (!EmiBddc::validDirectionalCompressionBits(bits)) return 1;
  for (int bits : {-2,0,33})
    if (EmiBddc::validDirectionalCompressionBits(bits)) return 2;
  std::cout << "Directional precision checks passed (inheritance, overrides, validation).\n";
}
