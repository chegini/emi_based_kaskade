#ifndef EMI_BDDC_PRECISION_HH
#define EMI_BDDC_PRECISION_HH

#include <stdexcept>

namespace EmiBddc
{
  // EMI uses uint32_t quantized symbols. -1 inherits the common bit setting;
  // explicit directional settings must fit that symbol type.
  inline bool validDirectionalCompressionBits(int bits)
  {
    return bits == -1 || (bits >= 1 && bits <= 32);
  }

  /**
   * Configure optional restriction/prolongation quantization overrides.
   * The default retains the existing common precision in both directions.
   * Mixed precision is a separate numerical configuration used to isolate
   * transfer sensitivity; it must not be labelled as uniform low precision.
   */
  template <class Transfer, class Options>
  void configureTransferPrecision(Transfer& transfer, Options const& options)
  {
    if constexpr (requires { transfer.setQuantizationBits(options.bddcCompressionBits); })
    {
      transfer.setQuantizationBits(options.bddcCompressionBits);
      if constexpr (requires { options.bddcProlongateConstantSplit; })
      {
        if constexpr (requires { transfer.setProlongateConstantSplit(false); })
          transfer.setProlongateConstantSplit(options.bddcProlongateConstantSplit!=0);
        else if (options.bddcProlongateConstantSplit)
          throw std::runtime_error("Kaskade headers lack constant-split prolongation");
      }
      if constexpr (requires { options.bddcProlongateGroupSize; })
      {
        if constexpr (requires { transfer.setProlongateGroupSize(0); })
          transfer.setProlongateGroupSize(options.bddcProlongateGroupSize);
        else if (options.bddcProlongateGroupSize)
          throw std::runtime_error("Kaskade headers lack grouped prolongation");
      }
      if constexpr (requires { options.bddcProlongateZeroPreserving; })
      {
        if constexpr (requires { transfer.setProlongateZeroPreserving(false); })
          transfer.setProlongateZeroPreserving(options.bddcProlongateZeroPreserving!=0);
        else if (options.bddcProlongateZeroPreserving)
          throw std::runtime_error("Kaskade headers lack zero-preserving prolongation");
      }
      if constexpr (requires { options.bddcRestrictGroupSize; })
      {
        if constexpr (requires { transfer.setRestrictGroupSize(0); })
          transfer.setRestrictGroupSize(options.bddcRestrictGroupSize);
        else if (options.bddcRestrictGroupSize)
          throw std::runtime_error("Kaskade headers lack grouped restriction quantization");
      }
      if constexpr (requires { options.bddcRestrictZeroPreserving; })
      {
        if constexpr (requires { transfer.setRestrictZeroPreserving(false); })
          transfer.setRestrictZeroPreserving(options.bddcRestrictZeroPreserving!=0);
        else if (options.bddcRestrictZeroPreserving)
          throw std::runtime_error("Kaskade headers lack the zero-preserving restriction experiment");
      }
      if constexpr (requires {
                      options.bddcRestrictCompressionBits;
                      transfer.setRestrictQuantizationBits(1);
                    })
        if (options.bddcRestrictCompressionBits > 0)
          transfer.setRestrictQuantizationBits(options.bddcRestrictCompressionBits);
      if constexpr (requires {
                      options.bddcProlongateCompressionBits;
                      transfer.setProlongateQuantizationBits(1);
                    })
        if (options.bddcProlongateCompressionBits > 0)
          transfer.setProlongateQuantizationBits(options.bddcProlongateCompressionBits);
    }
  }
}

#endif
