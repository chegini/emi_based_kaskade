#ifndef EMI_BDDC_PRECISION_HH
#define EMI_BDDC_PRECISION_HH

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
